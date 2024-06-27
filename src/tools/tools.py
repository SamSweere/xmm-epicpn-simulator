import os
import shutil
import tarfile
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from functools import partial
from operator import countOf
from pathlib import Path
from tempfile import TemporaryDirectory, mkdtemp

import numpy as np
import requests
from loguru import logger
from tqdm import tqdm

from src.config import DownloadCfg, EnergyCfg, EnvironmentCfg, MultiprocessingCfg, SimputCfg
from src.illustris_tng.fits import cutout_to_xray_fits
from src.illustris_tng.web_api import (
    get_available_simulations,
    get_cutouts,
    get_subhalos,
)
from src.simput.tools import get_spectrumfile
from src.tools.files import compress_gzip, decompress_targz
from src.xmm.tools import get_spectrum_file


@contextmanager
def _tmp_chdir(path: Path):
    old_dir = os.getcwd()
    os.chdir(path)

    try:
        yield
    finally:
        os.chdir(old_dir)


def download_cloudy_emissivity(environment: EnvironmentCfg):
    cloudy_emissivity = environment.working_dir / "cloudy_emissivity_v2.h5"
    if not cloudy_emissivity.exists():
        logger.info(f"Downloading cloudy_emissivity_v2.h5 to {cloudy_emissivity.resolve()}")
        retries = 3
        while retries > 0:
            try:
                with requests.get(
                    "http://yt-project.org/data/cloudy_emissivity_v2.h5",
                    stream=True,
                ) as r:
                    r.raise_for_status()
                    with open(cloudy_emissivity, "wb") as f:
                        for chunk in r.iter_content(chunk_size=int(1e6)):
                            f.write(chunk)
                retries = 0
            except:  # noqa
                retries = retries - 1

        if not cloudy_emissivity.exists():
            raise FileNotFoundError(f"Failed to load cloudy_emissivity_v2.h5 {cloudy_emissivity}!")
    return cloudy_emissivity


def download_data(
    download_cfg: DownloadCfg,
    energies: EnergyCfg,
    env_cfg: EnvironmentCfg,
    mp_cfg: MultiprocessingCfg,
    api_key: str,
    delete_product: bool,
) -> None:
    decompress_fs: dict[str, Future] = {}
    _get_cutouts = partial(
        get_cutouts,
        api_key=api_key,
        cutout_datafolder=download_cfg.cutouts_path,
        fail_on_error=env_cfg.fail_on_error,
    )
    _get_subhalos = partial(
        get_subhalos,
        api_key=api_key,
        params={
            "limit": download_cfg.top_n,
            "primary_flag": 1,
            "order_by": "-mass_gas",
        },
    )
    cloudy_emissivity = download_cloudy_emissivity(env_cfg)
    _cutout_to_xray_fits = partial(
        cutout_to_xray_fits,
        output_dir=download_cfg.fits_path,
        mode_dict=download_cfg.modes,
        cloudy_emissivity_root=cloudy_emissivity.parent,
        energies=energies,
        resolutions=download_cfg.resolutions,
        environment=env_cfg,
    )

    with ProcessPoolExecutor(max_workers=mp_cfg.num_cores, max_tasks_per_child=100) as executor:
        if download_cfg.cutouts_compressed.exists():
            logger.info(f"Found compressed cutouts in {download_cfg.cutouts_compressed}. Decompressing...")
            decompress_fs["cutouts"] = executor.submit(
                decompress_targz,
                in_file_path=download_cfg.cutouts_compressed,
                out_file_dir=download_cfg.cutouts_path,
            )
        if download_cfg.fits_compressed.exists():
            logger.info(f"Found compressed FITS in {download_cfg.fits_compressed.resolve()}. Decompressing...")
            decompress_fs["fits"] = executor.submit(
                decompress_targz,
                in_file_path=download_cfg.fits_compressed,
                out_file_dir=download_cfg.fits_path,
            )

        logger.info("START\tGetting simulations.")
        simulations: list[str] = []
        subhalos_fs = []
        for simulation_name, simulation_url in get_available_simulations(api_key=api_key):
            if simulation_name in download_cfg.simulations:
                simulations.append(simulation_name)
                for snapshot in download_cfg.snapshots:
                    subhalos_fs.append(
                        executor.submit(_get_subhalos, simulation_url=simulation_url, snapshot_num=snapshot)
                    )

        if not simulations:
            raise ValueError("No simulations found! Please check your config file.")

        logger.success(f"DONE\tFound {len(simulations)} available simulations. ({', '.join(simulations)})")

        # Wait for the decompression of the cutouts to finish if it has been started
        _ = decompress_fs.pop("cutouts").result() if "cutouts" in decompress_fs else None

        logger.info("START\tGetting subhalos.")
        cutouts_fs = []
        with tqdm(total=len(subhalos_fs), desc="Getting subhalos") as pbar:
            for future in as_completed(subhalos_fs):
                for subhalo_url in future.result():
                    logger.success(f"Got subhalo {subhalo_url} -> Start getting cutouts.")
                    cutouts_fs.append(
                        executor.submit(
                            _get_cutouts,
                            subhalo_url=subhalo_url,
                        )
                    )
                pbar.update()
            elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
        if not cutouts_fs:
            raise ValueError("No subhalos found! Please check your config file.")

        logger.success(f"DONE\tGot {len(cutouts_fs)} subhalos. Duration: {elapsed_time}")

        # Wait for the decompression of the FITS to finish if it has been started
        _ = decompress_fs.pop("fits").result() if "fits" in decompress_fs else None

        logger.info("START\tGetting cutouts.")
        tars = {}
        if env_cfg.tar_and_compress:
            tars["cutouts"] = (tarfile.open(download_cfg.cutouts_tar, "a"), download_cfg.cutouts_tar)
            tars["fits"] = (tarfile.open(download_cfg.fits_tar, "a"), download_cfg.fits_tar)
        cutouts_fits_fs = {}
        with tqdm(total=len(cutouts_fs), desc="Getting cutouts") as pbar:
            for future in as_completed(cutouts_fs):
                cutout = future.result()
                path: Path = cutout["file"]
                logger.success(f"Got cutout {path} -> Start generation of X-ray FITS.")
                fs = executor.submit(
                    _cutout_to_xray_fits,
                    cutout=path,
                    sub={"pos_x": float(cutout["x"]), "pos_y": float(cutout["y"]), "pos_z": float(cutout["z"])},
                    widths=download_cfg.simulations[cutout["file"].parts[-3]],
                    redshift=download_cfg.snapshots[int(cutout["file"].parts[-2])],
                )
                cutouts_fits_fs[fs] = path
                if "cutouts" in tars:
                    tar, tar_path = tars["cutouts"]
                    tar.add(path, path.relative_to(download_cfg.cutouts_path))
                    logger.success(f"Added {path} to {tar_path}.")
                pbar.update()
            elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])

        if "cutouts" in tars:
            tar, tar_path = tars["cutouts"]
            tar.close()
            executor.submit(
                compress_gzip,
                in_file_path=tar_path,
                out_file_path=download_cfg.cutouts_compressed,
                remove_file=True,
            )

        logger.success(f"DONE\tGetting cutouts. Duration: {elapsed_time}")

        logger.info("START\tGenerating FITS from cutouts.")
        with tqdm(total=len(cutouts_fits_fs), desc="Creating FITS from cutouts") as pbar:
            for future in as_completed(cutouts_fits_fs):
                cutout = cutouts_fits_fs[future]
                logger.success(f"Converted {cutout} to X-ray FITS.")
                if "fits" in tars:
                    tar, tar_path = tars["fits"]
                    for fits in future.result():
                        tar.add(fits, fits.relative_to(download_cfg.fits_path))
                        logger.success(f"Added {fits} to {tar_path}.")
                        if delete_product:
                            fits.unlink()
                            logger.success(f"Deleted {fits}.")
                pbar.update()
            elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
        logger.success(f"DONE\tCreated FITS from cutouts. Duration: {elapsed_time}")

        if "fits" in tars:
            tar, tar_path = tars["fits"]
            tar.close()
            executor.submit(
                compress_gzip,
                in_file_path=tar_path,
                out_file_path=download_cfg.fits_compressed,
                remove_file=True,
            )


def generate_simput(
    simput_cfg: SimputCfg,
    energies: EnergyCfg,
    env_cfg: EnvironmentCfg,
    mp_cfg: MultiprocessingCfg,
    satellites: list,
    agn_counts_file: Path | None,
    delete_product: bool,
) -> None:
    generator_fs = {}
    with TemporaryDirectory(prefix="simput_") as tmp_dir:
        tmp_dir = Path(tmp_dir)

        with ProcessPoolExecutor(max_workers=mp_cfg.num_cores, max_tasks_per_child=100) as executor:
            if simput_cfg.img.n_gen != 0:
                from src.simput.image import simput_image

                logger.info("START\tGenerating SIMPUT for mode 'img'...")
                img_path = simput_cfg.simput_dir / "img"
                img_path.mkdir(parents=True, exist_ok=True)
                rng = np.random.default_rng()
                xspec_file = get_spectrumfile(run_dir=tmp_dir)

                _simput_image = partial(
                    simput_image,
                    emin=energies.emin,
                    emax=energies.emax,
                    run_dir=Path(mkdtemp(dir=tmp_dir, prefix="img_")),
                    xspec_file=xspec_file,
                    consume_data=env_cfg.consume_data,
                )

                if simput_cfg.fits_compressed.exists():
                    logger.info(f"Found compressed FITS files in {simput_cfg.fits_compressed}. Decompressing...")
                    decompress_targz(
                        in_file_path=simput_cfg.fits_compressed,
                        out_file_dir=simput_cfg.fits_dir,
                    )

                amount_img = simput_cfg.img.n_gen
                for in_file in simput_cfg.fits_dir.rglob("*.fits"):
                    if amount_img > 0:
                        amount_img = amount_img - 1

                    tng_set, snapshot_num = in_file.parts[-3], in_file.parts[-2]
                    # Check how many files have already been generated and how many are left to generate
                    sample_num = simput_cfg.num_img_sample
                    for _ in (img_path / tng_set / snapshot_num).glob(f"{in_file.stem}*"):
                        sample_num = sample_num - 1
                        if sample_num == 0:
                            break

                    # None are left to be generated => Skip
                    if sample_num == 0:
                        logger.debug(f"Won't generate any images for {in_file.name}.")
                        if env_cfg.consume_data:
                            in_file.unlink()
                        continue

                    logger.info(f"Will generate {sample_num} images for {in_file.name}.")

                    zoom = np.round(
                        rng.uniform(
                            simput_cfg.zoom_range[0],
                            simput_cfg.zoom_range[1],
                            sample_num,
                        ),
                        2,
                    )
                    sigma_b = np.round(
                        rng.uniform(
                            simput_cfg.sigma_b_range[0],
                            simput_cfg.sigma_b_range[1],
                            sample_num,
                        ),
                        2,
                    )

                    offset_x = np.round(
                        rng.normal(
                            -simput_cfg.offset_std,
                            simput_cfg.offset_std,
                            sample_num,
                        ),
                        2,
                    )
                    offset_y = np.round(
                        rng.normal(
                            -simput_cfg.offset_std,
                            simput_cfg.offset_std,
                            sample_num,
                        ),
                        2,
                    )

                    fs = executor.submit(
                        _simput_image,
                        img_path_in=in_file,
                        zooms=zoom,
                        sigmas_b=sigma_b,
                        offsets_x=offset_x,
                        offsets_y=offset_y,
                        output_dir=img_path / tng_set / snapshot_num,
                    )

                    generator_fs[fs] = "img"

                    if amount_img == 0:
                        break

            if simput_cfg.bkg.n_gen:
                from src.simput.background import create_background

                logger.info("START\tGenerating SIMPUT for mode BKG...")
                bkg_path = simput_cfg.simput_dir / "bkg"
                bkg_path.mkdir(parents=True, exist_ok=True)

                _background = partial(
                    create_background,
                    run_dir=Path(mkdtemp(dir=tmp_dir, prefix="bkg_")),
                    output_dir=bkg_path,
                    emin=energies.emin,
                    emax=energies.emax,
                )

                spectrum_fs = {}
                for sat in satellites:
                    for name, instrument in sat:
                        if not instrument.use:
                            continue

                        filter_abbrv = instrument.filter_abbrv
                        fs = executor.submit(
                            get_spectrum_file,
                            instrument_name=name,
                            spectrum_dir=Path.cwd() / "res" / "spectrums",
                            filter_abbr=filter_abbrv,
                        )
                        spectrum_fs[fs] = name

                logger.info("START\tGetting spectrum files.")
                with tqdm(total=len(spectrum_fs), desc="Getting spectrum files") as pbar:
                    for future in as_completed(spectrum_fs):
                        name = spectrum_fs[future]
                        spectrum_file = future.result()
                        fs = executor.submit(
                            _background,
                            spectrum_file=spectrum_file,
                            instrument_name=name,
                        )
                        generator_fs[fs] = "bkg"
                        pbar.update()
                    elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
                logger.success(f"DONE\tGetting spectrum files. Duration: {elapsed_time}")

            if simput_cfg.agn.n_gen > 0:
                skip_agn = agn_counts_file is None or (agn_counts_file.exists() and agn_counts_file.is_dir())
                if skip_agn:
                    logger.warning(f"{agn_counts_file} does not exist! Won't create any AGN SIMPUTs.")
                else:
                    from src.simput.agn import create_agn
                    from src.xmm.tools import get_fov

                    logger.info("START\tGenerating SIMPUT for mode 'agn'...")
                    agn_path = simput_cfg.simput_dir / "agn"
                    agn_path.mkdir(parents=True, exist_ok=True)

                    logger.info(f"Will generate {simput_cfg.agn.n_gen} AGNs")
                    # Get the spectrum file
                    spectrum_file = get_spectrumfile(run_dir=tmp_dir, norm=0.001)

                    for _ in range(simput_cfg.agn.n_gen):
                        fs = executor.submit(
                            create_agn,
                            agn_counts_file=agn_counts_file,
                            emin=energies.emin,
                            emax=energies.emax,
                            fov=get_fov("epn"),
                            run_dir=Path(mkdtemp(dir=tmp_dir, prefix="bkg_")),
                            output_dir=agn_path,
                            xspec_file=spectrum_file,
                        )

                        generator_fs[fs] = "agn"

            img_total = countOf(generator_fs.values(), "img")
            bkg_total = countOf(generator_fs.values(), "bkg")
            agn_total = countOf(generator_fs.values(), "agn")

            pbars = {}
            tars = {}
            if img_total:
                pbars["img"] = tqdm(total=img_total, desc="Creating SIMPUTs for IMG")
                if env_cfg.tar_and_compress:
                    tars["img"] = (tarfile.open(simput_cfg.img_tar, "a"), simput_cfg.img_tar)

            if bkg_total:
                pbars["bkg"] = tqdm(total=bkg_total, desc="Creating SIMPUTs for BKG")
                if env_cfg.tar_and_compress:
                    tars["bkg"] = (tarfile.open(simput_cfg.bkg_tar, "a"), simput_cfg.bkg_tar)

            if agn_total:
                pbars["agn"] = tqdm(total=agn_total, desc="Creating SIMPUTs for AGN")
                if env_cfg.tar_and_compress:
                    tars["agn"] = (tarfile.open(simput_cfg.agn_tar, "a"), simput_cfg.agn_tar)

            for future in as_completed(generator_fs):
                mode = generator_fs[future]
                pbar: tqdm = pbars[mode]

                if mode in tars:
                    tar, tar_path = tars[mode]
                    for file in future.result():
                        logger.success(f"Created SIMPUT {file} for mode {mode.upper()}.")
                        tar.add(file, file.relative_to(simput_cfg.simput_dir))
                        logger.success(f"Added {file} to {tar_path}.")
                        if delete_product:
                            file.unlink()
                            logger.success(f"Deleted {file}")

                pbar.update()
                if pbar.format_dict["n"] == pbar.format_dict["total"]:
                    elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
                    logger.info(f"DONE\tGenerating SIMPUT for mode {mode.upper()}. Duration: {elapsed_time}")
                    pbar.close()
                    shutil.rmtree(simput_cfg.simput_dir / mode)
                    if mode in tars:
                        tar, _ = tars[mode]
                        tar.close()

            for mode, tar in tars.items():
                tar, tar_path = tar
                compressed = env_cfg.output_dir / "simput" / f"{mode}.tar.gz"
                compressed.parent.mkdir(parents=True, exist_ok=True)
                compressed.unlink(missing_ok=True)
                compress_gzip(tar_path, compressed, remove_file=True)
