import shutil
import tarfile
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from functools import partial
from itertools import islice, repeat
from pathlib import Path
from tempfile import TemporaryDirectory, mkdtemp

import numpy as np
import requests
from loguru import logger
from tqdm import tqdm

from src.config import DownloadCfg, EnergyCfg, EnvironmentCfg, MultiprocessingCfg, SimputCfg, SimulationCfg
from src.illustris_tng.fits import cutout_to_xray_fits
from src.illustris_tng.web_api import (
    get_available_simulations,
    get_cutouts,
    get_subhalos,
)
from src.simput.tools import get_spectrumfile
from src.sixte.simulator import run_xmm_simulation
from src.tools.files import compress_gzip, compress_targz, decompress_targz
from src.xmm.tools import create_mask, create_psf_file, create_vinget_file, create_xml_files, get_spectrum_file


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

    with ProcessPoolExecutor(max_workers=mp_cfg.num_cores) as executor:
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
                cutout: Path = cutouts_fits_fs[future]
                logger.success(f"Converted {cutout} to X-ray FITS.")
                if "fits" in tars:
                    tar, tar_path = tars["fits"]
                    for fits in future.result():
                        tar.add(fits, fits.relative_to(download_cfg.fits_path))
                        logger.success(f"Added {fits} to {tar_path}.")
                        fits.unlink()
                        path.unlink(missing_ok=True)
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
) -> None:
    with TemporaryDirectory(prefix="simput_") as tmp_dir:
        tmp_dir = Path(tmp_dir)

        with ProcessPoolExecutor(max_workers=mp_cfg.num_cores) as executor:
            if simput_cfg.img.n_gen != 0:
                from src.simput.image import simput_image

                img_fs = {}

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
                fits_glob = simput_cfg.fits_dir.rglob("*.fits")
                in_files = fits_glob if amount_img == -1 else islice(fits_glob, amount_img)
                for in_file in in_files:
                    tng_set, snapshot_num = in_file.parts[-3], in_file.parts[-2]
                    # Check how many files have already been generated and how many are left to generate
                    simput_glob = (img_path / tng_set / snapshot_num).glob(f"{in_file.stem}*")
                    already_created = len(list(islice(simput_glob, simput_cfg.num_img_sample)))
                    missing = simput_cfg.num_img_sample - already_created

                    # None are left to be generated => Skip
                    if missing == 0:
                        logger.debug(f"Won't generate any images for {in_file.name}.")
                        if env_cfg.consume_data:
                            in_file.unlink()
                        continue

                    logger.info(f"Will generate {missing} images for {in_file.name}.")

                    zoom = np.round(
                        rng.uniform(
                            simput_cfg.zoom_range[0],
                            simput_cfg.zoom_range[1],
                            missing,
                        ),
                        2,
                    )
                    sigma_b = np.round(
                        rng.uniform(
                            simput_cfg.sigma_b_range[0],
                            simput_cfg.sigma_b_range[1],
                            missing,
                        ),
                        2,
                    )

                    offset_x = np.round(
                        rng.normal(
                            -simput_cfg.offset_std,
                            simput_cfg.offset_std,
                            missing,
                        ),
                        2,
                    )
                    offset_y = np.round(
                        rng.normal(
                            -simput_cfg.offset_std,
                            simput_cfg.offset_std,
                            missing,
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

                    img_fs[fs] = in_file

                img_tar = tarfile.open(simput_cfg.img_tar, "a") if env_cfg.tar_and_compress else None

                with tqdm(total=len(img_fs), desc="Creating SIMPUTs for IMG") as pbar:
                    for future in as_completed(img_fs):
                        out_files = future.result()
                        in_file = img_fs[future]
                        logger.success(f"Created {len(out_files)} SIMPUTs for {in_file}.")
                        if img_tar is not None:
                            for out_file in out_files:
                                img_tar.add(out_file, out_file.relative_to(simput_cfg.simput_dir))
                                out_file.unlink()
                        pbar.update()
                    elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
                logger.success(f"DONE\tGenerating SIMPUT for mode IMG. Duration: {elapsed_time}")

                if img_tar is not None:
                    img_tar.close()
                    shutil.rmtree(img_path)
                    executor.submit(
                        compress_gzip,
                        in_file_path=simput_cfg.img_tar,
                        out_file_path=env_cfg.output_dir / "simput" / "img.tar.gz",
                        remove_file=True,
                    )

            if simput_cfg.bkg.n_gen:
                from src.simput.background import create_background

                bkg_fs = {}

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
                            spectrum_dir=env_cfg.working_dir / "spectrums",
                            filter_abbr=filter_abbrv,
                        )
                        spectrum_fs[fs] = name

                logger.info("START\tGetting spectrum files.")
                with tqdm(total=len(spectrum_fs), desc="Getting spectrum files") as pbar:
                    for future in as_completed(spectrum_fs):
                        spectrum_file = future.result()
                        name = spectrum_fs[future]
                        fs = executor.submit(
                            _background,
                            spectrum_file=spectrum_file,
                            instrument_name=name,
                        )
                        bkg_fs[fs] = name
                        pbar.update()
                    elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
                logger.success(f"DONE\tGetting spectrum files. Duration: {elapsed_time}")

                bkg_tar = tarfile.open(simput_cfg.bkg_tar, "a") if env_cfg.tar_and_compress else None

                with tqdm(total=len(bkg_fs), desc="Creating SIMPUTs for BKG") as pbar:
                    for future in as_completed(bkg_fs):
                        out_files = future.result()
                        name = bkg_fs[future]
                        logger.success(f"Created BKG SIMPUT for {name}.")
                        if bkg_tar is not None:
                            for out_file in out_files:
                                bkg_tar.add(out_file, out_file.relative_to(simput_cfg.simput_dir))
                                out_file.unlink()
                        pbar.update()
                    elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
                logger.success(f"DONE\tGenerating SIMPUT for mode BKG. Duration: {elapsed_time}")

                if bkg_tar is not None:
                    bkg_tar.close()
                    shutil.rmtree(bkg_path)
                    executor.submit(
                        compress_gzip,
                        in_file_path=simput_cfg.bkg_tar,
                        out_file_path=env_cfg.output_dir / "simput" / "bkg.tar.gz",
                        remove_file=True,
                    )

            if simput_cfg.agn.n_gen > 0:
                skip_agn = agn_counts_file is None or (agn_counts_file.exists() and agn_counts_file.is_dir())
                if skip_agn:
                    logger.warning(f"{agn_counts_file} does not exist! Won't create any AGN SIMPUTs.")
                else:
                    from src.simput.agn import create_agn, get_s_n_from_file
                    from src.xmm.tools import get_fov

                    agn_fs = []

                    logger.info("START\tGenerating SIMPUT for mode 'agn'...")
                    agn_path = simput_cfg.simput_dir / "agn"
                    agn_path.mkdir(parents=True, exist_ok=True)

                    logger.info(f"Will generate {simput_cfg.agn.n_gen} AGNs")
                    rng = np.random.default_rng()
                    # Get the spectrum file
                    spectrum_file = get_spectrumfile(run_dir=tmp_dir, norm=0.001)
                    s, n = get_s_n_from_file(agn_counts_file)
                    n = n * np.pi * 0.25**2
                    d = np.flip(np.ediff1d(np.flip(n)))
                    d_sum = np.sum(d)
                    p = d / d_sum

                    star_counts = np.round(d_sum + rng.uniform(-1, 1, simput_cfg.agn.n_gen) * np.sqrt(d_sum)).astype(
                        int
                    )
                    fov = get_fov("epn")

                    for star_count in star_counts:
                        counts = np.bincount(rng.choice(range(len(p)), size=star_count, p=p), minlength=len(p))
                        indices = np.flatnonzero(counts)
                        fluxes = [rng.uniform(low=s[i], high=s[i + 1], size=counts[i]) for i in indices]
                        fluxes = np.concatenate(fluxes)
                        offsets = rng.uniform(low=-fov / 2.0, high=fov / 2.0, size=(fluxes.shape[0], 2))
                        fs = executor.submit(
                            create_agn,
                            fluxes=fluxes,
                            offsets=offsets,
                            emin=energies.emin,
                            emax=energies.emax,
                            run_dir=Path(mkdtemp(dir=tmp_dir, prefix="agn_")),
                            output_dir=agn_path,
                            xspec_file=spectrum_file,
                        )

                        agn_fs.append(fs)

                    agn_tar = tarfile.open(simput_cfg.agn_tar, "a") if env_cfg.tar_and_compress else None

                    with tqdm(total=len(agn_fs), desc="Creating SIMPUTs for AGN") as pbar:
                        for future in as_completed(agn_fs):
                            out_files = future.result()
                            for out_file in out_files:
                                logger.success(f"Created AGN SIMPUT {out_file}.")
                                if agn_tar is not None:
                                    agn_tar.add(out_file, out_file.relative_to(simput_cfg.simput_dir))
                                    out_file.unlink()
                            pbar.update()
                        elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
                    logger.success(f"DONE\tGenerating SIMPUT for mode AGN. Duration: {elapsed_time}")

                    if agn_tar is not None:
                        agn_tar.close()
                        shutil.rmtree(bkg_path)
                        executor.submit(
                            compress_gzip,
                            in_file_path=simput_cfg.agn_tar,
                            out_file_path=env_cfg.output_dir / "simput" / "agn.tar.gz",
                            remove_file=True,
                        )


def run_simulations(
    sim_cfg: SimulationCfg,
    energies: EnergyCfg,
    env_cfg: EnvironmentCfg,
    mp_cfg: MultiprocessingCfg,
    satellites: list,
) -> None:
    root_dir = env_cfg.working_dir
    with (
        TemporaryDirectory(prefix="xml_", dir=root_dir) as xml_dir,
        TemporaryDirectory(prefix="sim_", dir=root_dir) as sim_dir,
    ):
        xml_dir = Path(xml_dir)
        sim_dir = Path(sim_dir)

        # Create all needed directories
        for sat in satellites:
            for name, instrument in sat:
                if not instrument.use:
                    continue
                filter_dir: Path = xml_dir / name / instrument.filter
                for res_mult in sim_cfg.res_mults:
                    (filter_dir / f"{res_mult}x").mkdir(exist_ok=True, parents=True)

        with ProcessPoolExecutor(max_workers=mp_cfg.num_cores) as executor:
            # Decompress SIMPUT files if needed
            if env_cfg.working_dir != env_cfg.output_dir:
                for mode, amount in sim_cfg.modes:
                    if amount == 0:
                        logger.debug(f"Skipping {mode} since simulation amount is set to 0.")
                        continue

                    simput_dir = env_cfg.output_dir / "simput"

                    simput_compressed_files = [next(simput_dir.rglob("*.tar.gz"))]

                    for simput_compressed in simput_compressed_files:
                        if simput_compressed.exists():
                            logger.info(f"START\tDecompressing SIMPUT files in {simput_compressed.resolve()}.")
                            executor.submit(
                                decompress_targz,
                                in_file_path=simput_compressed,
                                out_file_dir=sim_cfg.simput_dir / mode,
                                tar_options="--strip-components=1",
                            )
            emask_fs = []
            logger.info("START\tCreating all the needed files.")
            for sat in satellites:
                for name, instrument in sat:
                    if instrument.use:
                        # Vignetting files
                        executor.submit(
                            create_vinget_file,
                            instrument_name=name,
                            xml_dir=xml_dir,
                        )
                        # EMASKS
                        fs = executor.submit(
                            create_mask,
                            instrument_name=name,
                            observation_id="0935190401",
                            mask_level=instrument.mask_level,
                            energies=energies,
                            out_dir=sim_dir,
                            res_mults=sim_cfg.res_mults,
                        )
                        emask_fs.append(fs)
                        for res_mult in sim_cfg.res_mults:
                            # PSF files
                            executor.submit(
                                create_psf_file,
                                instrument_name=name,
                                xml_dir=xml_dir,
                                res_mult=res_mult,
                            )
                            # XML files
                            executor.submit(
                                create_xml_files,
                                instrument_name=name,
                                xml_dir=xml_dir,
                                res_mult=res_mult,
                                xmm_filter=instrument.filter,
                                sim_separate_ccds=instrument.sim_separate_ccds,
                                wait_time=sim_cfg.wait_time,
                            )

        emasks = {}
        for fs in emask_fs:
            for key, value in fs.result().items():
                emasks[key] = value

        for sat in satellites:
            for name, instrument in sat:
                if instrument.use:
                    max_workers = mp_cfg.ram_gb // 8 if name == "epn" else mp_cfg.num_cores
                    xmm_filter_dir = sim_cfg.out_dir / name / instrument.filter
                    xmm_filter_dir.mkdir(exist_ok=True, parents=True)
                    with ProcessPoolExecutor(max_workers=max_workers) as executor:
                        for mode, amount in sim_cfg.modes:
                            if amount == 0:
                                logger.info(f"Skipping {mode.upper()} simulation since amount is 0.")
                                continue
                            mode_fs = {}
                            logger.info(f"START\tSimulating {name} for {mode.upper()}.")

                            # Find the simput files
                            mode_dir = sim_cfg.simput_dir / mode
                            if mode != "bkg":
                                mode_glob = mode_dir.rglob("*.simput.gz")
                                simputs = mode_glob if amount == -1 else islice(mode_glob, amount)
                            else:
                                simputs = repeat(next(mode_dir.rglob(f"*{name}.simput.gz")), amount)

                            for simput in simputs:
                                for res_mult in sim_cfg.res_mults:
                                    fs: Future = executor.submit(
                                        run_xmm_simulation,
                                        instrument_name=name,
                                        xml_dir=xml_dir,
                                        simput_file=simput,
                                        mode=mode,
                                        tmp_dir=sim_dir,
                                        out_dir=xmm_filter_dir,
                                        res_mult=res_mult,
                                        max_event_pattern=instrument.max_event_pattern,
                                        exposure=sim_cfg.max_exposure,
                                        xmm_filter=instrument.filter,
                                        sim_separate_ccds=instrument.sim_separate_ccds,
                                        consume_data=env_cfg.consume_data,
                                        emask=emasks[name][res_mult],
                                    )
                                    mode_fs[fs] = {"simput": simput, "res_mult": res_mult}

                            for future in tqdm(
                                as_completed(mode_fs), total=len(mode_fs), desc=f"Simulating {name} for {mode.upper()}"
                            ):
                                # Since this feature should be done, add a small timeout
                                outfiles = future.result(10)
                                simput = mode_fs[future]["simput"]
                                res_mult = mode_fs[future]["res_mult"]
                                logger.success(f"Simulated {name} for {simput} with res_mult {res_mult}.")
                                logger.info(f"Created {len(outfiles)} images")
                            logger.success(f"DONE\tSimulating {name} for {mode.upper()}. Duration: elapsed_time")

                            if env_cfg.tar_and_compress:
                                mode_compressed = (
                                    env_cfg.output_dir / "xmm_sim_dataset" / name / instrument.filter / f"{mode}.tar.gz"
                                )
                                mode_compressed.parent.mkdir(parents=True, exist_ok=True)
                                executor.submit(
                                    compress_targz,
                                    in_path=xmm_filter_dir / mode,
                                    out_file_path=mode_compressed,
                                    remove_files=True,
                                )
