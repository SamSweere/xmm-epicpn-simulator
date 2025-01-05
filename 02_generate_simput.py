import tomllib
from argparse import ArgumentParser
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from datetime import datetime
from itertools import islice
from pathlib import Path
from tempfile import TemporaryDirectory, gettempdir

import numpy as np
from loguru import logger
from tqdm import tqdm
from xspec import Model, Xset

from src.config import EnergyCfg, EnvironmentCfg, MultiprocessingCfg, SimputCfg
from src.tools.files import compress_targz, decompress_targz
from src.tools.run_utils import configure_logger, load_satellites
from src.xmm.tools import get_background_spectrum


def _generate_img_simputs(
    executor: ProcessPoolExecutor,
    simput_cfg: SimputCfg,
    env_cfg: EnvironmentCfg,
    energies: EnergyCfg,
) -> None:
    from src.simput.image import simput_image

    img_fs: dict[Future, Path] = {}

    if simput_cfg.fits_compressed:
        logger.info(f"Found compressed FITS files in {simput_cfg.fits_compressed}")
        decompress_targz(
            in_file_path=simput_cfg.fits_compressed,
            out_file_dir=simput_cfg.fits_dir,
            tar_options="--strip-components=1",
        )

    xspec_file = Path(gettempdir()) / "img_spectrum.xcm"
    Model("phabs*power", setPars={1: 0.04, 2: 2.0, 3: 0.01})
    Xset.save(f"{xspec_file}")
    fits_glob = simput_cfg.fits_dir.rglob("*.fits.gz")
    in_files = fits_glob if simput_cfg.img.n == -1 else islice(fits_glob, simput_cfg.img.n)
    for in_file in in_files:
        tng_set, snapshot_num = in_file.parts[-3], in_file.parts[-2]
        # Check how many files have already been generated and how many are left to generate
        simput_glob = (simput_cfg.img_dir / tng_set / snapshot_num).glob(f"{in_file.stem}*")
        already_created = len(list(islice(simput_glob, simput_cfg.img.num_samples)))
        missing = simput_cfg.img.num_samples - already_created

        # None are left to be generated => Skip
        if missing < 1:
            logger.info(
                f"Won't generate any images for {in_file.name}. "
                f"(Requested/Found: {simput_cfg.img.num_samples}/{already_created})"
            )
            if env_cfg.consume_data:
                in_file.unlink()
            continue

        logger.info(f"START\tGenerate {missing} images for {in_file.relative_to(simput_cfg.fits_dir)}")

        fs = executor.submit(
            simput_image,
            energies=energies,
            xspec_file=xspec_file,
            consume_data=env_cfg.consume_data,
            img_path_in=in_file,
            amount=missing,
            cfg=simput_cfg.img,
            output_dir=simput_cfg.img_dir / tng_set / snapshot_num,
        )

        img_fs[fs] = in_file

    for future in tqdm(as_completed(img_fs), desc="Creating SIMPUTs for IMG", total=len(img_fs)):
        out_files = future.result()
        in_file = img_fs[future]
        logger.success(f"Created {len(out_files)} SIMPUTs for {in_file.relative_to(simput_cfg.fits_dir)}")
    xspec_file.unlink()
    logger.success("DONE\tGenerating SIMPUT for mode IMG.")


def _generate_bkg_simputs(
    executor: ProcessPoolExecutor,
    simput_cfg: SimputCfg,
    env_cfg: EnvironmentCfg,
    energies: EnergyCfg,
    satellites: list,
) -> None:
    from src.simput.background import create_background

    bkg_fs: dict[Future, Path] = {}

    spectrum_fs = {}
    for sat in satellites:
        for name, instrument in sat:
            if not instrument.use:
                continue

            filter_abbrv = instrument.filter_abbrv
            fs = executor.submit(
                get_background_spectrum,
                instrument_name=name,
                spectrum_dir=env_cfg.working_dir / "spectrums",
                filter_abbr=filter_abbrv,
            )
            spectrum_fs[fs] = name

    logger.info("START\tGetting spectrum files.")
    for future in tqdm(as_completed(spectrum_fs), desc="Getting spectrum files", total=len(spectrum_fs)):
        spectrum_file = future.result()
        name = spectrum_fs[future]
        logger.success(f"Got spectrum file for {name} --> START Generation of BKG simputs.")
        fs = executor.submit(
            create_background,
            spectrum_file=spectrum_file,
            instrument_name=name,
            output_dir=simput_cfg.bkg_dir,
            energies=energies,
        )
        bkg_fs[fs] = name
    logger.success("DONE\tGetting spectrum files.")

    for _ in tqdm(as_completed(bkg_fs), desc="Creating SIMPUTs for BKG", total=len(bkg_fs)):
        # Do nothing. Logging is already taken care of and there is nothing else to do.
        pass
    logger.success("DONE\tGenerating SIMPUT for mode BKG.")


def _generate_agn_simputs(
    executor: ProcessPoolExecutor,
    simput_cfg: SimputCfg,
    energies: EnergyCfg,
    agn_counts_file: Path | None,
) -> None:
    skip_agn = agn_counts_file is None or (agn_counts_file.exists() and agn_counts_file.is_dir())
    if skip_agn:
        logger.warning(f"{agn_counts_file} does not exist! Won't create any AGN SIMPUTs.")
    else:
        from src.simput.agn import create_agn, get_s_n_from_file
        from src.xmm.tools import get_fov

        agn_fs = []

        logger.info("START\tGenerating SIMPUT for mode 'agn'...")

        logger.info(f"Will generate {simput_cfg.agn.n} AGNs")
        rng = np.random.default_rng()
        xspec_file = Path(gettempdir()) / "agn_spectrum.xcm"
        Model("phabs*power", setPars={1: 0.04, 2: 2.0, 3: 0.001})
        Xset.save(f"{xspec_file}")
        s, n = get_s_n_from_file(agn_counts_file)
        n = n * np.pi * 0.25**2
        d = np.flip(np.ediff1d(np.flip(n)))
        d_sum = np.sum(d)
        p = d / d_sum

        star_counts = np.round(d_sum + rng.uniform(-1, 1, simput_cfg.agn.n) * np.sqrt(d_sum)).astype(int)
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
                energies=energies,
                output_dir=simput_cfg.agn_dir,
                xspec_file=xspec_file,
            )

            agn_fs.append(fs)

        for future in tqdm(as_completed(agn_fs), desc="Creating SIMPUTs for AGN", total=len(agn_fs)):
            out_files = future.result()
            logger.success(f"Created AGN SIMPUT: {out_files}")
        xspec_file.unlink()
        logger.success("DONE\tGenerating SIMPUT for mode AGN.")


def generate_simputs(
    simput_cfg: SimputCfg,
    env_cfg: EnvironmentCfg,
    energies: EnergyCfg,
    mp_cfg: MultiprocessingCfg,
    satellites: list,
    agn_counts_file: Path | None,
) -> None:
    logger.info(f"Found satellites with instruments: {satellites}")

    starttime = datetime.now()
    with TemporaryDirectory(prefix="simput_") as tmp_dir:
        tmp_dir = Path(tmp_dir)

        with ProcessPoolExecutor(max_workers=mp_cfg.num_cores) as executor:
            # TODO Move download of background spectrums to the beginning
            #  That way this rather lengthy process can finish before
            #  the background simputs are to be created

            if simput_cfg.img.n != 0:
                # TODO Check if ccf files have been downloaded
                _generate_img_simputs(
                    executor=executor,
                    simput_cfg=simput_cfg,
                    env_cfg=env_cfg,
                    energies=energies,
                )
                if env_cfg.tar_and_compress:
                    logger.info(
                        f"IMG SIMPUTs will be compressed and moved to {simput_cfg.img_tar}. "
                        f"Existing file will be overwritten."
                    )
                    simput_cfg.img_tar.unlink(missing_ok=True)
                    executor.submit(
                        compress_targz,
                        in_path=simput_cfg.img_dir,
                        out_file_path=simput_cfg.img_tar,
                        remove_files=True,
                    )

            if simput_cfg.bkg.n:
                _generate_bkg_simputs(
                    executor=executor,
                    simput_cfg=simput_cfg,
                    env_cfg=env_cfg,
                    energies=energies,
                    satellites=satellites,
                )
                if env_cfg.tar_and_compress:
                    logger.info(
                        f"BKG SIMPUTs will be compressed and moved to {simput_cfg.bkg_tar}. "
                        f"Existing file will be overwritten."
                    )
                    simput_cfg.bkg_tar.unlink(missing_ok=True)
                    executor.submit(
                        compress_targz,
                        in_path=simput_cfg.bkg_dir,
                        out_file_path=simput_cfg.bkg_tar,
                        remove_files=True,
                    )

            if simput_cfg.agn.n > 0:
                _generate_agn_simputs(
                    executor=executor,
                    simput_cfg=simput_cfg,
                    energies=energies,
                    agn_counts_file=agn_counts_file,
                )
                if env_cfg.tar_and_compress:
                    logger.info(
                        f"AGN SIMPUTs will be compressed and moved to {simput_cfg.agn_tar} . "
                        f"Existing file will be overwritten."
                    )
                    simput_cfg.agn_tar.unlink(missing_ok=True)
                    executor.submit(
                        compress_targz,
                        in_path=simput_cfg.agn_dir,
                        out_file_path=simput_cfg.agn_tar,
                        remove_files=True,
                    )
    endtime = datetime.now()
    logger.info(f"Duration: {endtime - starttime}")


if __name__ == "__main__":
    parser = ArgumentParser(prog="", description="")
    parser.add_argument(
        "-a",
        "--agn_counts_file",
        default=Path(__file__).parent.resolve() / "res" / "agn_counts.cgi",
        type=Path,
        help="Path to agn_counts_cgi.",
    )
    parser.add_argument(
        "-p",
        "--config_path",
        type=Path,
        default=Path(__file__).parent.resolve() / "config.toml",
        help="Path to config file.",
    )

    args = parser.parse_args()

    with open(args.config_path, "rb") as file:
        cfg: dict[str, dict] = tomllib.load(file)

    environment = EnvironmentCfg(**cfg.pop("environment"))
    simput = SimputCfg(
        **cfg.pop("simput"),
        working_dir=environment.working_dir,
        output_dir=environment.output_dir,
    )
    energy = EnergyCfg(**cfg.pop("energy"))
    multiprocessing = MultiprocessingCfg(**cfg.pop("multiprocessing"))
    sats = load_satellites(cfg.pop("instruments"))

    logger.remove()

    configure_logger(
        log_dir=environment.log_dir,
        log_name="02_generate_simput_{time}.log",
        enqueue=True,
        debug=environment.debug,
        verbose=environment.verbose,
    )

    generate_simputs(
        simput_cfg=simput,
        env_cfg=environment,
        energies=energy,
        mp_cfg=multiprocessing,
        satellites=sats,
        agn_counts_file=args.agn_counts_file,
    )
