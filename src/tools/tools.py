from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from itertools import islice, repeat
from pathlib import Path
from tempfile import TemporaryDirectory

from loguru import logger
from requests import Session
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util import Retry

from src.config import EnergyCfg, EnvironmentCfg, MultiprocessingCfg, SimulationCfg
from src.sixte.simulator import run_simulation
from src.tools.files import compress_targz, decompress_targz
from src.xmm.tools import create_mask, create_psf_file, create_vinget_file, create_xml_files


def download_file(url: str, out_path: Path) -> Path:
    retry_strategy = Retry(total=10, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session = Session()
    prefix = "https://" if url.startswith("https://") else "http://"
    session.mount(prefix, adapter)

    with session.get(url, stream=True) as response, open(out_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=int(1e6)):
            f.write(chunk)

    return out_path


def run_simulations(
    sim_cfg: SimulationCfg,
    energies: EnergyCfg,
    env_cfg: EnvironmentCfg,
    mp_cfg: MultiprocessingCfg,
    satellites: list,
) -> None:
    with TemporaryDirectory(prefix="xml_") as xml_dir, TemporaryDirectory(prefix="sim_") as sim_dir:
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

                    simput_compressed = next(simput_dir.rglob(f"{mode}.tar.gz"))

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

        for mode, amount in sim_cfg.modes:
            if amount == 0:
                logger.info(f"Skipping {mode.upper()} simulation since amount is 0.")
                continue

            mode_dir = sim_cfg.simput_dir / mode
            pattern = "*.simput.gz" if mode != "bkg" else f"*_{{0}}_{energies.emin}keV_{energies.emax}keV.simput.gz"

            for sat in satellites:
                for name, instrument in sat:
                    if not instrument.use:
                        logger.info(f"Skipping {name} simulation for {mode.upper()} since 'use' is False")
                        continue

                    if mode != "bkg":
                        mode_glob = mode_dir.rglob(pattern)
                        simputs = mode_glob if amount == -1 else islice(mode_glob, amount)
                    else:
                        simputs = repeat(next(mode_dir.rglob(pattern.format(name))), amount)

                    mode_fs = {}
                    logger.info(f"START\tSimulating {name} for {mode.upper()}.")
                    ram_max = mp_cfg.ram_gb // 4 if name == "epn" else mp_cfg.ram_gb // 2
                    max_workers = min(mp_cfg.num_cores, ram_max)

                    logger.debug(f"Will use {max_workers} processes.")

                    xmm_filter_dir = sim_cfg.out_dir / name / instrument.filter
                    xmm_filter_dir.mkdir(exist_ok=True, parents=True)
                    with ProcessPoolExecutor(max_workers=max_workers) as executor:
                        for simput in simputs:
                            fs: Future = executor.submit(
                                run_simulation,
                                tmp_dir=sim_dir,
                                out_dir=xmm_filter_dir / mode,
                                xml_dir=xml_dir,
                                instrument_name=name,
                                xmm_filter=instrument.filter,
                                simput_path=simput,
                                res_mults=sim_cfg.res_mults,
                                exposure=sim_cfg.max_exposure,
                                max_event_pattern=instrument.max_event_pattern,
                                mode=mode,
                                sim_separate_ccds=instrument.sim_separate_ccds,
                                consume_data=env_cfg.consume_data,
                                emasks=emasks[name],
                            )
                            mode_fs[fs] = {"simput": simput}

                        for future in tqdm(
                            as_completed(mode_fs), total=len(mode_fs), desc=f"Simulating {name} for {mode.upper()}"
                        ):
                            exception = future.exception()

                            if exception:
                                print(exception)
                                logger.exception(exception)
                                executor.shutdown(cancel_futures=True)
                                raise exception

                            outfiles = future.result()
                            simput = mode_fs[future]["simput"]
                            logger.success(f"Simulated {name} for {simput}.")
                            logger.info(f"Created {len(outfiles)} images")
                        logger.success(f"DONE\tSimulating {name} for {mode.upper()}. Duration: elapsed_time")

                        mode_fs.clear()

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
