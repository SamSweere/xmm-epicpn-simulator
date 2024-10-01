import pathlib
import tomllib
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from tempfile import TemporaryDirectory

from loguru import logger
from tqdm import tqdm

from src.config import EnvironmentCfg, MultiprocessingCfg, SimulationCfg
from src.heasoft import heasoft as hsp
from src.tools.run_utils import configure_logger, load_satellites

logger.remove()


def _do(infile: Path, outfile: Path, exposure: int, rebin: bool) -> None:
    outfile.parent.mkdir(exist_ok=True, parents=True)
    with TemporaryDirectory(prefix="pre_") as tmp_dir:
        tmp_dir = Path(tmp_dir)
        infile = hsp.ftcopy(infile=infile, outfile=tmp_dir / f"in{''.join(infile.suffixes)}")

        if rebin:
            infile = hsp.fimgbin(infile=infile, outfile=infile, xbinsize=2)

        hsp.ftimgcalc(outfile=outfile, expr=f"A / {exposure}", A=infile)


def preprocess(
    env_cfg: EnvironmentCfg,
    satellites: list,
):
    out_dir = env_cfg.output_dir / "xmm_sim_dataset" / "scaled"

    for satellite in satellites:
        for name, instrument in satellite:
            instrument_dir = env_cfg.output_dir / "xmm_sim_dataset" / name / instrument.filter
            out_instrument_dir = out_dir / name / instrument.filter

            rebin = name != "epn"

            for mode_dir in instrument_dir.iterdir():
                if not mode_dir.is_dir():
                    continue

                out_mode_dir = out_instrument_dir / mode_dir.name

                for subdir in mode_dir.iterdir():
                    if not subdir.is_dir():
                        continue

                    exposure = int(subdir.name.split("ks")[0]) * 1000

                    with ProcessPoolExecutor(max_workers=mp_cfg.num_cores) as executor:
                        fs = []
                        for fits in subdir.rglob("*.fits.gz"):
                            out_fits = out_mode_dir / subdir.name / fits.relative_to(subdir)

                            if out_fits.exists():
                                continue

                            fs.append(
                                executor.submit(
                                    _do,
                                    infile=fits,
                                    outfile=out_fits,
                                    exposure=exposure,
                                    rebin=rebin,
                                )
                            )

                        if not fs:
                            continue

                        description = f"Preprocessing {name} for {mode_dir.name.upper()} with exposure {exposure}s"

                        for future in tqdm(as_completed(fs), total=len(fs), desc=description):
                            exception = future.exception()

                            if exception:
                                print(exception)
                                logger.exception(exception)
                                executor.shutdown(cancel_futures=True)
                                raise exception

                        del fs


if __name__ == "__main__":
    parser = ArgumentParser(prog="", description="")
    parser.add_argument(
        "-p",
        "--config_path",
        type=Path,
        default=pathlib.Path(__file__).parent.resolve() / "config.toml",
        help="Path to config file.",
    )

    args = parser.parse_args()

    with open(args.config_path, "rb") as file:
        cfg: dict[str, dict] = tomllib.load(file)

    env_cfg = EnvironmentCfg(**cfg.pop("environment"))
    sim_cfg = SimulationCfg(
        **cfg.pop("simulation"),
        simput_dir=env_cfg.working_dir / "simput",
        out_dir=env_cfg.working_dir / "xmm_sim_dataset",
    )
    mp_cfg = MultiprocessingCfg(**cfg.pop("multiprocessing"))

    satellites = load_satellites(cfg.pop("instruments"))

    del cfg

    configure_logger(
        log_dir=env_cfg.log_dir,
        log_name="04_preprocess.log",
        enqueue=True,
        debug=env_cfg.debug,
        verbose=env_cfg.verbose,
    )

    preprocess(
        env_cfg=env_cfg,
        satellites=satellites,
    )
