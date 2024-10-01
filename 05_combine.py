import pathlib
import tomllib
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from tempfile import TemporaryDirectory

from loguru import logger
from tqdm import tqdm

from src.config import EnergyCfg, EnvironmentCfg, MultiprocessingCfg, SimulationCfg
from src.heasoft import heasoft as hsp
from src.tools.run_utils import configure_logger, load_satellites
from src.xmm.epn import get_shift_xy
from src.xmm.tools import create_mask

logger.remove()


def _do(
    epn_file: Path,
    epn_dir: Path,
    mask: Path,
    shift_x: float,
    shift_y: float,
    tmp_root: Path,
    out_dir: Path,
) -> None:
    logger.info(f"Working on {epn_file}")
    emos1_dir: Path = env_cfg.output_dir / "xmm_sim_dataset" / "scaled" / "emos1" / "thin"
    emos2_dir: Path = env_cfg.output_dir / "xmm_sim_dataset" / "scaled" / "emos2" / "thin"
    parts = list(epn_file.relative_to(epn_dir).parts)

    final_out = out_dir / Path(*parts)

    parts[-2] = "1x"
    new = parts[-1].replace("_mult_2_", "_mult_1_")
    parts[-1] = new

    emos1_file = emos1_dir / Path(*parts)
    emos2_file = emos2_dir / Path(*parts)
    assert emos1_file.exists()
    # logger.info(f"Found {emos1_file}")
    assert emos2_file.exists()
    # logger.info(f"Found {emos2_file}")

    with TemporaryDirectory(dir=tmp_root) as tmp_dir:
        tmp_dir = Path(tmp_dir)

        epn_file = hsp.ftcopy(epn_file, tmp_dir / "epn.fits.gz")
        emos1_file = hsp.ftcopy(emos1_file, tmp_dir / "emos1.fits.gz")
        emos2_file = hsp.ftcopy(emos2_file, tmp_dir / "emos2.fits.gz")

        # Merge images
        xoffset_emos1 = -99 + shift_x
        yoffset_emos1 = 82 + shift_y
        xoffset_emos2 = -95.5 + shift_x
        yoffset_emos2 = -89.5 + shift_y

        outfile = hsp.fimgmerge(
            epn_file,
            [emos1_file, emos2_file],
            tmp_dir / "out.fits.gz",
            [xoffset_emos1, xoffset_emos2],
            [yoffset_emos1, yoffset_emos2],
        )
        # logger.info(f"Added {emos1_file} and {emos2_file}")

        # Add detmask of epn and save image
        final_out.parent.mkdir(exist_ok=True, parents=True)
        hsp.ftimgcalc(final_out, "A * B", a=outfile, b=mask)
        # logger.info(f"Added detmask")


def combine() -> None:
    out_dir = env_cfg.output_dir / "xmm_sim_dataset" / "combined"

    epn_cfg = satellites[0].epn
    epn_dir: Path = env_cfg.output_dir / "xmm_sim_dataset" / "scaled" / "epn" / epn_cfg.filter
    assert epn_dir.exists()

    with TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        # Create mask
        masks = create_mask(
            instrument_name="epn",
            observation_id="0935190401",
            mask_level="expmap",
            energies=energies,
            out_dir=tmp_dir,
            res_mults=[2],
        )
        mask = masks["epn"][2]

        shift_x, shift_y = get_shift_xy(2)

        for mode, _ in sim_cfg.modes:
            if mode == "bkg":
                # Do not merge background images. These should be randomly picked during training
                continue

            # Find EPN images
            epn_files = epn_dir.rglob(f"{mode}/**/2x/*fits.gz")

            with ProcessPoolExecutor(max_workers=mp_cfg.num_cores) as executor:
                fs = []
                for epn_file in epn_files:
                    fs.append(
                        executor.submit(
                            _do,
                            epn_file=epn_file,
                            epn_dir=epn_dir,
                            mask=mask,
                            tmp_root=tmp_dir,
                            out_dir=out_dir,
                            shift_x=shift_x,
                            shift_y=shift_y,
                        )
                    )

                description = f"Merging HR images for mode {mode.upper()}"

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
    energies = EnergyCfg(**cfg.pop("energy"))
    satellites = load_satellites(cfg.pop("instruments"))

    del cfg

    configure_logger(
        log_dir=env_cfg.log_dir,
        log_name="05_comb_sim.log",
        enqueue=True,
        debug=env_cfg.debug,
        verbose=env_cfg.verbose,
    )

    combine()
