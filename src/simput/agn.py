from pathlib import Path
from uuid import uuid4

import matplotlib.pyplot as plt
import numpy as np
from loguru import logger

from src.simput.pointsource import create_pointsource
from src.simput.tools import merge_simputs
from src.tools.files import compress_gzip


def get_s_n_from_file(file_path: Path) -> tuple[np.ndarray, np.ndarray]:
    with open(file_path) as f:
        lines = f.readlines()
    # This is very specific for the file format, if the files changes recheck this
    lines = lines[7:27]
    vals = [line.strip().split() for line in lines]

    # S [erg / cm ** 2 / s]
    s = [float(x[0]) for x in vals]
    # N( > S) [deg ** -2]
    n = [float(x[1]) for x in vals]

    return np.array(s), np.array(n)


def get_fluxes(file_path: Path) -> np.ndarray:
    # S [erg / cm ** 2 / s]
    # N( > S) [deg ** -2]
    s, n = get_s_n_from_file(file_path)

    n = n * np.pi * 0.25**2  # correct for the xmm fov, xmm fov has r=15 arcmin = 0.25 degrees

    # Calculate the differences of the bins
    d = np.flip(np.ediff1d(np.flip(n)))
    d_sum = np.sum(d)

    # Take the number of stars and add some poisson error
    num_stars = round(d_sum + np.random.uniform(-1, 1) * np.sqrt(d_sum))

    # Normalize to get the probability per bin
    p = d / d_sum

    # Based on:
    # https://stackoverflow.com/questions/31730028/how-can-i-generate-a-random-sample-of-bin-counts-given-a-sequence-of-bin-probabi
    # Random sample_num from the bins based on the probabilities of getting a certain bin
    counts = np.bincount(np.random.choice(range(len(p)), size=num_stars, p=p), minlength=len(p))

    # Linearly interpolate and save the fluxes
    fluxes = []
    rng = np.random.default_rng()
    for i, count in enumerate(counts):
        s_min = s[i]
        s_max = s[i + 1]

        fluxes.append(rng.uniform(low=s_min, high=s_max, size=count))

    return np.concatenate(fluxes)


def create_agn(
    fluxes,
    offsets,
    emin: float,
    emax: float,
    run_dir: Path,
    output_dir: Path,
    xspec_file: Path,
) -> list[Path]:
    unique_id = uuid4().int
    final_name = f"agn_{unique_id}_p0_{emin}ev_p1_{emax}ev.simput.gz"
    out_file = output_dir / final_name
    simput_files: list[Path] = []

    for i, (flux, offset) in enumerate(zip(fluxes, offsets, strict=False)):
        logger.debug(f"Creating AGN with flux={flux}")
        output_file = run_dir / f"ps_{unique_id}_{i}.simput"
        compressed = output_file.with_suffix(".simput.gz")
        output_file = create_pointsource(
            emin=emin,
            emax=emax,
            output_file=output_file,
            src_flux=flux,
            xspec_file=xspec_file,
            offset=offset,
        )
        compress_gzip(output_file, compressed, remove_file=True)
        simput_files.append(compressed)
    merged = merge_simputs(simput_files=simput_files, output_file=run_dir / f"merged_{unique_id}.simput")
    compress_gzip(in_file_path=merged, out_file_path=out_file, remove_file=True)

    for file in simput_files:
        file.unlink(missing_ok=True)

    return [out_file]


def plot(file_path: Path, out_dir: Path):
    s, n = get_s_n_from_file(file_path)

    plt.plot(s, n, "g", label="data")

    plt.xlabel("S [erg / cm ** 2 / s]")
    plt.ylabel("N( > S) [deg ** -2]")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.savefig(out_dir / "xray_agn_number_count.pdf")
