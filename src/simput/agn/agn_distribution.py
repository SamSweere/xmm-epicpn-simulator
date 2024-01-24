from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np


def get_s_n_from_file(file_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    with open(file_path, "r") as f:
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

    n = (
        n * np.pi * 0.25**2
    )  # correct for the xmm fov, xmm fov has r=15 arcmin = 0.25 degrees

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
    counts = np.bincount(
        np.random.choice(range(len(p)), size=num_stars, p=p), minlength=len(p)
    )

    # Linearly interpolate and save the fluxes
    fluxes = []
    rng = np.random.default_rng()
    for i, count in enumerate(counts):
        s_min = s[i]
        s_max = s[i + 1]

        fluxes.append(rng.uniform(low=s_min, high=s_max, size=count))

    return np.concatenate(fluxes)


def plot(file_path: Path, out_dir: Path):
    s, n = get_s_n_from_file(file_path)

    plt.plot(s, n, "g", label="data")

    plt.xlabel("S [erg / cm ** 2 / s]")
    plt.ylabel("N( > S) [deg ** -2]")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.savefig(out_dir / "xray_agn_number_count.pdf")
