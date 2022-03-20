import os

import numpy as np


def get_S_N_from_file():
    f = open(os.path.join(os.path.dirname(__file__), "agn_counts.cgi"), "r")

    # This is very specific for the file format, if the files changes recheck this
    lines = f.readlines()
    lines = lines[7:27]
    lines = [x[2:] for x in lines]
    lines = [x.replace('\n', '') for x in lines]
    vals = [x.split(" ") for x in lines]

    # S [erg / cm ** 2 / s]
    S = [float(x[0]) for x in vals]
    # N( > S) [deg ** -2]
    N = [float(x[-1]) for x in vals]

    return np.array(S), np.array(N)


def get_fluxes():
    # S [erg / cm ** 2 / s]
    # N( > S) [deg ** -2]
    S, N = get_S_N_from_file()

    N = N * np.pi * 0.25 ** 2  # correct for the xmm fov, xmm fov has r=15 arcmin = 0.25 degrees

    # Calculate the differences of the bins
    D = []
    for i in range(len(N) - 1):
        D.append(N[i] - N[i + 1])

    # Take the number of stars and add some poisson error
    num_stars = round(sum(D) + np.random.uniform(-1, 1) * np.sqrt(sum(D)))

    # Normalize to get the probability per bin
    P = D / np.sum(D)

    # Based on: https://stackoverflow.com/questions/31730028/how-can-i-generate-a-random-sample-of-bin-counts-given-a-sequence-of-bin-probabi
    # Random sample_num from the bins based on the probabilities of getting a certain bin
    counts = np.bincount(np.random.choice(range(len(P)), size=num_stars, p=P), minlength=len(P))

    # Linearly interpolate and save the fluxes
    fluxes = []
    for i in range(len(counts)):
        c = counts[i]
        s_min = S[i]
        s_max = S[i + 1]

        for j in range(c):
            s = np.random.uniform(s_min, s_max)
            fluxes.append(s)

    # print("Fluxes:",fluxes)
    # print("Max flux:",max(fluxes))

    return fluxes
