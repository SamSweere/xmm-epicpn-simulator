import os

from utils.external_run import run_headas_command


def get_spectrumfile(run_dir, norm=0.01, verbose=True):
    # TODO: put this somewhere else
    # spectrum_file = "/home/sam/Documents/ESA/data/sim/test_spectrum.xcm"
    command = "xspec"
    norm = norm
    spectrum_file = os.path.join(run_dir, "spectrum.xcm")
    input = f"model phabs*power\n0.04\n2.0\n{str(norm)}\n save model {spectrum_file}"
    run_headas_command(command, input=input, verbose=verbose)

    return spectrum_file
