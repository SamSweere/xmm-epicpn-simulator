import random
import shutil
from pathlib import Path
from typing import List

from src.utils.external_run import run_headas_command


def get_spectrumfile(run_dir: Path, norm=0.01, verbose=True):
    # TODO: put this somewhere else
    # spectrum_file = "/home/sam/Documents/ESA/data/sim/test_spectrum.xcm"
    command = "xspec"
    spectrum_file = run_dir / "spectrum.xcm"
    input = f"model phabs*power\n0.04\n2.0\n{str(norm)}\n save model {spectrum_file.resolve()}"
    run_headas_command(command, input=input, verbose=verbose)

    return spectrum_file


def merge_and_save_simputs(
        run_dir: Path,
        simput_files: List[Path],
        output_file: Path,
        verbose=True
):
    # Combine the simput point sources
    if len(simput_files) == 1:
        shutil.copy(simput_files[0], output_file)
    else:
        files = [(simput_files[0], simput_files[1])]
        files.extend([(simput_file, output_file) for simput_file in simput_files[2:]])
        for infile1, infile2 in files:
            # Merge the simput files
            merge_command = f"simputmerge FetchExtensions=yes Infile1={infile1.resolve()} Infile2={infile2.resolve()}" \
                            f" Outfile={output_file.resolve()}"

            if verbose:
                print(infile1)
                print(infile2)
                print(f"{output_file.resolve()}")
                print(merge_command)
                print("------------------------")

            run_headas_command(merge_command, run_dir=run_dir, verbose=verbose)
    return output_file


def get_simputs(
        simput_path: Path,
        mode,
        amount=-1,
        order='normal'
):
    # Order options: normal (front to back), reversed (back to front), random
    # Put the mode and amount in a list if they are passed as single values
    if type(mode) != list:
        mode = [mode]

    if type(amount) != list:
        amount = [amount]

    if not len(amount) == len(mode):
        raise AssertionError(f"the number of modes ({len(mode)}) should be equal "
                             f"to the number of amounts ({len(amount)})")

    simput_files = []
    for m, n in zip(mode, amount):
        if n == 0:
            continue

        mode_path = simput_path / m
        files = mode_path.glob("*.simput.gz")

        if n == -1:
            # Do all
            tmp = list(files)
        else:
            tmp = []
            for file_count, file in enumerate(files):
                if file_count < n:
                    tmp.append(file)
                else:
                    break

        simput_files.append(tmp)

    if order == "normal":
        simput_files = simput_files
    elif order == "reversed":
        simput_files = simput_files.reverse()
    elif order == "random":
        random.shuffle(simput_files)
    else:
        raise ValueError(f'Order: {order} not in known options list of "normal", "reversed", "random"')

    return simput_files
