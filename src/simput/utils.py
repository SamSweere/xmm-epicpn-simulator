import random
import shutil
from pathlib import Path
from typing import List, Dict
from src.xmm_utils.external_run import run_headas_command
from xspec import Model, Xset


def get_spectrumfile(run_dir: Path, norm=0.01, verbose=True):
    spectrum_file = run_dir / "spectrum.xcm"
    if not spectrum_file.exists():
        Model("phabs*power", setPars={1: 0.04, 2: 2.0, 3: norm})
        Xset.save(f"{spectrum_file.resolve()}")
    return spectrum_file


def _simput_merge(
        infiles: List[Path],
        outfile: Path,
        verbose: bool = True
) -> None:
    str_infiles = [str(infile.resolve()) for infile in infiles]
    str_infiles = ",".join(str_infiles)

    merge_command = f"simputmerge FetchExtensions=yes Infiles={str_infiles} Outfile={outfile.resolve()}"

    run_headas_command(merge_command, verbose=verbose)


def merge_simputs(
        simput_files: List[Path],
        output_file: Path,
        keep_files: bool = False,
        verbose=True
) -> Path:
    # Combine the simput point sources
    if len(simput_files) == 1:
        file = simput_files[0]
        if keep_files:
            shutil.copy2(file, output_file)
        else:
            file.rename(output_file)
    else:
        _simput_merge(simput_files, output_file, verbose=verbose)

        if not keep_files:
            for file in simput_files:
                file.unlink()

    return output_file


def _order(iterable: List, order: str) -> List:
    if order == "normal":
        pass
    elif order == "reversed":
        iterable.reverse()
    elif order == "random":
        random.shuffle(iterable)
    else:
        raise ValueError(f'Order: {order} not in known options list of "normal", "reversed", "random"')

    return iterable


def get_simputs(
        simput_path: Path,
        mode_amount_dict: Dict[str, int],
        order='normal'
) -> Dict[str, List[Path]]:
    # Order options: normal (front to back), reversed (back to front), random

    simput_files: Dict[str, List[Path]] = {}
    for mode, amount in mode_amount_dict.items():
        if amount == 0:
            continue

        mode_path = simput_path / mode
        files = mode_path.glob("*.simput.gz")

        if amount == -1:
            # Do all
            simput_files[mode] = _order(list(files), order)
        else:
            tmp = []
            for file_count, file in enumerate(files):
                if file_count < amount:
                    tmp.append(file)
                else:
                    break
            simput_files[mode] = _order(tmp, order)

    return simput_files
