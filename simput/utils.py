import os
import random
import shutil

from utils.external_run import run_headas_command
from utils.log import plog


def copy_to_save_folder(simput_file_path, run_dir, output_file_path, verbose=True):
    shutil.copy(os.path.join(run_dir, simput_file_path), output_file_path)
    plog(f"Saved final simput file at: {output_file_path}", verbose=verbose)

    return output_file_path


def merge_and_save_simputs(run_dir, simput_files, output_file_path, verbose=True):
    # Combine the simput point sources
    merged_files = []

    for i in range(1, len(simput_files)):
        if i == 1:
            infile1 = simput_files[0]
            infile2 = simput_files[i]
        else:
            infile1 = simput_files[i]
            infile2 = merged_files[-1]

        merge_file = os.path.join(run_dir, f"merged_{i}.simput")
        merged_files.append(merge_file)

        # Merge the simput files
        # The cd is necessary for an unkown reason. Otherwise it will a file access error
        merge_command = f"cd {run_dir} && simputmerge FetchExtensions=yes Infile1={infile1} Infile2={infile2} " \
                        f"Outfile={merge_file}"
        run_headas_command(merge_command, verbose=verbose)

        if verbose:
            print(infile1)
            print(infile2)
            print(merge_file)
            print(merge_command)
            print("------------------------")

    # If there is only one simput file the for loop won't do anything
    if len(simput_files) == 1:
        merged_files.append(simput_files[0])

    # Move the last merged file to the save folder
    final_simput_path = copy_to_save_folder(simput_file_path=merged_files[-1], run_dir=run_dir,
                                            output_file_path=output_file_path,
                                            verbose=verbose)

    return final_simput_path


def combine_simput_back_agn(run_dir, simput_image_path, background=None, agn=None):
    simput_files_combine = [simput_image_path]

    if agn:
        simput_files_combine.append(agn)

    if background:
        simput_files_combine.append(background)

    # in order to not have the names be too long, move the simputs to the rundir with a short name

    for i in range(len(simput_files_combine)):
        org_file = simput_files_combine[i]
        new_file = os.path.join(run_dir, f"{i}.simput")
        shutil.copy(org_file, new_file)
        simput_files_combine[i] = new_file

    final_simput_path = merge_and_save_simputs(run_dir, simput_files=simput_files_combine,
                                               output_file_path=os.path.join(run_dir, "final.simput"))

    return final_simput_path


def get_simputs(simput_base_path, mode, amount=-1, order='normal'):
    # Order options: normal (fron to back), reversed (back to front), random
    # Put the mode and amount in a list if they are passed as single values
    if type(mode) != list:
        mode = [mode]

    if type(amount) != list:
        amount = [amount]
    assert len(mode) == len(
        amount), f"the number of modes ({len(mode)}) should be equal to the number of amounts ({len(amount)})"

    simput_files = []
    counter = 0

    for m, n in zip(mode, amount):
        if amount == 0:
            continue

        simput_path = os.path.join(simput_base_path, m)
        files = os.listdir(simput_path)

        if order == 'normal':
            # Already in normal option
            pass
        elif order == 'reversed':
            files.reverse()
        elif order == 'random':
            random.shuffle(files)
        else:
            raise ValueError(f'Order: {order} not in known options list of "normal", "reversed", "random"')

        # Select only simput files
        for file in files:
            if n != -1:  # -1 means do all
                if counter >= n:  # we reached the desired amount of simputs of this mode
                    break
            if file.endswith(".simput.gz"):
                simput_files.append(os.path.join(simput_path, file))
                counter += 1

    return simput_files
