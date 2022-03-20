import os
import shutil

import numpy as np
from shutil import copy2

def get_fits_files(dataset_dir):
    if not os.path.exists(dataset_dir):
        raise OSError(f"Dataset directory {dataset_dir} does not exists")

    fits_files = []
    file_names = []

    # Only save the files that end with .fits
    for file in os.listdir(dataset_dir):
        if file.endswith(".fits") or file.endswith(".fits.gz"):
            fits_files.append(file)
            file_name = os.path.splitext(file)[0]
            file_names.append(file_name)

    print(f"Detected {len(fits_files)} fits files in {dataset_dir}")

    return fits_files, file_names

def match_file_list(filelist1, filelist2, split_key):
    """

    Filters and groups two filelists on having the same base name when split with the split_key.

    Args:
        filelist1: List of files1
        filelist2: List of files2
        split_key: The string where the filenames should be split on

    Returns:
        A list of lists of filelist1, a list of lists of filelist2 and list of the base_filenames all grouped on the
        base filenames.

    """

    base_names1 = [x.split(split_key)[0] for x in filelist1]
    base_names2 = [x.split(split_key)[0] for x in filelist2]

    # Get the intersection between the lists
    name_intersect_set = set(base_names1) & set(base_names2)
    name_intersect = list(name_intersect_set)

    filtered_filelist1 = [[] for i in range(len(name_intersect))] # Note that every list has to be created separately
    # to be separate items
    filtered_filelist2 = [[] for i in range(len(name_intersect))]

    # Only keep the ones that the files in both resolutions
    # Since we can have multiple files with the same base_name we want to cluster these too
    for file1 in filelist1:
        base_name = file1.split(split_key)[0]
        if base_name in name_intersect_set:  # Make use of the set speed
            index = name_intersect.index(base_name)
            filtered_filelist1[index].append(file1)

    for file2 in filelist2:
        base_name = file2.split(split_key)[0]
        if base_name in name_intersect_set:  # Make use of the set speed
            index = name_intersect.index(base_name)
            filtered_filelist2[index].append(file2)

    return filtered_filelist1, filtered_filelist2, name_intersect


source_path = '/home/sam/Documents/ESA/data/sim/xmm_sim_dataset'
dest_path = '/home/sam/Documents/ESA/data/sim/xmm_dev_big_dataset'

if not os.path.exists(dest_path):
    os.makedirs(dest_path)

# Copy the detector mask files
shutil.copytree(os.path.join(source_path, "detector_mask"), os.path.join(dest_path, "detector_mask"))

copy_num = 1000
jump_size = 10

base_file_names = None

dataset_dir = source_path

for exp_n in np.arange(10,110,10):
    exp = str(exp_n)
    for mode in ['img', 'agn', 'background', 'test_grid']:
        rel_lr_img_dir = os.path.join(str(exp) + 'ks', mode, str(1) + 'x')
        rel_hr_img_dir = os.path.join(str(exp) + 'ks', mode, str(2) + 'x')
        lr_img_source_path = os.path.join(dataset_dir, rel_lr_img_dir)
        hr_img_source_path = os.path.join(dataset_dir, rel_hr_img_dir)
        lr_img_dest_path = os.path.join(dest_path, rel_lr_img_dir)
        hr_img_dest_path = os.path.join(dest_path, rel_hr_img_dir)
        os.makedirs(lr_img_dest_path)
        os.makedirs(hr_img_dest_path)


        lr_img_files, _ = get_fits_files(dataset_dir=lr_img_source_path)
        hr_img_files, _ = get_fits_files(dataset_dir=hr_img_source_path)

        if mode == "background":
            # Just copy some they do not match
            for lr_img_name in lr_img_files[:copy_num]:
                copy2(os.path.join(lr_img_source_path, lr_img_name), os.path.join(lr_img_dest_path, lr_img_name))

            for hr_img_name in hr_img_files[:copy_num]:
                copy2(os.path.join(hr_img_source_path, hr_img_name), os.path.join(hr_img_dest_path, hr_img_name))

            continue


        # Filter the files such that only files that have a match in both the resolution are present
        # The lr_img_files and hr_img_files will be a list of list containing all the files with the same base_img name
        lr_img_files, hr_img_files, base_img_files = match_file_list(lr_img_files,
                                                                                    hr_img_files,
                                                                                    split_key="_mult_")

        base_img_files.sort()
        lr_img_files.sort()
        hr_img_files.sort()

        if len(base_img_files) < copy_num * jump_size:
            # Copy the first copy_num
            base_img_files = base_img_files[:copy_num]
        else:
            # Copy every tenth in order to prevent same samples
            tmp_fnames = []
            for i in np.arange(0, copy_num * jump_size, jump_size):
                f_name = base_img_files[i]
                tmp_fnames.append(f_name)
            base_img_files = tmp_fnames

        for f_name in base_img_files:
            for lr_img_name in lr_img_files:
                lr_img_name = lr_img_name[0]
                if lr_img_name.split("_mult_")[0] == f_name:
                    copy2(os.path.join(lr_img_source_path, lr_img_name), os.path.join(lr_img_dest_path, lr_img_name))

            for hr_img_name in hr_img_files:
                hr_img_name = hr_img_name[0]
                if hr_img_name.split("_mult_")[0] == f_name:
                    copy2(os.path.join(hr_img_source_path, hr_img_name), os.path.join(hr_img_dest_path, hr_img_name))




#
#             for base_name in base_file_names:
#                 for f_name in f_names:
#                     print("base name:", base_name)
#                     print("f_name", f_name)
#                     base_high_name = f_name.split("_mult_")[0]
#                     if base_high_name == base_name:
#                         # Same sim file copy
#                         copy2(os.path.join(root, f_name), os.path.join(dest_root, f_name))
#                         print("Copy", base_high_name)
#                         break
#
#             base_file_names.append(f_name.split("_mult_")[0])
#             copy2(os.path.join(root, f_name), os.path.join(dest_root, f_name))
#
#
#
#
#
# # Get all the image directories
# img_dir = os.path.join(self.dataset_dir, str(self.lr_exp) + 'ks', self.mode,
#                                str(self.lr_res_mult) + 'x')
# self.hr_img_source_path = os.path.join(self.dataset_dir, str(self.hr_exp) + 'ks', self.mode,
#                                        str(self.hr_res_mult) + 'x')
# # Get the fits files and file names
# self.lr_img_files, _ = get_fits_files(dataset_dir=self.lr_img_source_path)
# self.hr_img_files, _ = get_fits_files(dataset_dir=self.hr_img_source_path)
#
# # Filter the files such that only files that have a match in both the resolution are present
# # The lr_img_files and hr_img_files will be a list of list containing all the files with the same base_img name
# self.lr_img_files, self.hr_img_files, self.base_img_files = match_file_list(self.lr_img_files,
#                                                                             self.hr_img_files,
#                                                                             split_key="_mult_")
# print(f"Found {len(self.base_img_files)} image pairs (lr and hr simulation matches)")
#
# if self.agn:
#     self.lr_agn_dir = os.path.join(self.dataset_dir, str(self.lr_exp) + 'ks', 'agn',
#                                    str(self.lr_res_mult) + 'x')
#     self.hr_agn_dir = os.path.join(self.dataset_dir, str(self.hr_exp) + 'ks', 'agn',
#                                    str(self.hr_res_mult) + 'x')
#
#     self.lr_agn_files, _ = get_fits_files(dataset_dir=self.lr_agn_dir)
#     self.hr_agn_files, _ = get_fits_files(dataset_dir=self.hr_agn_dir)
#
#     # Filter the files such that only files that have a match in both the resolution are present
#     self.lr_agn_files, self.hr_agn_files, self.base_agn_files = match_file_list(self.lr_agn_files,
#                                                                                 self.hr_agn_files,
#                                                                                 split_key="_mult_")
#     print(f"Found {len(self.base_agn_files)} agn image pairs (lr and hr simulation matches)")
#
# if self.lr_background:
#     self.lr_background_dir = os.path.join(self.dataset_dir, str(self.lr_exp) + 'ks', 'background',
#                                           str(self.lr_res_mult) + 'x')
#
#     self.lr_background_files, _ = get_fits_files(dataset_dir=self.lr_background_dir)
#
# if self.hr_background:
#     self.hr_background_dir = os.path.join(self.dataset_dir, str(self.hr_exp) + 'ks', 'background',
#                                           str(self.hr_res_mult) + 'x')
#
#     self.hr_background_files, _ = get_fits_files(dataset_dir=self.hr_background_dir)




# for root, d_names, f_names in os.walk(source_path):
#     dest_root = root.replace(source_path, dest_path)
#     for d_name in d_names:
#         print(os.path.join(dest_root, d_name))
#         os.makedirs(os.path.join(dest_root, d_name))
#
#         # Sort the f_names to speed up the search
#         f_names.sort()
#
#         last_d_name = os.path.split(d_name)[-1]
#         print("Last_d_name", last_d_name)
#
#         if (last_d_name == "1x" or last_d_name == "2x" or last_d_name == "4x") and (len(f_names) != 0):
#             if base_file_names is None:
#                 # Save the base names
#                 base_file_names = []
#
#                 if len(f_names) < copy_num * 100:
#                     # Copy the first copy_num
#                     f_names = f_names[:copy_num]
#                 else:
#                     # Copy every tenth in order to prevent same samples
#                     tmp_fnames = []
#                     for i in np.arange(0, copy_num * 100, 100):
#                         f_name = f_names[i]
#                         tmp_fnames.append(f_name)
#                     f_names = tmp_fnames
#
#                 for f_name in f_names:
#                     base_file_names.append(f_name.split("_mult_")[0])
#                     copy2(os.path.join(root, f_name), os.path.join(dest_root, f_name))
#             else:
#                 # Find the base names of the 1x
#                 for base_name in base_file_names:
#                     for f_name in f_names:
#                         print("base name:", base_name)
#                         print("f_name", f_name)
#                         base_high_name = f_name.split("_mult_")[0]
#                         if base_high_name == base_name:
#                             # Same sim file copy
#                             copy2(os.path.join(root, f_name), os.path.join(dest_root, f_name))
#                             print("Copy", base_high_name)
#                             break
#         else:
#             # Reset base names
#             base_file_names = None