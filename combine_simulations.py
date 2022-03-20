# This file should not depend on other in project dependencies
import os
import random
import shutil
from datetime import datetime

from astropy.io import fits

config = {
    'docker': False,
    'home': os.path.join(os.path.expanduser("~"), 'Documents/ESA/data/sim'),
    # The home directory when not running inside docker
    'debug': False,
    'verbose': True,
    'run_dir': None,
    'data_dir': None,
    'dataset_dir': None,
    'dataset_name': 'xmm_demo_dataset',

    'mode': ['img', 'agn', 'background', 'test_grid', 'exposure_map'],
    'num': [-1, -1, -1, -1, 0],  # -1 is everything
    'agn': [True, False, False, False, False],  # Add agn
    'background': [True, True, False, True, False],  # Add backgound
    'det_mask': [True, True, True, True, True],  # Apply detector mask
    'sample_n': [1, 1, 0, 1, 0],  # Number of sampled agns and backgrounds
    'res_mult': [1, 2],
    'exposure': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
}

if config['docker']:
    config['data_dir'] = '/home/heasoft/data/sim'
    config['run_dir'] = '/home/heasoft/tmp_comb'

else:
    # Running on laptop/pc
    config['data_dir'] = config['home']
    config['run_dir'] = os.path.join(config['data_dir'], 'tmp_comb')

config['dataset_dir'] = os.path.join(config['data_dir'], config['dataset_name'])


def load_simulation_file(file_path, run_dir="."):
    if file_path.endswith('.fits') or file_path.endswith('.fits.gz'):
        hdu = fits.open(file_path)
    else:
        raise AssertionError(f"File {file_path} is not a .fits file or a .fits.gz file")

    return hdu


def process_one_simulation(file_path, data, general_header_dict, run_dir='.'):
    hdu = load_simulation_file(file_path, run_dir)
    hdu_data = hdu['PRIMARY'].data
    header = hdu['PRIMARY'].header

    hdu_general_header_dict = {
        'TELESCOP': (header['TELESCOP'], header.comments['TELESCOP']),
        'WCSAXES': (header['WCSAXES'], header.comments['WCSAXES']),
        'CRPIX1': (header['CRPIX1'], header.comments['CRPIX1']),
        'CRPIX2': (header['CRPIX2'], header.comments['CRPIX2']),
        'CDELT1': (header['CDELT1'], header.comments['CDELT1']),
        'CDELT2': (header['CDELT2'], header.comments['CDELT2']),
        'CUNIT1': (header['CUNIT1'], header.comments['CUNIT1']),
        'CUNIT2': (header['CUNIT2'], header.comments['CUNIT2']),
        'CTYPE1': (header['CTYPE1'], header.comments['CTYPE1']),
        'CTYPE2': (header['CTYPE2'], header.comments['CTYPE2']),
        'CRVAL1': (header['CRVAL1'], header.comments['CRVAL1']),
        'CRVAL2': (header['CRVAL2'], header.comments['CRVAL2']),
        'LONPOLE': (header['LONPOLE'], header.comments['LONPOLE']),
        'RADESYS': (header['RADESYS'], header.comments['RADESYS']),
        'EXPOSURE': (header['EXPOSURE'], header.comments['EXPOSURE']),
        'RESMULT': (header['RESMULT'], header.comments['RESMULT'])
    }

    if data is None:
        data = hdu_data
    else:
        data += hdu_data

    if general_header_dict is None:
        general_header_dict = hdu_general_header_dict
    else:
        if general_header_dict != hdu_general_header_dict:
            raise ValueError("General header values are not the same, something is wrong")

    return data, header, general_header_dict


def apply_det_mask(data, det_mask_path):
    det_mask = fits.open(det_mask_path)[0].data
    data = data * det_mask

    return data


def combine_simulations(img_path=None, agn_path=None, background_path=None, det_mask_path=None, run_dir='.'):
    if not img_path and not agn_path and not background_path:
        raise AssertionError("Not all inputfiles can be None")

    header = fits.Header()

    data = None
    general_header_dict = None

    if img_path is not None:
        data, img_header, general_header_dict = process_one_simulation(img_path, data, general_header_dict, run_dir)

        header["IMG_FILE"] = (os.path.basename(img_path), "Input simulation file used as img")
        header["IMG_SIMP"] = (img_header['SIMPUT'], "Input img simput used for simulation")

    if agn_path is not None:
        data, agn_header, general_header_dict = process_one_simulation(agn_path, data, general_header_dict, run_dir)

        header["AGN_FILE"] = (os.path.basename(agn_path), "Input simulation file used as agn")
        header["AGN_SIMP"] = (agn_header['SIMPUT'], "Input agn simput used for simulation")

    if background_path is not None:
        data, background_header, general_header_dict = process_one_simulation(background_path, data,
                                                                              general_header_dict, run_dir)

        header["BCK_FILE"] = (os.path.basename(background_path), "Input simulation file used as background")
        header["BCK_SIMP"] = (background_header['SIMPUT'], "Input background simput used for simulation")

    if det_mask_path is not None:
        # Apply the detector mask
        data = apply_det_mask(data, det_mask_path)
        header["DETMASK"] = (os.path.basename(det_mask_path), "Detector mask used on this image")

    # Add the general background to the header
    for key, value in general_header_dict.items():
        header[key] = value

    header['COMMENT'] = f"Created by Sam Sweere (samsweere@gmail.com) for ESAC at " \
                        f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}"

    hdu = fits.PrimaryHDU(data, header=header)
    return hdu


def random_sample_sim_files(sample_dir, num):
    files = []
    for file in os.listdir(sample_dir):
        if file.endswith(".fits") or file.endswith(".fits.gz"):
            files.append(file)

    if len(files) < num:
        print(f"WARNING: sampling {num} files from a source of {len(files)} files. This will cause (many) duplicates")

    sampled_files = random.sample(files, num)

    sampled_file_paths = []
    for file in sampled_files:
        sampled_file_paths.append(os.path.join(sample_dir, file))

    return sampled_file_paths


def combine_and_save_sim(output_dir, img_path=None, agn_path=None, background_path=None, det_mask_path=None,
                         run_dir='.'):
    print(f"Combining: {img_path}, {agn_path} and {background_path}")

    if not os.path.exists(output_dir):
        raise NotADirectoryError(f"Output directory {output_dir} does not exist")

    hdu = combine_simulations(img_path=img_path, agn_path=agn_path, background_path=background_path,
                              det_mask_path=det_mask_path, run_dir=run_dir)

    # Determine the output name
    out_name = None
    if img_path is not None:
        out_name = os.path.basename(img_path).replace(".fits.gz", "").replace(".fits", "")
    if agn_path is not None:
        agn_name = os.path.basename(agn_path).replace(".fits.gz", "").replace(".fits", "")
        if out_name is not None:
            out_name += "_agn_" + agn_name.split("_")[1]  # "_agn"
        else:
            out_name = agn_name
    if background_path is not None:
        background_name = os.path.basename(background_path).replace(".fits.gz", "").replace(".fits", "")
        if out_name is not None:
            out_name += "_back_" + background_name.split("_")[-1]
        else:
            out_name = background_name
    if out_name is None:
        raise AssertionError("Not all inputfiles can be None")

    compressed_out_path = os.path.join(output_dir, out_name + ".fits.gz")

    hdu.writeto(compressed_out_path, overwrite=True)

    # compressed_out_path = os.path.join(output_dir, out_name + ".fits.gz")
    # compress_file(in_file_path=tmp_out_path, out_file_path=compressed_out_path)

    # # Remove the tmp fits file
    # os.remove(tmp_out_path)


def change_mult_path(filepath, base_res_mult, res_mult):
    if type(filepath) == list:
        new_filepaths = []
        for filename in filepath:
            filename = filename.replace(f'/{base_res_mult}x/', f'/{res_mult}x/') \
                .replace(f'_mult_{base_res_mult}_', f'_mult_{res_mult}_')
            new_filepaths.append(filename)
        return new_filepaths
    elif type(filepath) == str:
        filepath = filepath.replace(f'/{base_res_mult}x/', f'/{res_mult}x/') \
            .replace(f'_mult_{base_res_mult}_', f'_mult_{res_mult}_')
        return filepath
    else:
        raise TypeError(f"{filepath} does not have type list or string but {type(filepath)}")


# output_dir = '/home/sam/Documents/ESA/data/sim/dev_dataset/combine_test'
# img_path = '/home/sam/Documents/ESA/data/sim/dev_dataset/50ks/img/1x/TNG300_1_z_99_subhalos_0_m_slice_r_2048_w_1000_n_x_p0_0_5ev_p1_2_0ev_sb_50_zoom_1_offx_0_0_offy_0_0_mult_1_50ks.tar.gz'
# agn_path = '/home/sam/Documents/ESA/data/sim/dev_dataset/50ks/agn/1x/agn_16263641327751281_p0_0.5ev_p1_2.0ev_mult_1_50ks.tar.gz'
# background_path = '/home/sam/Documents/ESA/data/sim/dev_dataset/50ks/background/1x/background_mult_1_50ks_16263951074803240.tar.gz'


# Create the rundir if it does not exits
run_dir = config['run_dir']
if not os.path.exists(run_dir):
    os.makedirs(run_dir)

for exposure in config['exposure']:
    exp_path_part = str(exposure) + "ks"
    print(f"Processing exposure: {exp_path_part}")
    for mode, mode_num, agn, background, sample_n, det_mask in zip(config['mode'], config['num'], config['agn'],
                                                                   config['background'], config['sample_n'],
                                                                   config['det_mask']):
        if mode_num == 0:
            continue  # Do not sample_num from this mode

        mode_path_part = mode
        print(f"Processing mode: {mode}")

        # Get the input images from the lowest res mult given
        sample_res_mult = config['res_mult'][0]
        sample_res_mult_part = str(sample_res_mult) + 'x'
        in_dir = os.path.join(config['dataset_dir'], exp_path_part, mode_path_part, sample_res_mult_part)
        img_in = []
        for file in os.listdir(in_dir):
            if file.endswith(".fits") or file.endswith(".fits.gz"):
                img_in.append(os.path.join(in_dir, file))

        if mode_num != -1:
            img_in = img_in[:mode_num]  # Limit the amount

        for img_path in img_in:
            # Start with no agn and background files, this can be replaced when sampling. Always have one in order to process no agn and no background
            agn_files = [None] * max(sample_n, 1)
            background_files = [None] * max(sample_n, 1)

            # Sample from the agns
            if agn:
                agn_dir = os.path.join(config['dataset_dir'], exp_path_part, 'agn', sample_res_mult_part)
                agn_files = random_sample_sim_files(sample_dir=agn_dir, num=sample_n)

            for res_mult in config['res_mult']:
                res_mult_path_part = str(res_mult) + "x"
                print(f"Processing res mult: {res_mult}")

                # Get the detector mask path, these should be present in the simulated dataset
                if det_mask:
                    det_mask_path = os.path.join(config['dataset_dir'], 'detector_mask', f'{res_mult}x',
                                                 f'pn_mask_500_2000_detxy_{res_mult}x.ds')
                else:
                    det_mask_path = None  # Apply no det mask

                    # Create the output dir:
                output_dir = os.path.join(config['dataset_dir'], 'combined', exp_path_part, mode_path_part,
                                          res_mult_path_part)
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)

                if mode != 'background':
                    # Get the img, agn and background files corresponding to this resolution
                    img_path_mult = change_mult_path(img_path, sample_res_mult, res_mult)
                else:
                    # If the image mode is background we need to just sample_num a background
                    background_dir = os.path.join(config['dataset_dir'], exp_path_part, 'background',
                                                  res_mult_path_part)
                    img_path_mult = random_sample_sim_files(sample_dir=background_dir, num=1)[0]

                agn_files_mult = [None] * max(sample_n, 1)
                if agn:
                    agn_files_mult = change_mult_path(agn_files, sample_res_mult, res_mult)

                # Since the background is the same input for everything we cannot match it
                if background:
                    background_dir = os.path.join(config['dataset_dir'], exp_path_part, 'background',
                                                  res_mult_path_part)
                    background_files = random_sample_sim_files(sample_dir=background_dir, num=sample_n)

                for agn_path, background_path in zip(agn_files_mult, background_files):
                    try:
                        combine_and_save_sim(output_dir=output_dir, img_path=img_path_mult, agn_path=agn_path,
                                             background_path=background_path, det_mask_path=det_mask_path,
                                             run_dir=run_dir)
                    except Exception as e:
                        print("Error:", e)

if not config['debug']:
    # Remove the temporary run dir
    shutil.rmtree(config['run_dir'])
