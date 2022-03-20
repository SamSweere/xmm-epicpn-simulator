import datetime
import json
import os


def write_dict_to_file(d, path):
    file = open(path, "w")
    json.dump(d, file)
    file.close()


def load_dict_from_file(path):
    file = open(path, 'r')
    d = json.load(file)
    file.close()
    return d


# Since the cutout files are quite large (+- 1gb) we do not want to download them if we already have them
# in order to do this we check

def get_saved_file(cutout_url, datafolder_name):
    reg_filename = 'reg.json'
    filename_dict = {}  # Start with an empty dict
    reg_path = os.path.join(datafolder_name, reg_filename)
    if not os.path.exists(datafolder_name):
        # If the folder does not exist make it and the registry file
        os.makedirs(datafolder_name)
        write_dict_to_file(filename_dict, reg_path)
        return ""

    # The folder exists, load the filenames
    filename_dict = load_dict_from_file(reg_path)
    if cutout_url in filename_dict:
        return filename_dict[cutout_url]
    else:
        return ""


def save_cutout(cutout_url, content, datafolder_name):
    reg_filename = 'reg.json'
    reg_path = os.path.join(datafolder_name, reg_filename)
    # Create a unique filename, use the url as base if possible
    try:
        parts = cutout_url.split('/')
        filename = parts[-6].replace('-', '_') + "_z_" + parts[-4] + "_" + parts[-3] + "_" + parts[-2] + ".hdf5"
        # url = 'http://www.tng-project.org/api/TNG300-1/snapshots/99/subhalos/538308/cutout.hdf5'
        # becomes: TNG300-1_z_99_subhalos_538308.hdf5
    except:
        # The url is not the right format, use the datetime timestamp as name
        d = datetime.datetime.utcnow()
        epoch = datetime.datetime(1970, 1, 1)
        t = (d - epoch).total_seconds()
        filename = str(t).replace(".", "") + ".hdf5"

    # Write the file
    with open(os.path.join(datafolder_name, filename), 'wb') as f:
        f.write(content)

    # Update the reg dict
    filename_dict = load_dict_from_file(reg_path)
    filename_dict[cutout_url] = filename
    write_dict_to_file(filename_dict, reg_path)

    return filename
