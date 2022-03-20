import requests
from tqdm import tqdm

from illustris_tng.data_handling import get_saved_file, save_cutout


def get(path, headers, cutout_datafolder, params=None):
    if 'cutout' in path:
        # If the path is a cutout getter, check if we already have it
        filename = get_saved_file(path, cutout_datafolder)
        if filename != "":
            # We have the file return the name
            return filename

    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json()  # parse json responses automatically

    if 'cutout' in path and 'content-disposition' in r.headers:
        filename = save_cutout(cutout_url=path, content=r.content, datafolder_name=cutout_datafolder)
        return filename  # return the filename string
    return r


def get_and_save_cutouts(subs, headers, cutout_datafolder):
    cutout_request = {
        'gas': 'Coordinates,Density,ElectronAbundance,GFM_Metallicity,InternalEnergy,Masses,NeutralHydrogenAbundance,Velocities'}

    for i in tqdm(range(len(subs))):
        try:
            sub = get(subs[i]['url'], headers=headers, cutout_datafolder=cutout_datafolder)
            cutout = get(sub['cutouts']['subhalo'], headers=headers, cutout_datafolder=cutout_datafolder,
                         params=cutout_request)
        except Exception as e:
            print(f"Failed to load number {i}, error:")
            print(e)


def get_cutout_from_sub(sub, headers, cutout_datafolder):
    print("Downloading cutout:", sub['cutouts']['subhalo'])
    cutout_request = {
        'gas': 'Coordinates,Density,ElectronAbundance,GFM_Metallicity,InternalEnergy,Masses,NeutralHydrogenAbundance,Velocities'}
    cutout = get(sub['cutouts']['subhalo'], headers=headers, cutout_datafolder=cutout_datafolder, params=cutout_request)
    return cutout
