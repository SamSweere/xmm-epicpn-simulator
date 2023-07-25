# Content

This file describes the individual parameters of each config file.

# illustris.json

- `multiprocessing`: Contains information needed for multiprocessing
    - `gb_per_process`:
    - `num_processes`: How many processes to use. `0` is as many as possible within cpu/memory capabilities.
- `environment`:
    - `working_directory`: Where heasoft is located at and where everything will be saved.
- `illustris`:
    - `api_key`: Your private API key. If you do not have one, you can request it
      at https://www.tng-project.org/users/register/. **Do not upload this to GitHub!**
    - `dataset_dir`: The name of directory the files will be saved in (relative to `working_directory`)
    - `simulation_names`:
        - Options: `['TNG50-1', 'TNG50-2', 'TNG50-3', 'TNG50-4', 'TNG100-1', 'TNG100-2', 'TNG100-3', 'TNG300-1', 'TNG300-2', 'TNG300-3']`
        - Illlustris TNG300 is the largest simulation (300Mpc box) and has 3 resolutions and 3 sub boxes of full
          physics.
    - `snapshot_nums`: The tng simulations have snapshots at different times, 99 is the last snapshot. Multiple numbers
      can be used.
    - `top_n`: Top number of each simulation
    - `emin`
    - `emax`
    - `redshift`: The redshift to the object
    - `overwrite`:
    - `create_preview`: Create preview images
    - `width`: How much to zoom by changing the width of the image
    - `resolution`:
    - `modes`:
    - `proj_normals`:
      - Only relevant for `mode == proj`
      - Take the normals from all directions
    - `slice_axes`:
      - Only relevant for `mode == slice`
      - Normals do not work for slices, thus make them axes