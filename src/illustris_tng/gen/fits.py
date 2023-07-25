import os
from datetime import datetime
from pathlib import Path

import yt


def cutout_to_xray_fits(
        cutout: Path,
        output_dir: Path,
        sub,
        mode='proj',
        emin=0.5,
        emax=2.0,
        normal='y',
        width=1000,
        resolution=2048,
        redshift=0.05,
        overwrite=False,
        create_preview=True
):
    filename = cutout.name.replace('.hdf5', '')

    if isinstance(normal, str):
        normal_str = normal
    elif isinstance(normal, list) and len(normal) == 3:
        normal_str = f"{normal[0]}_{normal[1]}_{normal[2]}"
    else:
        raise ValueError(f"Unexpected `normal`: {normal}")

    filename = f"{filename}_m_{mode}_r_{resolution}_w_{width}_n_{normal_str}"

    fits_filename = filename + '.fits'
    fits_path = output_dir / fits_filename

    if fits_path.exists() and not overwrite:
        raise FileExistsError(f"{fits_path} already exists and `overwrite` is False")

    try:

        # simput_file_name_final = filename + ".simput"  # parts of the simput file will be fits files, this will be the final filename

        print(f"Processing cutout: {cutout}. Arguments: mode={mode}; normal={normal}; emin={emin}; emax={emax}; "
              f"width={width}; resolution={resolution}; redshift={redshift}; overwrite={overwrite}; "
              f"create_preview={create_preview}, output_dir={output_dir}; cutout_datafolder={cutout.parent}")
        ds = yt.load(cutout, default_species_fields="ionized")

        # See documentation (https://yt-project.org/doc/reference/api/yt.fields.xray_emission_fields.html?highlight=xray_emission_fields#module-yt.fields.xray_emission_fields)
        yt.add_xray_emissivity_field(ds, emin, emax, redshift)  # , table_type='apec')

        if mode == "proj":
            yt_fits = yt.FITSOffAxisProjection(ds, normal=normal, fields=('gas', 'xray_photon_intensity_0.5_2.0_keV'),
                                               center=(sub['pos_x'], sub['pos_y'], sub['pos_z']),
                                               width=(width, "code_length"), image_res=resolution)
        elif mode == "slice":
            # TODO: there is a bug with FITSOffAxisSlice: yt.utilities.exceptions.YTNonIndexedDataContainer:
            #  The data container (<class 'yt.data_objects.index_subobjects.particle_container.ParticleContainer'>)
            #  is an unindexed type.  Operations such as ires, icoords, fcoords and fwidth will not work on it.

            # # TODO:  Doing a workaround
            # yt_fits = yt.FITSOffAxisSlice(ds, normal=normal, fields=('gas', 'xray_photon_intensity_0.5_2.0_keV'),
            #                                    center=(sub['pos_x'], sub['pos_y'], sub['pos_z']),
            #                                    width=(width, "code_length"), image_res=resolution)
            yt_fits = yt.FITSSlice(ds, axis=normal, fields=('gas', 'xray_photon_intensity_0.5_2.0_keV'),
                                   center=(sub['pos_x'], sub['pos_y'], sub['pos_z']),
                                   width=(width, "code_length"), image_res=resolution)
        else:
            raise ValueError(f"Mode {mode} is not suported. Suported modes: `projection` and 'slice'")

        yt_fits.update_header(field="all", key="NORMAL", value=str(normal))
        yt_fits.update_header(field="all", key="WIDTH", value=width)
        yt_fits.update_header(field="all", key="REDSHIFT", value=redshift)
        yt_fits.update_header(field="all", key="EMIN", value=emin)
        yt_fits.update_header(field="all", key="EMAX", value=emax)

        yt_fits.writeto(str(fits_path.resolve()), overwrite=overwrite)

        # if show_results:
        #     slc = yt.ParticlePlot(ds, "particle_position_x", "particle_position_y", "particle_mass",
        #                           center=(sub['pos_x'], sub['pos_y'], sub['pos_z']), width=(1000.0, 'code_length'))
        #     slc.set_unit('particle_mass', 'Msun')
        #     slc.show()

        # sp = ds.sphere((sub['pos_x'], sub['pos_y'], sub['pos_z']), (1000.0 / zoom, "code_length"))
        #
        # # Make a thermal source model from the temperature, metallicity, density using an APEC emission model, only using particles where kT>0.1keV and assuming 0.5-7keV energy band
        # thermal_model = pyxsim.ThermalSourceModel(
        #     "apec",
        #     emin=emin,
        #     emax=emax,
        #     nchan=int((emax - emin) * 1000),
        #     # The number of channels in the spectrum. If one is thermally broadening lines (the default), it is recommended that this number create an energy resolution per channel of roughly 1 eV.
        #     temperature_field=("PartType0", "Temperature")
        # )
        #
        # # Make a photon list from the thermal source model, assuming z=0.05 and NH = 4e20[cm-2]. Area and spectral response are taken from ARF/RMF files of e.g. XMM
        # # TODO: is the redshift interesting?
        # redshift = 0.05  # The redshift to the object.
        # area = (600., "cm**2")  # A constant collecting effective area to generate the photons with.
        # exp_time = (exp_time, "ks")  # The exposure time to generate the photons with.
        # center = sp.center  # A center in 3D for the photon positions. If not specified,
        # # the center of the `data_source` will be chosen.
        #
        # photons = pyxsim.PhotonList.from_data_source(sp, redshift, area, exp_time,
        #                                              thermal_model, center=center)
        #
        # # Next to create the events we project the photons through the foreground galactic absorption along the y axis, assuming a sky center (0,0)
        # nH = 0.04  # [10^22 cm-2]
        # sky_center = (0, 0)  # sky center in RA,DEC
        # # normal (character or array-like) – Normal vector to the plane of projection. If “x”, “y”, or “z”, will assume to be along that axis (and will probably be faster). Otherwise, should be an off-axis normal vector, e.g [1.0, 2.0, -3.0]
        # events = photons.project_photons(normal, sky_center, absorb_model="wabs", nH=nH)
        #
        # # We need to save the event file as a SIMPUT file
        # fov = (30 / zoom, "arcmin")  # the field of view / width of the image
        # nx = resolution  # The resolution of the image on a side, this is a bit larger than the expected resolution of 4x
        #
        #
        # events.write_fits_image(fits_path, fov, nx, overwrite=overwrite, emin=emin, emax=emax)
        #
        # # events.write_simput_file(filename, overwrite=True, emin=emin, emax=emax)
        #
        # # Rename the filename to end on .simput
        # # os.rename(filename + "_simput.fits", simput_file_name_final)

        print('')
        print(f'Converted {cutout} to {fits_filename}')

        if create_preview:
            fig = yt_fits.to_aplpy()

            # Create a png from the fits e
            preview_dir = os.path.join(output_dir, 'preview')
            if not os.path.exists(preview_dir):
                os.makedirs(preview_dir)

            # fig = aplpy.FITSFigure(fits_path)
            fig.show_colorscale(stretch='log', cmap="arbre")
            # fig.add_colorbar() # This can cause an error
            fig.set_title(f"log strechted, {filename}")
            preview_filepath = os.path.join(preview_dir, filename + '.png')
            fig.save(os.path.join(preview_dir, preview_filepath))
            fig.close()

            print("Saved preview to:", preview_filepath)

    except Exception as e:

        # Returns a datetime object containing the local date and time
        dateTimeObj = datetime.now()
        message = str(dateTimeObj) + f" ERROR, failed to process {filename} with errro: {e}"

        with open(os.path.join(output_dir, "illustristng_failed_processes.txt"), "a") as myfile:
            myfile.write(message + "\n")
