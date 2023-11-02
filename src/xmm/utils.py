import os
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, Literal, Tuple

import heasoftpy as hsp
import numpy as np
from astropy.io import fits
from loguru import logger
from pysas.wrapper import Wrapper as sas

available_instruments = ["epn", "emos1", "emos2"]
instrument_to_sixte_dir = {
    "epn": "epicpn",
    "emos1": "epicmos",
    "emos2": "epicmos"
}


def get_fov_for_instrument(
        instrument_name: Literal["epn", "emos1", "emos2"]
) -> float:
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    if instrument_name == "epn":
        from src.xmm.epn import get_fov
        return get_fov()

    if instrument_name == "emos1":
        from src.xmm.emos import get_fov
        return get_fov(1)

    if instrument_name == "emos2":
        from src.xmm.emos import get_fov
        return get_fov(2)


def get_cdelt_for_instrument(
        instrument_name: Literal["epn", "emos1", "emos2"],
        res_mult: int
) -> float:
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    if instrument_name == "epn":
        from src.xmm.epn import get_cdelt
        return get_cdelt(res_mult)

    if instrument_name == "emos1":
        from src.xmm.emos import get_cdelt
        return get_cdelt(1, res_mult)

    if instrument_name == "emos2":
        from src.xmm.emos import get_cdelt
        return get_cdelt(2, res_mult)


def get_pixel_size_for_instrument(
        instrument_name: Literal["epn", "emos1", "emos2"],
        res_mult: int
) -> float:
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    if instrument_name == "epn":
        from src.xmm.epn import get_pixel_size
        return get_pixel_size(res_mult)

    if instrument_name == "emos1":
        from src.xmm.emos import get_pixel_size
        return get_pixel_size(1, res_mult)

    if instrument_name == "emos2":
        from src.xmm.emos import get_pixel_size
        return get_pixel_size(2, res_mult)


def get_surface_for_instrument(
        instrument_name: Literal["epn", "emos1", "emos2"],
        res_mult: int
) -> float:
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    if instrument_name == "epn":
        from src.xmm.epn import get_surface
        return get_surface(res_mult=res_mult)

    if instrument_name == "emos1":
        from src.xmm.emos import get_surface
        return get_surface(emos_num=1, res_mult=res_mult)

    if instrument_name == "emos2":
        from src.xmm.emos import get_surface
        return get_surface(emos_num=2, res_mult=res_mult)


def get_width_height_for_instrument(
        instrument_name: Literal["epn", "emos1", "emos2"],
        res_mult: int
) -> Tuple[int, int]:
    """
    Returns:
        Tuple[int, int]: The width and height of the instrument in pixels.
    """
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    if instrument_name == "epn":
        from src.xmm.epn import get_img_width_height
        return get_img_width_height(res_mult=res_mult)

    if instrument_name == "emos1":
        from src.xmm.emos import get_img_width_height
        return get_img_width_height(emos_num=1, res_mult=res_mult)

    if instrument_name == "emos2":
        from src.xmm.emos import get_img_width_height
        return get_img_width_height(emos_num=2, res_mult=res_mult)


def get_crpix12_for_instrument(
        instrument_name: Literal["epn", "emos1", "emos2"],
        res_mult: int
) -> Tuple[float, float]:
    """
    Returns:
        Tuple[float, float]: The pixels corresponding to the focal point on the instrument.
    """
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    if instrument_name == "epn":
        from src.xmm.epn import get_crpix
        return get_crpix(res_mult=res_mult)

    if instrument_name == "emos1":
        from src.xmm.emos import get_crpix
        return get_crpix(emos_num=1, res_mult=res_mult)

    if instrument_name == "emos2":
        from src.xmm.emos import get_crpix
        return get_crpix(emos_num=2, res_mult=res_mult)


def get_focal_length(
        instrument_name: Literal["epn", "emos1", "emos2"]
) -> float:
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    if instrument_name == "epn":
        from src.xmm.epn import get_focal_length
        return get_focal_length()

    if instrument_name == "emos1":
        from src.xmm.emos import get_focal_length
        return get_focal_length(emos_num=1)

    if instrument_name == "emos2":
        from src.xmm.emos import get_focal_length
        return get_focal_length(emos_num=2)


def get_instrument_files(
        instrument_name: Literal["epn", "emos1", "emos2"]
) -> Path:
    """
    Returns:
        Path: Path to the SIXTE instruments directory.
    Raises:
        NotADirectoryError: If the instrument files have not been downloaded from
            https://www.sternwarte.uni-erlangen.de/sixte/instruments/
    """
    p = Path(os.environ["SIXTE"]) / "share" / "sixte" / "instruments" / "xmm" / instrument_to_sixte_dir[instrument_name]

    if not p.exists():
        raise NotADirectoryError(f"It looks like you haven't downloaded the instrument files provided by SIXTE "
                                 f"(see https://www.sternwarte.uni-erlangen.de/sixte/instruments/)! Please download "
                                 f"and extract them as given in their instructions.")

    return p


def create_xml_files(
        xml_dir: Path,
        instrument_name: Literal["epn", "emos1", "emos2"],
        res_mult: int,
        xmm_filter: Literal["thin", "med", "thick"],
        sim_separate_ccds: bool,
        wait_time: float,
) -> List[Path]:
    """
    Returns:
        List[Path]: A list containing paths to the corresponding CCDs. If sim_separate_ccds == True, then the list
            contains a single Path to /path/to/instrument/combined.xml
    """
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    instrument_files = get_instrument_files(instrument_name=instrument_name)
    instrument_dir = xml_dir / instrument_name
    out_dir = instrument_dir / xmm_filter / f"{res_mult}x"
    instrument_dir.mkdir(exist_ok=True, parents=True)
    out_dir.mkdir(exist_ok=True, parents=True)

    psf_name = get_psf_name(instrument_name, res_mult)
    psf_file = instrument_dir / psf_name

    if instrument_name == "epn":
        from src.xmm.xml import create_pn_xml
        # Add symbolic links to used instrument_files
        files = ["xmm_pn_vignet.fits", f"pn-{xmm_filter}-10.rmf", f"pn-{xmm_filter}-10.arf"]
        for file in files:
            tmp_link = out_dir / file
            tmp_link.symlink_to(instrument_files / file)
        # Add symbolic link to psf_file
        (out_dir / psf_name).symlink_to(psf_file)

        return create_pn_xml(out_dir=out_dir, res_mult=res_mult, xmm_filter=xmm_filter,
                             sim_separate_ccds=sim_separate_ccds,
                             wait_time=wait_time)

    if instrument_name == "emos1" or instrument_name == "emos2":
        from src.xmm.xml import create_mos_xml
        # Add symbolic links to used instrument_files
        files = [f"mos1-{xmm_filter}-10.rmf", f"mos1-{xmm_filter}-10.arf"]
        for file in files:
            tmp_link = out_dir / file
            tmp_link.symlink_to(instrument_files / file)
        # Add symbolic link to psf_file
        (out_dir / psf_name).symlink_to(psf_file)
        return create_mos_xml(out_dir=out_dir, emos_num=1, res_mult=res_mult, xmm_filter=xmm_filter,
                              sim_separate_ccds=sim_separate_ccds, wait_time=wait_time)

    if instrument_name == "emos2":
        from src.xmm.xml import create_mos_xml
        # Add symbolic links to used instrument_files
        files = [f"mos2-{xmm_filter}-10.rmf", f"mos2-{xmm_filter}-10.arf"]
        for file in files:
            tmp_link = out_dir / file
            tmp_link.symlink_to(instrument_files / file)
        # Add symbolic link to psf_file
        (out_dir / psf_name).symlink_to(psf_file)
        return create_mos_xml(out_dir=out_dir, emos_num=2, res_mult=res_mult, xmm_filter=xmm_filter,
                              sim_separate_ccds=sim_separate_ccds, wait_time=wait_time)


def get_xml_files(
        instrument_name: Literal["epn", "emos1", "emos2"],
        res_mult: int,
        xmm_filter: Literal["thin", "med", "thick"],
        sim_separate_ccds: bool
) -> List[Path]:
    """
    Returns:
        List[Path]: A list containing paths to the corresponding CCDs. If sim_separate_ccds == True, then the list
            contains a single Path to /path/to/instrument/combined.xml
    """
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    if instrument_name == "epn":
        from src.xmm.xml import get_pn_xml
        return get_pn_xml(res_mult=res_mult, xmm_filter=xmm_filter, sim_separate_ccds=sim_separate_ccds)

    if instrument_name == "emos1":
        from src.xmm.xml import get_mos_xml
        return get_mos_xml(emos_num=1, res_mult=res_mult, xmm_filter=xmm_filter, sim_separate_ccds=sim_separate_ccds)

    if instrument_name == "emos2":
        from src.xmm.xml import get_mos_xml
        return get_mos_xml(emos_num=2, res_mult=res_mult, xmm_filter=xmm_filter, sim_separate_ccds=sim_separate_ccds)


def create_psf_file(
        out: Path,
        instrument_name: Literal["epn", "emos1", "emos2"],
        res_mult: int,
        energy_list: List[int],
        theta_list: List[int]
) -> None:
    inst_name_map = {
        "epn": "EPN",
        "emos1": "EMOS1",
        "emos2": "EMOS2"
    }
    focal_length = get_focal_length(instrument_name)
    stretch_factor = 1.0 / res_mult

    cdelt = stretch_factor / 3600.0 / 180.0 * np.pi * focal_length

    new_items = [("HDUCLASS", "OGIP"), ("HDUDOC", "CAL/GEN/92-027"), ("HDUVERS", "1.0.0"),
                 ("HDUCLAS1", "IMAGE"), ("HDUCLAS2", "PSF"), ("HDUCLAS3", "PREDICTED"), ("HDUCLAS4", "NET"),
                 ("TELESCOP", "XMM-Newton"), ("INSTRUME", inst_name_map[instrument_name]), ("FILTER", "NONE"),
                 ("BACKGRND", 0.0), ('CTYPE1', 'DETX'), ('CTYPE2', 'DETY'), ('CUNIT1', 'm '),
                 ('CUNIT2', 'm'), ('CRPIX1', 61.0), ('CRPIX2', 61.0), ('CRVAL1', 0.0), ('CRVAL2', 0.0),
                 ('CDELT1', cdelt), ('CDELT2', cdelt), ('ENERG_LO', 1.0), ('ENERG_HI', 1.0), ('CHANMAX', -99),
                 ('CHANTYPE', 'PI'), ('SUMRCTS', 1.0), ('CCLS0001', 'BCF'), ('CDTP0001', 'DATA'),
                 ('CCNM0001', '2D_PSF'), ('CVSD0001', '2000-01-01'), ('CVST0001', '00:00:00'),
                 ('CDES0001', 'ELLBETA model PSF')]

    logger.info(f"Computing XMM PSF (ELLBETA) for instrument {inst_name_map[instrument_name]} in 1.0\" pixels.")
    logger.info(f"For a radius of 1', images will be 120 pixels wide, with a scale of {cdelt} m/pixel.")

    with TemporaryDirectory() as temp, open(os.devnull, "w") as dn, redirect_stdout(dn):
        run_dir = Path(temp)

        for energy in energy_list:
            energy_formatted = f"{energy:3f}"

            for theta in theta_list:
                for phi in range(0, 360, 4):
                    output = run_dir / f"PsfImg_{energy}eV_{theta}arcmin_{phi}deg.fits"
                    phi_sas = 90.0 - phi

                    if phi_sas < 0:
                        phi_sas = phi_sas + 360.0

                    phi_sas = phi_sas / 180.0 * np.pi

                    inargs = [f"instrument={inst_name_map[instrument_name]}", f"coordtype=TEL", f"x={theta}",
                              f"y={phi_sas}", f"output={output.resolve()}", f"level=ELLBETA", f"xsize=120",
                              f"ysize=120",
                              f"energy={energy}", "-V 0", "-w 0"]

                    psfgen = sas("psfgen", inargs)
                    psfgen.run()

                    with fits.open(output.resolve(), "update") as file:
                        data = file[0].data
                        header = file[0].header

                        data_sum = np.sum(data)
                        file[0].data = data / data_sum

                        header.remove("COMMENT", ignore_missing=True, remove_all=True)

                        header.extend(new_items)

                        header.insert("HDUCLASS", ("EXTNAME", f"${energy}eVthet${theta}arcsecphi${phi}deg"))
                        header.insert("BACKGRND", ("ENERGY", f"{energy}", "Energy in eV"), after=True)
                        header.insert("ENERGY", ("THETA", f"{theta / 60.0:.3f}", "Off-axis angle in arcmin"),
                                      after=True)
                        header.insert("THETA", ("PHI", f"{phi}", "Azimuth in degree"), after=True)
                        header.insert("CUNIT2", ("CRPIX1", 61.), after=True)
                        header.insert("CRPIX1", ("CRPIX2", 61.), after=True)
                        header.insert("CCNM0001", ("CBD10001", f"ENERGY( {energy_formatted})kEV"), after=True)
                        header.insert("CBD10001", ("CBD20001", f"THETA( {theta}.000000)arcsec"), after=True)
                        header.insert("CBD20001", ("CBD30001", f"PHI( {phi}.000000)deg"), after=True)

                    if out.exists():
                        with hsp.utils.local_pfiles_context(temp):
                            hsp.ftappend(infile=f"{output.resolve()}", outfile=f"{out.resolve()}", history="no",
                                         chatter=0)
                        output.unlink()
                    else:
                        output.rename(out.resolve())


def get_psf_name(
        instrument_name: Literal["epn", "emos1", "emos2"],
        res_mult: int,
) -> str:
    return f"{instrument_name}_psf_{1.0 / res_mult}x_e_0.5_2.0_kev.fits"


def add_ccdnr_and_xy(hdu, xmm, ccdnr):
    # This function adds a ccdnr to an hdu
    # Based off https://docs.astropy.org/en/stable/io/fits/usage/table.html
    data = hdu['EVENTS'].data

    # Create an numpy array of the same length with the ccdnr
    event_len = data.shape[0]
    ccd_data = np.full(event_len, ccdnr, dtype=np.uint8)

    # Convert the RAWX and RAWY data to DETX and DETY using the previously created xmm
    rawx_data = hdu['EVENTS'].data['RAWX']
    rawy_data = hdu['EVENTS'].data['RAWY']

    detx_data = np.zeros(rawx_data.shape[0], dtype=np.uint16)
    dety_data = np.zeros(rawy_data.shape[0], dtype=np.uint16)

    for i in range(rawx_data.shape[0]):
        rawx = rawx_data[i]
        rawy = rawy_data[i]

        detx, dety = xmm.ccds[ccdnr - 1].convert_raw_to_relative(rawx=rawx,
                                                                 rawy=rawy)  # remember, ccdnrs are from 1 to 12 not indexes
        detx_data[i] = detx
        dety_data[i] = dety

    # Create a fits record of the new ccnr
    ccd_col = fits.Column(name='CCDNR', format='B', array=ccd_data)

    # Create a fits record of the new DETX and DETY
    ccd_detx = fits.Column(name='DETX', format='I',
                           array=detx_data)  # name = 'DETX'; format = 'I'; unit = '0.05 arcsec'; coord_ref_point = 0; coord_inc = 1.38888888888889e-05
    ccd_dety = fits.Column(name='DETY', format='I', array=dety_data)

    coldefs = fits.ColDefs([ccd_col, ccd_detx, ccd_dety])
    hdu_ccdnr = fits.BinTableHDU.from_columns(coldefs)

    # Merge the two fits records to create the new record with the ccnrs
    new_columns = data.columns + hdu_ccdnr.columns

    header = hdu['EVENTS'].header
    # See https://docs.astropy.org/en/stable/io/fits/api/headers.html#astropy.io.fits.Header
    header['TFIELDS'] = 17
    header['TTYPE17'] = 'CCDNR'
    header.comments['TTYPE17'] = 'CCD number of the telescope'
    header['TFORM17'] = 'B'
    header.comments['TFORM17'] = 'data format of field: BYTE'

    header['TFIELDS'] = 18
    header['TTYPE18'] = 'DETX'
    header.comments['TTYPE18'] = 'Linearised Camera X-Coordinate'
    header['TFORM18'] = 'I'
    header.comments['TFORM18'] = 'data format of field: 2-byte INTEGER'

    header['TFIELDS'] = 19
    header['TTYPE19'] = 'DETY'
    header.comments['TTYPE19'] = 'Linearised Camera Y-Coordinate'
    header['TFORM19'] = 'I'
    header.comments['TFORM19'] = 'data format of field: 2-byte INTEGER'

    header['NAXIS1'] = header['NAXIS1'] + 1 + 4  # Add the extra byte of the ccdnr and two for the detx and dety
    # header['NAXIS1'] = xmm.total_width
    # header['NAXIS2'] = xmm.total_height

    new_hdu = fits.BinTableHDU.from_columns(new_columns, header=header)

    # Replace the old events with the new data
    hdu['EVENTS'] = new_hdu

    return hdu


if __name__ == '__main__':
    create_psf_file(Path("/home/bojantodorkov/Projects/xmm-epicpn-simulator/test.fits"), "epn", 2, [500, 1000, 2000],
                    [0, 210, 420, 600, 720, 900, 1200])
