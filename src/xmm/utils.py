import os
import shutil
from pathlib import Path
from typing import Literal

import numpy as np
from astropy.io import fits

from src.config import EnergySettings
from src.xmm.ccf import get_xrt_xareaef

available_instruments = ["epn", "emos1", "emos2"]
instrument_to_sixte_dir = {"epn": "epicpn", "emos1": "epicmos", "emos2": "epicmos"}


def get_fov(instrument_name: str) -> float:
    if instrument_name == "epn":
        from src.xmm.epn import get_fov

        return get_fov()

    if instrument_name == "emos1" or instrument_name == "emos2":
        from src.xmm.emos import get_fov

        return get_fov(int(instrument_name[-1]))

    raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")


def get_cdelt(instrument_name: str, res_mult: int) -> tuple[float, float]:
    if instrument_name == "epn":
        from src.xmm.epn import get_cdelt

        return get_cdelt(res_mult)

    if instrument_name == "emos1" or instrument_name == "emos2":
        from src.xmm.emos import get_cdelt

        return get_cdelt(res_mult)

    raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")


def get_pixel_size(instrument_name: str, res_mult: int) -> float:
    if instrument_name == "epn":
        from src.xmm.epn import get_pixel_size

        return get_pixel_size(res_mult)

    if instrument_name == "emos1" or instrument_name == "emos2":
        from src.xmm.emos import get_pixel_size

        return get_pixel_size(int(instrument_name[-1]), res_mult)

    raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")


def get_naxis12(instrument_name: str, res_mult: int) -> tuple[int, int]:
    if instrument_name == "epn":
        from src.xmm.epn import get_naxis12

        return get_naxis12(res_mult=res_mult)

    if instrument_name == "emos1" or instrument_name == "emos2":
        from src.xmm.emos import get_naxis12

        return get_naxis12(emos_num=int(instrument_name[-1]), res_mult=res_mult)

    raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")


def get_crpix12(instrument_name: str, res_mult: int):
    naxis1, naxis2 = get_naxis12(instrument_name=instrument_name, res_mult=res_mult)

    if instrument_name == "epn":
        from src.xmm.epn import get_shift_xy

        shift_y, shift_x = get_shift_xy(res_mult=res_mult)
        crpix1 = round(((naxis1 + 1) / 2.0) - shift_y, 6)
        crpix2 = round(((naxis2 + 1) / 2.0) + shift_x, 6)
    else:
        crpix1 = round(((naxis1 + 1) / 2.0), 6)
        crpix2 = round(((naxis2 + 1) / 2.0), 6)

    return crpix1, crpix2


def get_focal_length(instrument_name: str) -> float:
    if instrument_name == "epn":
        from src.xmm.epn import get_focal_length

        return get_focal_length()

    if instrument_name == "emos1" or instrument_name == "emos2":
        from src.xmm.emos import get_focal_length

        return get_focal_length(emos_num=int(instrument_name[-1]))

    raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")


def get_instrument_files(instrument_name: str) -> Path:
    """
    Returns:
        Path: Path to the SIXTE instruments directory.
    Raises:
        NotADirectoryError: If the instrument files have not been downloaded from
            https://www.sternwarte.uni-erlangen.de/sixte/instruments/
    """
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    p = Path(os.environ["SIXTE"]) / "share" / "sixte" / "instruments" / "xmm" / instrument_to_sixte_dir[instrument_name]

    if not p.exists():
        raise NotADirectoryError(
            "It looks like you haven't downloaded the instrument files provided by SIXTE "
            "(see https://www.sternwarte.uni-erlangen.de/sixte/instruments/)! Please download "
            "and extract them as given in their instructions."
        )

    return p


def create_vinget_file(instrument_name: str, xml_dir: Path) -> None:
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    out_file = get_vignet_file(xml_dir=xml_dir.resolve(), instrument_name=instrument_name)
    xrt_xareaef = get_xrt_xareaef(instrument_name=instrument_name)

    with fits.open(name=xrt_xareaef, mode="readonly") as file:
        vignetting = file["VIGNETTING"]
        energy = vignetting.data["ENERGY"]
        theta_stepsize = vignetting.header["D_THETA"]
        theta_bins = vignetting.data["VIGNETTING_FACTOR"].shape[1]
        ccf_vignet_data = vignetting.data["VIGNETTING_FACTOR"]

    # Transform them to the sixte format
    energy_lo = np.array(energy[:-1] / 1000.0).reshape(1, len(energy[:-1]))  # E to keV
    energy_hi = np.array(energy[1:] / 1000.0).reshape(1, len(energy[1:]))  # # E to keV

    # Create the theta step_size bins
    theta = [0.0]
    for _ in range(theta_bins - 1):
        theta.append(theta[-1] + theta_stepsize)

    theta = np.array(theta).reshape(1, len(theta))

    # We do not have a phi, thus an array of zero
    phi = np.zeros((1, 1))

    # Create the vignet array
    vignet = np.zeros((1, theta.shape[1] - 1, energy_lo.shape[1]))

    for i in range(ccf_vignet_data.shape[0] - 1):
        # i indexes the theta variable
        for j in range(ccf_vignet_data.shape[1] - 1):
            # j intexes the energy bin variable
            # Convert to keV
            vignet[0][j][i] = ccf_vignet_data[i][j]

    # The format is based on: https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_021/cal_gen_92_021.html

    energy_lo_col = fits.Column(
        "ENERG_LO",
        format=f"{energy_lo.shape[1]}E",
        unit="keV",
        dim=f"{energy_lo.shape[1]}",
        array=energy_lo,
    )
    energy_hi_col = fits.Column(
        "ENERG_HI",
        format=f"{energy_hi.shape[1]}E",
        unit="keV",
        dim=f"{energy_hi.shape[1]}",
        array=energy_hi,
    )
    theta_col = fits.Column(
        "THETA",
        format=f"{theta.shape[1]}E",
        unit="degree",
        dim=f"{theta.shape[1]}",
        array=theta,
    )
    phi_col = fits.Column(
        "PHI",
        format=f"{phi.shape[1]}E",
        unit="degree",
        dim=f"{phi.shape[1]}",
        array=phi,
    )
    vignet_col = fits.Column(
        "VIGNET",
        format=f"{vignet.shape[0] * vignet.shape[1] * vignet.shape[2]}E",
        dim=f"{vignet.shape[::-1]}",
        array=vignet,
    )

    coldefs = fits.ColDefs([energy_lo_col, energy_hi_col, theta_col, phi_col, vignet_col])
    vignet_hdr = fits.Header()
    # Include the required headers
    vignet_hdr["HDUCLASS"] = "OGIP"
    vignet_hdr["HDUCLAS1"] = "RESPONSE"
    vignet_hdr["HDUVERS1"] = "1.0.0"
    vignet_hdr["HDUCLAS2"] = "VIGNET"
    vignet_hdr["HDUVERS2"] = "1.1.0"
    vignet_hdr["VERSION"] = "20171016"
    vignet_hdr["MISSION"] = "XMM"
    vignet_hdr["TELESCOP"] = "XMM"
    vignet_hdr["DETNAM"] = "XMM"
    vignet_hdr["INSTRUME"] = "EPN"

    vignet_hdr["TUNIT1"] = "keV"
    vignet_hdr["TUNIT2"] = "keV"
    vignet_hdr["TUNIT3"] = "degree"
    vignet_hdr["TUNIT4"] = "degree"

    # vignet_hdr['INSTRUME'] = 'EPIC-PN'
    vignet_hdr.add_history(
        "Produced by Sam Sweere (ESAC Trainee) according to data from the XMM-PN calibration file XRT3_XAREAEF_0012.CCF"
    )
    vignet_hdu = fits.BinTableHDU.from_columns(coldefs, header=vignet_hdr)
    vignet_hdu.name = "VIGNET"

    # Create the final hdul
    primary_hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([primary_hdu, vignet_hdu])

    # Save the new fits file
    hdul.writeto(out_file, overwrite=True)


def create_psf_file(instrument_name: str, xml_dir: Path, res_mult: int) -> None:
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    out = get_psf_file(xml_dir=xml_dir.resolve(), instrument_name=instrument_name, res_mult=res_mult)
    inst_name_map = {"epn": "pn", "emos1": "mos1", "emos2": "mos2"}

    instrument_files = get_instrument_files(instrument_name=instrument_name)
    psf_file = list(instrument_files.glob(f"*_{inst_name_map[instrument_name]}_psf.fits"))[0]

    shutil.copy(src=psf_file, dst=out)

    if res_mult != 1:
        with fits.open(out, mode="update") as hdu_list:
            for primary_hdu in hdu_list:
                primary_hdu.header["CDELT1"] = primary_hdu.header["CDELT1"] / res_mult
                primary_hdu.header["CDELT2"] = primary_hdu.header["CDELT2"] / res_mult


def create_mask(
    instrument_name: str,
    observation_id: str,
    mask_level: str | None,
    energies: EnergySettings,
    res_mults: list[int] = None,
) -> dict[str, dict[int, Path]] | None:
    if mask_level is None:
        return None

    if res_mults is None:
        res_mults = [1]

    out_dir = Path("res").resolve()
    if instrument_name == "epn":
        from src.xmm.epn import create_mask

        return {
            instrument_name: create_mask(energies.emin, energies.emax, observation_id, out_dir, mask_level, res_mults)
        }

    if instrument_name == "emos1" or instrument_name == "emos2":
        from src.xmm.emos import create_mask

        return {
            instrument_name: create_mask(
                energies.emin, energies.emax, instrument_name[-1], observation_id, out_dir, mask_level, res_mults
            )
        }

    raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")


def create_xml_files(
    instrument_name: str,
    xml_dir: Path,
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
    wait_time: float,
) -> list[Path]:
    """
    Returns:
        List[Path]: A list containing paths to the corresponding CCDs. If sim_separate_ccds == True, then the list
            contains a single Path to /path/to/instrument/combined.xml
    """
    instrument_files = get_instrument_files(instrument_name=instrument_name)
    instrument_dir = xml_dir / instrument_name
    out_dir = instrument_dir / xmm_filter / f"{res_mult}x"
    instrument_dir.mkdir(exist_ok=True, parents=True)
    out_dir.mkdir(exist_ok=True, parents=True)

    vignet_file = get_vignet_file(xml_dir=xml_dir, instrument_name=instrument_name)
    vignet_name = vignet_file.name
    psf_file = get_psf_file(xml_dir=xml_dir, instrument_name=instrument_name, res_mult=res_mult)
    psf_name = psf_file.name

    prefix = f"{instrument_name[1:]}"
    files = {
        instrument_files / f"{prefix}-{xmm_filter}-10.rmf": out_dir / f"{prefix}-{xmm_filter}-10.rmf",
        instrument_files / f"{prefix}-{xmm_filter}-10.arf": out_dir / f"{prefix}-{xmm_filter}-10.arf",
        psf_file: out_dir / psf_name,
        vignet_file: out_dir / vignet_name,
    }
    # Add symbolic links
    for src, dest in files.items():
        dest.unlink(missing_ok=True)
        dest.symlink_to(src)

    if instrument_name == "epn":
        from src.xmm.epn import create_xml

        return create_xml(
            out_dir=out_dir,
            res_mult=res_mult,
            xmm_filter=xmm_filter,
            sim_separate_ccds=sim_separate_ccds,
            wait_time=wait_time,
        )

    if instrument_name == "emos1" or instrument_name == "emos2":
        from src.xmm.emos import create_xml

        return create_xml(
            out_dir=out_dir,
            emos_num=int(instrument_name[-1]),
            res_mult=res_mult,
            xmm_filter=xmm_filter,
            sim_separate_ccds=sim_separate_ccds,
            wait_time=wait_time,
        )

    raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")


def get_xml_file(
    instrument_name: str,
    xml_dir: Path,
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
) -> Path:
    """
    Returns:
        List[Path]: A list containing paths to the corresponding CCDs. If sim_separate_ccds == True, then the list
            contains a single Path to /path/to/instrument/combined.xml
    """
    if instrument_name == "epn":
        from src.xmm.epn import get_xml

        return get_xml(
            xml_dir=xml_dir,
            res_mult=res_mult,
            xmm_filter=xmm_filter,
            sim_separate_ccds=sim_separate_ccds,
        )

    if instrument_name == "emos1" or instrument_name == "emos2":
        from src.xmm.emos import get_xml

        return get_xml(
            xml_dir=xml_dir,
            emos_num=int(instrument_name[-1]),
            res_mult=res_mult,
            xmm_filter=xmm_filter,
            sim_separate_ccds=sim_separate_ccds,
        )

    raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")


def get_psf_file(
    instrument_name: str,
    xml_dir: Path,
    res_mult: int,
) -> Path:
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    return xml_dir / f"{instrument_name}_psf_{1.0 / res_mult}x.fits"


def get_vignet_file(instrument_name: str, xml_dir: Path) -> Path:
    if instrument_name not in available_instruments:
        raise ValueError(f"Unknown instrument '{instrument_name}'! Available instruments: {available_instruments}.")

    return xml_dir / f"{instrument_name}_vignet.fits"


def add_ccdnr_and_xy(hdu, xmm, ccdnr):
    # This function adds a ccdnr to an hdu
    # Based off https://docs.astropy.org/en/stable/io/fits/usage/table.html
    data = hdu["EVENTS"].data

    # Create an numpy array of the same length with the ccdnr
    event_len = data.shape[0]
    ccd_data = np.full(event_len, ccdnr, dtype=np.uint8)

    # Convert the RAWX and RAWY data to DETX and DETY using the previously created xmm
    rawx_data = hdu["EVENTS"].data["RAWX"]
    rawy_data = hdu["EVENTS"].data["RAWY"]

    detx_data = np.zeros(rawx_data.shape[0], dtype=np.uint16)
    dety_data = np.zeros(rawy_data.shape[0], dtype=np.uint16)

    for i in range(rawx_data.shape[0]):
        rawx = rawx_data[i]
        rawy = rawy_data[i]

        detx, dety = xmm.ccds[ccdnr - 1].convert_raw_to_relative(
            rawx=rawx, rawy=rawy
        )  # remember, ccdnrs are from 1 to 12 not indexes
        detx_data[i] = detx
        dety_data[i] = dety

    # Create a fits record of the new ccnr
    ccd_col = fits.Column(name="CCDNR", format="B", array=ccd_data)

    # Create a fits record of the new DETX and DETY
    ccd_detx = fits.Column(
        name="DETX", format="I", array=detx_data
    )  # name = 'DETX'; format = 'I'; unit = '0.05 arcsec'; coord_ref_point = 0; coord_inc = 1.38888888888889e-05
    ccd_dety = fits.Column(name="DETY", format="I", array=dety_data)

    coldefs = fits.ColDefs([ccd_col, ccd_detx, ccd_dety])
    hdu_ccdnr = fits.BinTableHDU.from_columns(coldefs)

    # Merge the two fits records to create the new record with the ccnrs
    new_columns = data.columns + hdu_ccdnr.columns

    header = hdu["EVENTS"].header
    # See https://docs.astropy.org/en/stable/io/fits/api/headers.html#astropy.io.fits.Header
    header["TFIELDS"] = 17
    header["TTYPE17"] = "CCDNR"
    header.comments["TTYPE17"] = "CCD number of the telescope"
    header["TFORM17"] = "B"
    header.comments["TFORM17"] = "data format of field: BYTE"

    header["TFIELDS"] = 18
    header["TTYPE18"] = "DETX"
    header.comments["TTYPE18"] = "Linearised Camera X-Coordinate"
    header["TFORM18"] = "I"
    header.comments["TFORM18"] = "data format of field: 2-byte INTEGER"

    header["TFIELDS"] = 19
    header["TTYPE19"] = "DETY"
    header.comments["TTYPE19"] = "Linearised Camera Y-Coordinate"
    header["TFORM19"] = "I"
    header.comments["TFORM19"] = "data format of field: 2-byte INTEGER"

    header["NAXIS1"] = header["NAXIS1"] + 1 + 4  # Add the extra byte of the ccdnr and two for the detx and dety
    # header['NAXIS1'] = xmm.total_width
    # header['NAXIS2'] = xmm.total_height

    new_hdu = fits.BinTableHDU.from_columns(new_columns, header=header)

    # Replace the old events with the new data
    hdu["EVENTS"] = new_hdu

    return hdu
