import numpy as np
from astropy.io import fits


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
