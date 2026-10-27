
def dend_to_fits(d, filename):
    """Export the dendrogram 'd' to the FITS file 'filename'"""
    from astropy.io import fits
    import numpy as np

    try:
        primary_hdu = fits.PrimaryHDU(header=d.wcs.to_header())
    except AttributeError:
        primary_hdu = fits.PrimaryHDU()

    primary_hdu.header["MIN_NPIX"] = (d.params['min_npix'],
                                      "Minimum number of pixels in a leaf.")
    primary_hdu.header["MIN_DELT"] = (d.params['min_delta'],
                                       "Minimum branch height.")
    primary_hdu.header["MIN_VAL"] = (d.params['min_value'],
                                       "Minimum intensity value.")

    hdus = [primary_hdu,
            fits.ImageHDU(d.data, name='data'),
            fits.ImageHDU(d.index_map, name='index_map'),
            fits.ImageHDU(np.array([ord(x) for x in d.to_newick()]), name='newick')]

    hdulist = fits.HDUList(hdus)

    hdulist.writeto(filename, overwrite=True)