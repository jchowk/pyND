import numpy as np
from astropy.io import fits

def read_mods(input_file):

    hdr = fits.getheader(input_file)
    flux = (fits.getdata(input_file,0))[1]
    err = (fits.getdata(input_file,1))[1]

    # Create the wavelength array [VACUUM]
    wave = np.arange(np.size(flux))*hdr.get('CDELT1')+hdr.get('CRVAL1')

    def _ret():  return (wave, flux, err)

    return _ret()

