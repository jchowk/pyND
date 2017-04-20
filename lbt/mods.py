import numpy as np
from astropy.io import fits

def read_mods1d(input_file,header=False):
    """Read data from mods 1D output format
    
    :param input_file: 
    :return: 
    """

    hdr = fits.getheader(input_file)
    flux = (fits.getdata(input_file,0))
    err = (fits.getdata(input_file,1))

    # How many spectra, data points?
    num_arrays = (np.shape(flux))[0]
    num_points = (np.shape(flux))[1]

    # Exclude the calibration spec [index 0]
    flux = flux[1:num_arrays]*1.e-17
    err = err[1:num_arrays] * 1.e-17

    # Create the wavelength array [VACUUM]
    wave = []
    wave_base = np.arange(num_points)*hdr.get('CDELT1')+hdr.get('CRVAL1')
    for j in np.arange(1,num_arrays):
        wave.append(wave_base)

    if header == False:
        def _ret():  return (wave, flux, err)
    else:
        def _ret():  return (wave, flux, err, hdr)

    return _ret()

def join_mods1d(blue_file,red_file,object_numbers=None,header=False):

    # Object extraction:
    blue_wave,blue_flux,blue_err, blue_hdr = read_mods1d(blue_file,header=True)
    red_wave,red_flux,red_err, red_hdr = read_mods1d(red_file,header=True)

    if object_numbers != None:
        i=object_number
        blue_wave, blue_flux, blue_err = blue_wave[i],blue_flux[i],blue_err[i]
        red_wave, red_flux, red_err = red_wave[i], red_flux[i], red_err[i]

    # How many objects:
    num_obj = (np.shape(red_wave))[0]

    # Create output arrays. MODS is on an evenly-spaced linear grid.
    obj_wave = []
    obj_wave_init = np.arange(3200.,10000.5,0.5)
    for j in np.arange(num_obj):
        obj_wave.append(obj_wave_init)

    # Weighting vectors
    obj_inv_var = np.zeros_like(obj_wave)
    obj_weighted_flux = np.zeros_like(obj_wave)

    # Output vectors
    obj_flux = []
    obj_err = []

    # Loop through each object
    for j in np.arange(num_obj):
        # Put in the blue data:
        gBlue = np.where(obj_wave[j] <= np.max(blue_wave[j]))
        obj_inv_var[j][gBlue] += 1./blue_err[j]**2
        obj_weighted_flux[j][gBlue] += blue_flux[j]/blue_err[j]**2

        # Put in the red data:
        iRed = np.where(red_wave[j] >= 5770.)
        gRed = np.where(obj_wave[j] >= np.min(red_wave[j][iRed]))
        obj_inv_var[j][gRed] += 1./red_err[j][iRed]**2
        obj_weighted_flux[j][gRed] += red_flux[j][iRed]/red_err[j][iRed]**2

        obj_flux.append(obj_weighted_flux[j] / obj_inv_var[j])
        obj_err.append(np.sqrt(1./obj_inv_var[j]))

    # Keep the arrays simple if there's only one object:
    if num_obj == 1:
        obj_wave = obj_wave[0]
        obj_flux = obj_flux[0]
        obj_err = obj_err[0]

    if header == False:
        def _ret():
            return (obj_wave, obj_flux, obj_err)
    else:
        def _ret():
            return (obj_wave, obj_flux, obj_err, blue_hdr)

    return _ret()