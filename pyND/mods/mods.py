"""Routines for accessing, analyzing MODS spectroscopy.

wave,flux,err = read_mods1d(input_file,header=False)
                -- Read individual MODS 1D spectra.

wave, flux, err = join_mods1d(blue_file, red_file, object_number=None, header=False):
                -- Join blue+red spectra.
"""
from __future__ import print_function, absolute_import, division, unicode_literals
from ..plotting import plotzero, plotaxes

def read_mods1d(input_file,header=False):
    """Read data from mods 1D output format

    :param input_file: Input filename.
    :param header(=False): return the header with True.
    :return: wave,flux,err [optional: header]
    """

    import numpy as np
    from astropy.io import fits

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

def join_mods1d(blue_file, red_file, object_number=None, header=False):
    """Combine red and blue MODS spectra into a single spectrum from 3,200 to 10,000 Ang."""

    import numpy as np
    from astropy.io import fits

    # Read in the blue, red spectra separately.
    blue_wave,blue_flux,blue_err, blue_hdr = read_mods1d(blue_file,header=True)
    red_wave,red_flux,red_err, red_hdr = read_mods1d(red_file,header=True)

    if object_number != None:
        i=object_number
        blue_wave, blue_flux, blue_err = blue_wave[i],blue_flux[i],blue_err[i]
        red_wave, red_flux, red_err = red_wave[i], red_flux[i], red_err[i]

    # How many objects:
    num_obj = (np.shape(red_wave))[0]

    # TODO: Consider allowing wavelength shifts:

    # join_mods1d(blue_file, red_file, object_number=None, header=False, blue_shift=None, red_shift=None):

    # Check to see if any shifts are input:
    # if blue_shift != None:
    #     num_shifts = np.size(blue_shift)
    #     if num_shifts == 1:
    #         blue_wave = np.array(blue_wave)+blue_shift
    #     else:
    #         blue_wave = np.array(blue_wave)
    #         for j in np.arange(num_shifts):
    #             blue_wave[j,:] = blue_wave[j,:]+blue_shift[j]
    #
    # if red_shift != None:
    #     num_shifts = np.size(red_shift)
    #     if num_shifts == 1:
    #         red_wave = np.array(red_wave) + red_shift
    #     else:
    #         red_wave = np.array(red_wave)
    #         for j in np.arange(num_shifts):
    #             red_wave[j, :] = red_wave[j, :] + red_shift[j]

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
        iRed = np.where((red_wave[j] >= 5770.) & (red_wave[j] <= np.max(obj_wave[j])))
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
