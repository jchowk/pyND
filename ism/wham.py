from __future__ import print_function, absolute_import, division, unicode_literals
from pyND.abs import plotzero, plotaxes

# TODO: Allow an averaging of nearby spectra.
# TODO: Create environment variable for data location.
def get_spectrum(longitude, latitude,silent=False):
    """Extract a single WHAM spectrum from the WHAM-SS.

    :returns: (lsr_velocity[km/s], intensity[R/(km/s)], error[R/(km/s)])
    """

    import numpy as np

    from astropy.coordinates import SkyCoord, match_coordinates_sky
    from astropy.io import fits
    import astropy.units as u

    import pyND.ism as ism


    # Set up the location of the data:
    wham_dir = ism.__path__[0]+'/data/'
    FITSfile = 'wham-ss-DR1-v161116-170316.fits.gz'
    # Load WHAM
    wham_structure = fits.open(wham_dir+FITSfile)
    wham_survey = wham_structure[1].data

    # Create our coordinate lists
    input_coords = SkyCoord(longitude,latitude,frame='galactic',unit=u.deg)
    wham_coords = SkyCoord(wham_survey['GAL-LON'],wham_survey['GAL-LAT'],
                           frame='galactic',unit=u.deg)

    # Find the closest pointing:
    idx, d2d, d3d = match_coordinates_sky(input_coords,wham_coords)

    vel_out = wham_survey['VELOCITY'][idx]
    int_out = wham_survey['DATA'][idx]
    err_out = np.sqrt(wham_survey['VARIANCE'][idx])

    coords_out = wham_coords[idx]
    totintens_out = wham_survey['INTEN'][idx]

    if silent != None:
        print("Pointing centered at (l,b) = ({0:0.1f}, {1:0.1f}) -- {2:0.1f} from target".format(
            coords_out.l.value,coords_out.b.value,d2d[0]))
        print("Total intensity = {0:0.1f} mR".format(totintens_out*1.e3))

    if wham_survey['STAR'][idx] != False:
            print('Possible stellar contamination.')

    def _ret():
        return (vel_out,int_out,err_out)

    return _ret()


def get_intensity(longitude, latitude, silent=False):
    """Extract the integrated intensity for a longitude, latitude.

    :returns: Intensity[R]
    """

    import numpy as np
    from astropy.io import fits
    import pyND.ism as ism

    # Set up the location of the data:
    wham_dir = ism.__path__[0] + '/data/'
    FITSfile = 'wham-ss-DR1-v161116-170316-int-grid.fits.gz'
    # Load WHAM
    wham_structure = fits.open(wham_dir + FITSfile)
    header = wham_structure[0].header
    wham_intensity = wham_structure[0].data

    # Calculate the pixel coordinates for our position.
    #   NOTE: the CD1_2, CD2_1 values are known to be zero!
    xpos = np.int((longitude-header['CRVAL1'])/header['CD1_1']+header['CRPIX1'])
    ypos = np.int((latitude-header['CRVAL2'])/header['CD2_2']+header['CRPIX2'])

    intensity_out = wham_intensity[xpos,ypos]

    if silent != None:
        print("Total intensity = {0:0.1f} mR".format(intensity_out * 1.e3))

    def _ret():
        return intensity_out

    return _ret()
