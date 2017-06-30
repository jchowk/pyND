import numpy as np
from astropy.io import fits
import os

def solarabundance(input, error=False, best=True, photo=False, meteor=False):
    """
    eps=solarabundance(input,error=False,photo=False,meteor=False)
    
    :param input: Input elements (by symbol or Z) 
    :param error: If True, return the error in abundance [default: False]    
    :param photo: If True, force photospheric abundances [default: False]
    :param meteor: If True, force meteoritic abundances [default: False] 
    :return: abundance or abundance, error
    """

    def _ret():
        if error:
            return bestabundance, besterr
        else:
            return bestabundance

    # Where to find the abundance data:
    # Could have used: data_dir = analysis.__path__[0]+'/data/'
    data_dir = os.path.join(os.path.dirname(__file__), 'data/')
    abundancefile = data_dir+'asplund2009_abundances.fits'

    # Read in the data:
    a = fits.getdata(abundancefile)
    
    ##Define output variables:
    bestabundance = np.zeros(np.size(input))
    besterr = np.zeros(np.size(input))

    index = np.zeros(np.size(input),dtype=int)

    # Inputs can be either string of element name or nuclear charge:
    if type(input[0]) == str:
        for j in np.arange(np.size(input)):
            temp = np.where(a['ELEMENT'] == input[j].lower())
            if np.size(temp) == 0:
                index[j]=999
            else:
                index[j] = temp[0][0]
    else:
        for j in np.arange(np.size(input)):
            temp = np.where(a['Z'] == input[j])
            if size(temp) == 0:
                index[j] = 999
            else:
                index[j] = temp[0][0]

    for k in np.arange(np.size(index)):
        # If the input is crap, return crap.
        if index[k] == 999:
            bestabundance[k] = -30.
            besterr[k] = -30.
        else:
            # Did the user force photospheric results?
            if (photo is True):
                bestabundance[k] = a['PHOTO'][index[k]]
                besterr[k] = a['PHOTO_ERR'][index[k]]
            # Did the user force photospheric results?
            elif (meteor is True):
                bestabundance[k] = a['METEOR'][index[k]]
                besterr[k] = a['METEOR_ERR'][index[k]]
            # If neither of the above, adopt the best abundances
            else:
                bestabundance[k] = a['BEST'][index[k]]
                besterr[k] = a['ERR'][index[k]]


    bestabundance = np.around(bestabundance,3)
    besterr= np.around(besterr, 3)

    return _ret()

