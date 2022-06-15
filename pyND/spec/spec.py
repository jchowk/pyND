from ..plotting import plotzero,plotaxes
from numpy import mean

def rebin(myarr,factor,estimator=mean):
    """

    rebin(myarray, binfactor, estimator=mean)

    Rebin an input spectral element. Each element (wave, flux, variance)
    must be entered separately. A 1D spectrum or 2D collection of spectra
    [n_spectra, n_wavelengths] may be passed.

    KEYWORDS:
        estimator - default to meanp.  You can downsample by summing or
            something else if you want a different estimator
            (e.g., downsampling error: you want to sum & divide by sqrt(n))
    """
    import numpy as np

    shape = myarr.shape
    #Test for 1D shape:
    if len(shape)==1:
        xs=shape[0]
        crarr = myarr[:xs-(xs % int(factor))]
        dsarr = np.zeros(int(xs/factor))
        dsarr = estimator(
            np.concatenate([[crarr[j::factor]
                                for j in range(factor)]]),
            axis=0)
    else:
        ys,xs=shape

        crarr = myarr[:,:xs-(xs % int(factor))]
        dsarr = np.zeros([ys,int(xs/factor)])
        for i in np.arange(ys):
            dsarr[i,:] = estimator(
                np.concatenate([[crarr[i,j::factor]
                                    for j in range(factor)]]), axis=0)

    return dsarr

def congrid(a, newdims, method='linear', centre=True, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.

    Adopted from sci-py cookbook: http://scipy-cookbook.readthedocs.io/items/Rebinning.html .
    '''
    import numpy as np

    import scipy.interpolate
    import scipy.ndimage

    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print("[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions.")
        return None
    newdims = np.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = np.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = np.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = np.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = np.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = np.mgrid[nslices]

        newcoords_dims = range(np.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print("Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported.")
        return None




def resample(new_xxx, xxx, flux, err=None,
            do_sig=False,grow_bad_sig=True,
            fill_value=None,
            **kwargs):
    """ ADAPTED FROM LINETOOLS REBIN CODE.
    [https://github.com/linetools/linetools]

    Resample a single spectrum to a new velocity or wavelength array.

    Uses simple linear interpolation.  The default (and only)
    option conserves counts (and flambda).

    WARNING: Do not trust either edge pixel of the new array.

    Also be aware that neighboring pixels are likely to be
    correlated in a manner that is not described by the error
    array.

    Parameters
    ----------
    new_xxx : array
      New velocity / wavelength array

    xxx : array
      Input velocity / wavelength array


    fill_value : float, optional
      Fill value at the edges
      Default = None [filling with NAN]

    grow_bad_sig : bool, optional
      Allow sig<=0. values and grow them

    """
    from scipy.interpolate import interp1d
    import numpy as np
    import warnings

    # Deal with nan
    badf = np.any([np.isnan(flux), np.isinf(flux)], axis=0)
    if np.sum(badf) > 0:
        warnings.warn("Ignoring pixels with NAN or INF in flux")
    # Select the good data
    gdf = ~badf
    flux = flux[gdf]

    # Check for bad pixels (not prepared for these)
    if err is not None:
        sig = err.copy()
        bad_sig = sig[gdf] <= 0.
        if np.sum(bad_sig) > 0:
            if not grow_bad_sig:
                raise IOError("Data contains rejected pixels (sig=0). Use grow_bad_sig to proceed and grow them.")
        bads = np.any([np.isnan(sig[gdf]), np.isinf(sig[gdf]**2)], axis=0)  # Latter is for way too large values
        bad_sig[bads] = True

    # Endpoints of original pixels
    npix = len(xxx)

    # Average velocity positions
    vlh = (xxx + np.roll(xxx, -1)) / 2.
    vlh[npix - 1] = xxx[npix - 1] + \
                    (xxx[npix - 1] - xxx[npix - 2]) / 2.
    # Delta xxx
    dvl = vlh - np.roll(vlh, 1)
    dvl[0] = 2 * (vlh[0] - xxx[0])
    med_dvl = np.median(dvl)

    # Select "good" data points – those not NAN or INF
    vlh = vlh[gdf]
    dvl = dvl[gdf]

    # We "bin" the variance rather than error
    if do_sig:
        if err is None:
            raise IOError("err must be set to rebin uncertainties")
        var = sig[gdf]**2
        var[bad_sig] = 0.
    else:
        var = np.ones_like(flux)



    # To conserve flux, use the cumulative sum as a function of xxx as
    # the basis for the interpolation.

    # Cumulative sum of the flux, variance
    cumsum = np.cumsum(flux * dvl)
    cumvar = np.cumsum(var * dvl, dtype=np.float64)


    # Interpolate the cumulative sum FILL_VALUE should probably be 0.
    fcum = interp1d(vlh, cumsum, fill_value=0., bounds_error=False)
    fvar = interp1d(vlh, cumvar, fill_value=0., bounds_error=False)

    # Create a reference interpolation to fill/flag pixels outside the range
    # of the original data.
    fcum_ref = interp1d(vlh, cumsum, fill_value=fill_value, bounds_error=False)

    # Endpoints of new pixels
    nnew = len(new_xxx)
    nvlh = (new_xxx + np.roll(new_xxx, -1)) / 2.
    nvlh[nnew - 1] = new_xxx[nnew - 1] + \
                     (new_xxx[nnew - 1] - new_xxx[nnew - 2]) / 2.
    # Pad starting point
    bvl = np.zeros(nnew + 1)
    bvl[0] = new_xxx[0] - (new_xxx[1] - new_xxx[0]) / 2.
    bvl[1:] = nvlh

    # Evaluate
    newcum = fcum(bvl)
    newvar = fvar(bvl)

    # Rebinned flux
    new_fx = (np.roll(newcum, -1) - newcum)[:-1]
    new_var = (np.roll(newvar, -1) - newvar)[:-1]

    # Normalize (preserve counts and flambda)
    new_dvl = bvl - np.roll(bvl, 1)
    new_fx = new_fx / new_dvl[1:]
    # Preserve S/N (crudely)
    med_newdvl = np.median(new_dvl)
    new_var = new_var / (med_newdvl/med_dvl) / new_dvl[1:]


    # Deal with the regions beyond the original data
    if fill_value == None:
        bd_vel = np.isnan(fcum_ref(bvl)[:-1])
        new_fx[bd_vel] = np.nan
        new_var[bd_vel] = np.nan
    else:
        bd_vel = (fcum_ref(bvl) == fill_value)[:-1]
        new_fx[bd_vel] = fill_value
        new_var[bd_vel] = fill_value

    # Return new spectrum
    if do_sig:
        # Create new_sig
        new_sig = np.zeros_like(new_var)
        gd = new_var > 0.
        new_sig[gd] = np.sqrt(new_var[gd])
        # Deal with bad pixels (grow_bad_sig should be True)
        bad = np.where(var <= 0.)[0]
        # Find nearby wavelengths in rebinned wavelength
        nearidxs = np.searchsorted(new_xxx, xxx[bad])
        # Pad arrays to enable vector operations
        pvl = np.concatenate([xxx,[xxx[-1]+dvl[-1]]])
        pndvl = np.concatenate([new_dvl,[new_dvl[-1]]])
        pnvl = np.concatenate([new_xxx,[new_xxx[-1]+new_dvl[-1]]])

        # Find distances between original bad wavelengths and nearby new ones
        ldiff = np.abs(new_xxx[nearidxs-1]-pvl[bad]) - \
                (pndvl[1:][nearidxs]+dvl[bad])/2
        rdiff = np.abs(pvl[bad]-pnwv[nearidxs]) - \
                (pndvl[1:][nearidxs] + dvl[bad]) / 2
        # Set errors to 0; we have to mind the padding above
        new_sig[nearidxs[(ldiff<0)&(nearidxs<len(new_xxx))]] = 0
        new_sig[nearidxs[(rdiff<0)&(nearidxs<len(new_xxx))]] = 0
        ### Old (very slow) way looping through bad pix
        #for ibad in bad:
        #    bad_new = np.where(np.abs(new_xxx-xxx[ibad]) <
        #                       (new_dvl[1:]+dvl[ibad])/2)[0]
        #    new_sig[bad_new] = 0.

        # Zero out edge pixels -- not to be trusted
        igd = np.where(gd)[0]
        if len(igd) == 0:  # Should not get here!
            raise ValueError("Not a single good pixel?!  Something went wrong...")
        new_sig[igd[0]] = 0.
        new_sig[igd[-1]] = 0.
    else:
        new_sig = None


    if do_sig:
        return new_fx, new_sig
    else:
        return new_fx


def vac_to_air(lam_vac):
    """
    Convert vacuum to air wavelengths using
    equation (1) of Ciddor 1996, Applied Optics 35, 1566
        http://dx.doi.org/10.1364/AO.35.001566

    :param lam_vac - Wavelength in Angstroms
    :return: lam_air - Wavelength in Angstroms

    ** Borrowed from M. Cappellari's ppxf_util software.

    """
    sigma2 = (1e4/lam_vac)**2
    fact = 1 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)

    return lam_vac/fact

###############################################################################

def air_to_vac(lam_air):
    """
    Convert air to vacuum wavelengths using
    equation (1) of Ciddor 1996, Applied Optics 35, 1566
        http://dx.doi.org/10.1364/AO.35.001566
    :param lam_air - Wavelength in Angstroms
    :return: lam_vac - Wavelength in Angstroms

    ** Borrowed from M. Cappellari's ppxf_util software.

    """
    sigma2 = (1e4/lam_air)**2
    fact = 1 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)

    return lam_air*fact
