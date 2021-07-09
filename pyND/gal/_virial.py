""" Calculate virial quantities for a dark matter halo given a
cosmology.
"""

# p2.6+ compatibility
from __future__ import division, print_function, unicode_literals
try:
    unicode
except NameError:
    unicode = basestring = str
    xrange = range

import astropy
from astropy import cosmology
from astropy.utils import isiterable
import astropy.constants as c
import astropy.units as u
import numpy as np

from collections import namedtuple
import pdb

km = 1e5
kpc = u.kpc.to('cm')

find_rvT_output = namedtuple('virial_rvT', 'r v T')

def deltavir(redshift, cosmo=None):
    """ The virial overdensity as a function redshift.
    This is an approximation parameterized by Bryan & Norman 98, ApJ,
    495, 80.  Good to 1% for omega(z) = 0.1-1, requires either omega =
    1 (flat universe) or omega_lambda = 0.
    This is given by dissipationless spherical top-hat models of
    gravitational collapse (Peebles 1990 'The Large Scale structure of
    the Universe', Eke et al. 1996 MNRAS, 282, 263).
    """
    if cosmo is None:
        cosmo = cosmology.get_current()

    if cosmo.Ok0 != 0:
        if cosmo.Ode0 == 0:
            x = cosmo.Om(redshift) - 1
            return 18*np.pi**2 + 60*x - 32*x**2
        else:
            raise ValueError("Can't compute deltavir for a non-flat cosmology "
                             "with a cosmological constant")
    else:
        x = cosmo.Om(redshift) - 1
        return 18*np.pi**2 + 82*x - 39*x**2

def _calc_rvT(M_g, rho_virial, mu=0.59):
    """Find the virial radius, circular velocity and temperature for
    a given dark matter halo mass and virial overdensity.
    Consider using find_rvT() instead - it is a higher level function
    that uses more convenient units and handles broadcasting.
    Parameters
    ----------
    M_g : float
      Halo mass in solar masses.
    rho_virial : float
      Virial overdensity.
    mu : float (default 0.59)
      Mean molecular weight in proton masses.
    Returns
    -------
    r,v,T : floats
      The virial radius (proper kpc), circular velocity (km/s) and
      temperature (K).
    """

    #Use constants from astropy
    kboltz = c.k_B.to('erg/K').value
    G = c.G.to("cm**3/(g s**2)").value
    mp = c.m_p.to('g').value

    rvir = ((3 * M_g) / (4 * np.pi * rho_virial))**(1./3)
    vcirc = np.sqrt(G * M_g / rvir)
    Tvir = mu * mp * vcirc * vcirc / (2 * kboltz)

    return rvir, vcirc, Tvir

def find_rvT(M, z, r200=False, cosmo=None, mu=0.59):
    """ Find the virial radius, circular velocity and temperature for
    a dark matter halo at a given mass and redshift.
    Parameters
    ----------
    mass : array_like, shape N
      Total halo mass (including dark matter) in solar masses.
    z : array_like, shape M
      Redshift.
    r200 : Use overdensity of 200 (optional; default = False)
    cosmo : cosmology (optional)
      The cosmology to use.
    mu : float (optional)
      The mean molecular weight. The virial temperature is
      proportional to mu.
    Returns
    -------
    vir : named tuple of ndarrays with shape (N, M)
      The virial radius (proper kpc), circular velocity (km/s) and
      temperature (K). If N or M is 1, that dimension is suppressed.
    Notes
    -----
    The value of mu depends on the ionization fraction of the gas; mu
    = 0.59 for a fully ionized primordial gas (the default), mu = 0.61
    for a gas with ionized hydrogen but only singly ionized helium,
    and mu = 1.22 for neutral primordial gas.
    Examples
    --------
    >>> vir = find_rvT(1e12, 0)
    >>> print vir.r, vir.v, vir.T
    261.195728743 128.338776885 588643.476006
    >>> r,v,T = find_rvT(1e12, [0.5, 1, 1.5, 2])
    >>> r
    array([ 198.57846074,  156.44358398,  127.91018732,  107.74327378])
    >>> v
    array([ 147.1888328 ,  165.82959995,  183.39536858,  199.82317152])
    >>> T
    array([  774259.0063081 ,   982789.81965216,  1202024.36495813,
           1427014.0148491 ])
    """

    if isiterable(M):
        M = np.asarray(M)

    if cosmo is None:
        cosmo = cosmology.default_cosmology.get()

    # convert to cgs
    M_g = M * c.M_sun

    crit_dens = cosmo.critical_density(z)
    if astropy.__version__ >= '0.3':
        crit_dens = crit_dens.value # in ('g/cm^3')

    if r200 == True:
        rho_virial = 200.*crit_dens
    else:
        rho_virial = deltavir(z, cosmo=cosmo) * crit_dens

    # deal with cases where we have two input arrays
    if isiterable(M) and isiterable(rho_virial):
        M_g = np.atleast_2d(M_g).T
        rho_virial = np.atleast_2d(rho_virial)

    # TODO: Identify alternate Rvir definitions:
    rvir, vcirc, Tvir = _calc_rvT(M_g.to('g').value, rho_virial, mu=mu)

    return find_rvT_output(r=np.round(rvir/u.kpc.to('cm'),2)*u.kpc,
                               v=np.round(vcirc/u.km.to('cm'),2)*u.km/u.s,
                               T=np.round(Tvir,2)*u.K)
