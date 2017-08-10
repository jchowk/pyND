def smhm_shan(logMstar, redshift):
    """Halo masses from the SMHM relations of Shan+ (2017)."""

    import numpy as np

    # Coefficients for M200c (+ scatter)
    if redshift < 0.2:
        logM1 = 12.52
        logM0 = 10.98
        beta = 0.47
        delta = 0.55
        gamma = 1.43
        print('Redshift z={0} is outside of the Shan+ (2017) range; '
              'assuming value for 0.2 < z < 0.4'.format(redshift))
    elif ((redshift > 0.2) and (redshift < 0.4)):
        logM1 = 12.52
        logM0 = 10.98
        beta = 0.47
        delta = 0.55
        gamma = 1.43
    elif ((redshift >= 0.4) & (redshift <= 0.6)):
        logM1 = 12.70
        logM0 = 11.11
        beta = 0.50
        delta = 0.54
        gamma = 1.72
    elif ((redshift > 0.6) & (redshift <= 1.0)):
        logM1 = 12.70
        logM0 = 11.11
        beta = 0.50
        delta = 0.54
        gamma = 1.72
        print('Redshift z={0} is outside of the Shan+ (2017) range; '
              'assuming value for 0.4 < z < 0.6'.format(redshift))

    # Calculate the M200
    logMstarM0 = logMstar - logM0
    logM200 = (logM1 + beta * logMstarM0 +
               10 ** (delta * logMstarM0) / (1 + 10 ** (-gamma * logMstarM0)) - 0.5)

    return logM200

def virial_radius(logMhalo,redshift,delta=200.,rhocrit=True,BryanNorman=False,WMAP=False,COSHalos=False):
    """Calculate the virial radius of a galaxy.

    :param logMhalo: Halo mass
    :param redshift: Galaxy redshift
    :param delta: Overdensity (Default=200)
    :param rhocrit: Use the critical density (Default=True); alternately the mean density.
    :param WMAP: Use WMAP9 cosmology. (Default=False)
    :param BryanNorman: Use the Bryan & Norman (1998) scaling (Default=False)
    :param COSHalos: Use the COS-Halos assumptions (Default=False)
    """

    import numpy as np
    import astropy.units as u
    import astropy.constants as c

    # Set some terms based on the COSHalos flag, which matches the calculations
    # from the COS-Halos work (Prochaska+ 2017).
    #
    # This overrides the other user-set flags.
    if COSHalos:
        BryanNorman=True
        WMAP=True
        rhocrit=True

    # Choose the cosmology. Default is Plank15.
    if WMAP:
        from astropy.cosmology import WMAP9 as cosmo
    else:
        from astropy.cosmology import Planck15 as cosmo

    # Choose whether to scale by mean or critical density. Default: Critical
    if rhocrit == True:
        rho = cosmo.critical_density(redshift)
    else:
        rho = cosmo.Om(redshift)*cosmo.critical_density(redshift)

    # Use the Bryan & Norman (2008) definition of Rvir? Default: False
    if BryanNorman:
        # Overdensity depends on redshift:
        x = cosmo.Om(redshift)-1.
        delta = 18.*np.pi**2+82.*x-39.*x**2

        # Also assume critical density scaling:
        rho = cosmo.critical_density(redshift)

    Mhalo = (10.**logMhalo)*u.M_sun
    Rvir3 = (3./(4.*np.pi))*(Mhalo/(delta*rho.to('Msun/kpc3')))
    Rvir = (Rvir3)**(1./3.)

    return Rvir
