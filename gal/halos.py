def smhm_shan(logMstar, redshift):
    """Calculate halo masses from the SMHM relations of Shan+ (2017).

    --- Inputs ---
      logMstar = log of stellar mass of galaxies (can be list/array)
      redshift = redshift of galaxies (can be list/array)

    --- Returns ---
      logM200c = M200 mass compared with critical density (numpy array)

    """

    import numpy as np

    def calc_m200(logMstar,redshift):
        """Calculate M_200c following Shan+ (2017)"""
        # Coefficients for M200c (+ scatter) from Shan+ (2017)
        if redshift < 0.2:
            logM1 = 12.52
            logM0 = 10.98
            beta = 0.47
            delta = 0.55
            gamma = 1.43
            print('Redshift z={0} is outside of the Shan+ (2017) range; '
                  'assuming value for 0.2 < z < 0.4'.format(redshift))
        elif ((redshift >= 0.2) and (redshift < 0.4)):
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

    #########################
    # Main part of the procedure just loops over the number of galaxies:
    num_gals = np.size(logMstar)

    # If we've go multiple galaxies but only one redshift, clone the redshift.
    if (num_gals != 1) & (np.size(redshift) == 1):
        zzz = np.full_like(logMstar,redshift)
        logM = np.array(logMstar)
    # If we've got multiple galaxies and redshifts
    else:
        zzz = np.array(redshift)
        logM = np.array(logMstar)

    logMhalo = []
    # Loop over the number of galaxies
    for j in np.arange(num_gals):
        logMhalo.append(calc_m200(logM[j],zzz[j]))

    # Return a numpy array:
    logMhalo = np.array(logMhalo)

    return logMhalo



def virial_radius(logMhalo,redshift,delta=200.,rhocrit=True,BryanNorman=False,WMAP=False,COSHalos=False):
    """Calculate the virial radius of a galaxy.

     --- Inputs: ---
      logMhalo: Halo mass(es).
      redshift: Galaxy redshift
      delta=200: Overdensity (Default=200)
      rhocrit: Use the critical density (Default=True); alternately the mean density of the universe.
      WMAP: Use WMAP9 cosmology. (Default=False)
      BryanNorman: Use the Bryan & Norman (1998) scaling (Default=False)
      COSHalos: Use the COS-Halos assumptions (Default=False)

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

    # Use the Bryan & Norman (2008) definition of Rvir? Default: False
    #
    # This overrides user-set flags for delta and rhocrit.
    if BryanNorman:
        # Overdensity depends on redshift:
        x = cosmo.Om(redshift)-1.
        delta = (18.*np.pi**2+82.*x-39.*x**2)

        # Also assume critical density scaling:
        rhocrit = True


    # Choose whether to scale by mean or critical density. Default: Critical
    if rhocrit == True:
        rho = cosmo.critical_density(redshift)
    else:
        rho = cosmo.Om(redshift)*cosmo.critical_density(redshift)


    # Linear halo mass (requires logMhalo in numpy array for calculations.
    Mhalo = (10.**np.array(logMhalo))*u.M_sun
    # Calculate the virial radius.
    Rvir3 = (3./(4.*np.pi))*(Mhalo/(delta*rho.to('Msun/kpc3')))
    Rvir = (Rvir3)**(1./3.)

    return Rvir

def calc_r200(logM200, redshift):
    """Convenience method for calculating R_200.
    ---Inputs---
     logM200 = Halo mass at 200x critical density (can be list/array)
     redshift = Halo redshift (can be list/array)
    """

    import numpy as np

    # Set the parameters for the virial radius calculation
    delta=200.
    rhocrit=True

    r200=virial_radius(logM200, redshift, delta=delta, rhocrit=rhocrit)

    return r200