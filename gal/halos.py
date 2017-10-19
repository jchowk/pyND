def smhm(logMstar, redshift):
    """Convenience method to calculate halo masses using default SMHM relations.

    Current default: smhm_rodriguez().

    --- Inputs ---
      logMstar = log of stellar mass of galaxies (can be list/array)
      redshift = redshift of galaxies (can be list/array)

    --- Returns ---
          logM200c = M200 mass compared with critical density (numpy array)
    """

    logMhalo = smhm_rodriguez(logMstar, redshift)

    return logMhalo


def smhm_tinker(logMstar,redshift,h=0.7):
    """Calculate Tinker halo masses for given stellar masses.

    This fit was done to Figure 10 of Tinker+ (2017) uses and assumed log M* = 0.18 dex. Below log M*=10.9 we match the power-law index from Velander+ (2014).

    Defaults to h=0.7. """

    import numpy as np
    from astropy.modeling import models

    # Tinker fit. This fit was done to Figure 10 of Tinker+ (2017) uses
    # an assumed log M* = 0.18 dex. Below log M*=10.9, we match the power-law
    # index from Velander+ (2014).
    tinker_smhm_fit = models.Chebyshev1D(5, c0=12.81806600627276,
                              c1=1.2018634412571902,
                              c2=0.013199452285390979, c3=0.01775015568831073,
                              c4=-0.029254096888480078, c5=-0.025509308396318747,
                              domain=[10.3, 12.19191])

    #########################
    # If only one redshift, duplicate it.
    num_gals = np.size(logMstar)
    if (num_gals != 1) & (np.size(redshift) == 1):
        zzz = np.full_like(logMstar,redshift)
    else:
        zzz = redshift

    # Calculate halo masses with the Tinker fit. The h^-1 correction
    #   is applied here following the discussion in their paper.
    logMhalo = np.array(tinker_smhm_fit(logMstar))-np.log10(h)

    # For stellar masses outside of the fit domain, replace masses with NaN
    bad = ((logMstar <= np.min(tinker_smhm_fit.domain)) |
            (logMstar >= np.max(tinker_smhm_fit.domain)))
    logMhalo[bad] = np.nan

    # Return a numpy array:
    logMhalo = np.array(logMhalo)

    return logMhalo


def smhm_behroozi(logMstar, redshift):
    """Calculate halo masses from the SMHM relations of Behroozi+ (2010).
       **Appropriate for redshifts z<=1.**

    --- Inputs ---
      logMstar = log of stellar mass of galaxies (can be list/array)
      redshift = redshift of galaxies (can be list/array)

    --- Returns ---
      logM200c = M200 mass compared with critical density (numpy array)

    """

    import numpy as np

    def _calc_m200_Behroozi2010(logMstar,redshift):
        """Calculate M_200c following Behroozi+ (2010)"""

        # Coefficients for M200c with redshift evolution from B10:
        M10 = 12.35
        M1a = 0.28
        M00 = 10.72
        M0a = 0.55
        beta0 = 0.44
        betaa = 0.18
        delta0 = 0.57
        deltaa = 0.17
        gamma0 = 1.56
        gammaa = 2.51

        # Scale factor (for z<1 this works fine).
        a = 1./(1.+redshift)

        logM1 = M10 + M1a*(a-1.)
        logM0 = M00 + M0a*(a-1.)
        beta  = beta0 + betaa*(a-1.)
        delta = delta0 + deltaa*(a-1.)
        gamma = gamma0 + gammaa*(a-1.)


        # Calculate the M200
        logMstarM0 = logMstar - logM0
        logM200 = (logM1 + beta * logMstarM0 +
                   10 ** (delta * logMstarM0) / (1 + 10 ** (-gamma * logMstarM0)) - 0.5)

        return logM200

    # TODO: Incorporate Behroozi+ (2013) as an option
    # def _calc_m200_Behroozi2013(logMstar,redshift):
    #     """Calculate M_200c following Behroozi+ (2013)"""
    #
    #     def _fff(x,alpha,delta,gamma)
    #         fx = -np.log10(10**(alpha*x)+1.)
    #         fx += delta*(np.log10(1.+np.exp(x)))**gamma/(1.+np.exp(10.**(-x)))
    #
    #         return fx
    #
    #
    #     # Coefficients for M200c with redshift evolution from B13:
    #     M10 = 12.35
    #     M1a = 0.28
    #     M1z =
    #
    #     delta0 = 0.57
    #     deltaa = 0.17
    #     deltaz =
    #     gamma0 = 1.56
    #     gammaa = 2.51
    #     gammaz =
    #
    #     alpha0 =
    #     alphaa =
    #     epsilon0  =
    #     epsilona  =
    #     epsilonz  =
    #     epsilona2 =
    #
    #     # Scale factor (for z<1 this works fine).
    #     z=redshift
    #     a = 1./(1.+z) # Scale factor...calculate this.
    #
    #     nu = np.exp(-4.*a**2)
    #
    #     logM1 = M10 + (M1a*(a-1.)+M1z*z)*nu
    #     logEpsilon = epsilon0+(epsilona*(a-1.)+epsilonz*z)*nu+epsilona2*(a-1.)
    #     alpha = alpha0+(alphaa*(a-1.))*nu
    #
    #     delta = delta0 + (deltaa*(a-1.)+deltaz*z)*nu
    #     gamma = gamma0 + (gammaa*(a-1.)+gammaz*z)*nu
    #
    #
    #     # Calculate the M200
    #     logMstarM0 = logMstar - logM0
    #     logM200 = (logM1 + beta * logMstarM0 +
    #                10 ** (delta * logMstarM0) / (1 + 10 ** (-gamma * logMstarM0)) - 0.5)
    #
    #     return logM200


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
    if num_gals == 1:
        logMhalo = _calc_m200_Behroozi2010(logM,zzz)
    else:
        for j in np.arange(num_gals):
            logMhalo.append(_calc_m200_Behroozi2010(logM[j],zzz[j]))

    # Return a numpy array:
    logMhalo = np.array(logMhalo)

    return logMhalo


def smhm_shan(logMstar, redshift):
    """Calculate halo masses from the SMHM relations of Shan+ (2017). At
    log M_star < 11, this relies on the results of M. Hudson+ (2015).

    --- Inputs ---
      logMstar = log of stellar mass of galaxies (can be list/array)
      redshift = redshift of galaxies (can be list/array)

    --- Returns ---
      logM200c = M200 mass compared with critical density (numpy array)

    """

    import numpy as np

    def _calc_m200(logMstar,redshift):
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
    if num_gals == 1:
        logMhalo = _calc_m200(logM,zzz)
    else:
        for j in np.arange(num_gals):
            logMhalo.append(_calc_m200(logM[j],zzz[j]))

    # Return a numpy array:
    logMhalo = np.array(logMhalo)

    return logMhalo

def smhm_rodriguez(logMstar, redshift):
        """Calculate halo masses from the SMHM relations of Rodriguez-Puebla+ (2017).

        --- Inputs ---
          logMstar = log of stellar mass of galaxies (can be list/array)
          redshift = redshift of galaxies (can be list/array)

        --- Returns ---
          logM200c = M200 mass compared with critical density (numpy array)

        """

        import numpy as np

        def _calc_mvir(logMstar,redshift):
            """Calculate M_vir following Rodriguez-Puebla+ (2017)"""
            RP17_redshifts=np.array([0.1,0.25,0.5,0.75,1.,1.25,1.5])
            logM1 = np.array([12.58,12.61,12.68,12.77,12.89,13.01,13.15])
            logM0 = np.array([10.9,10.93,10.99,11.08,11.19,11.31,11.47])
            beta = np.array([0.48,0.48,0.48,0.50,0.51,0.53,0.54])
            delta = np.array([0.29,0.27,0.23,0.18,0.12,0.03,-0.10])
            gamma = np.array([1.52,1.46,1.39,1.33,1.27,1.22,1.17])

            idx=(np.abs(redshift-RP17_redshifts)).argmin()

            # Calculate the M200
            logMstarM0 = logMstar - logM0[idx]
            logMvir = (logM1[idx] + beta[idx] * logMstarM0 +
                       10 ** (delta[idx] * logMstarM0) / (1 + 10 ** (-gamma[idx] * logMstarM0)) - 0.5)

            return logMvir

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
        if num_gals == 1:
            logMhalo = _calc_mvir(logM,zzz)
        else:
            for j in np.arange(num_gals):
                logMhalo.append(_calc_mvir(logM[j],zzz[j]))

        # Return a numpy array:
        logMhalo = np.array(logMhalo)

        return logMhalo


def virial_radius(logMhalo, redshift, delta=200., rhocrit=True,
                    BryanNorman=False, WMAP=False, COSHalos=False):
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
