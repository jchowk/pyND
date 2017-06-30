def ccm_filters(av=1, rv=3.1, showtable=False):
    """output = ccmfilters(av=1,rv=3.1,showtable=False)

    Calculate A_lambda/A_V for broadband filters given an R_V. Return extinctions for an assumed Av. Assumes CCM extinction curve.

    :param av: Input av [default: av=1]         
    :param rv: Input ratio of total to selective exinction [default: rv=3.1]
    :param showtable: Print the result? [default: False] 
    :return: Returns an astropy table of extinctions in Johnson-Cousins bands U-L. 
    """

    import numpy as np
    from astropy.table import Table
    import pdb

    #     Set up the CCM information.
    filters = np.array(['U', 'B', 'Hb', 'V', 'Ha','R', 'I', 'J', 'H', 'K', 'L'])
    ccm_a = np.array([0.953, 0.9982, 1.0118, 1.000, 1.0211, 0.8686, 0.680, 0.4008, 0.2694, 0.1615, 0.08])
    ccm_b = np.array([1.909, 1.0495, 0.4613, 0, -0.2825, -0.3660, -0.6239, -0.3679, -0.2463, -0.1483, -0.0734])

    extinct = (ccm_a + ccm_b / rv) * av
    extinct_table = Table(np.around(extinct, 3), names=filters)

    if showtable:
        print(extinct_table)

    def _ret():  return extinct_table
    return _ret()

def ccm_extinct(lam_in,av=1, rv=3.1):
    """Calculate A_lambda for lambda in Angstroms. 
    (Default for Av=1, Rv=3.1.)"""
    import astropy.units as u
    import numpy as np
    import pdb

    def a_opt(x):
        y = (x-1.82)
        ax = 1.+0.17699*y ** 2-0.02427*y ** 3+0.72085*y ** 4+0.01979*y ** 5-0.77530*y ** 6+0.32999*y ** 7
        return ax

    def b_opt(x):
        y = (x-1.82)
        bx = 1.41338*y+2.28305*y ** 2+1.07233*y ** 3-5.38434*y ** 4-0.62251*y ** 5+5.30260*y ** 6-2.09002*y ** 7
        return bx

    def a_nuv(x):
        if (x >=5.9):
            fa = -0.4473*(x-5.9)**2.-0.009779*(x-5.9)**3.
        else:
            fa=0.

        ax = 1.752-0.316*x-0.104 / ((x-4.67)**2.+0.341)+fa
        return ax
    def b_nuv(x):
        if (x >=5.9):
            fb = 0.2130*(x-5.9)**2.+0.1207*(x-5.9)**3
        else:
            fb = 0.
        bx = -3.090+1.825*x+1.206 / ((x-4.62)**2.+0.263)+fb
        return bx

    def a_fuv(x):
        ax = -1.073-0.628*(x-8.)+0.137*(x-8.)**2.-0.070*(x-8.)**3
        return ax
    def b_fuv(x):
        bx = 13.67+4.257*(x-8.)-0.42*(x-8.)**2.+0.374*(x-8.)**3.
        return bx

    def extinct_calc(x,ax,bx,Rv=3.1):
        out = ax + bx / Rv
        return out

    # Make sure our input is a proper array:
    lam_in = np.array(lam_in)

    # Convert input to microns
    lam_micron = lam_in / 1.e4
    x = (1./lam_micron)

    extinct_out = np.zeros_like(lam_micron)
    if np.size(x) != 1:
        for j in np.arange(np.size(x)):
            if  ((x[j] <= 3.3) & (x[j] >= 0.9)):
                aaa=a_opt(x[j])
                bbb=b_opt(x[j])
            elif ((x[j] > 3.3) & (x[j] <=8.0)):
                aaa = a_nuv(x[j])
                bbb = b_nuv(x[j])
            elif ((x[j] > 8.0) & (x[j] <= 11.0)):
                aaa = a_fuv(x[j])
                bbb = b_fuv(x[j])
            extinct_out[j] = av*extinct_calc(x[j],aaa,bbb,rv)
    else:
        if ((x <= 3.3) & (x >= 0.9)):
            aaa = a_opt(x)
            bbb = b_opt(x)
        elif ((x > 3.3) & (x <= 8.0)):
            aaa = a_nuv(x)
            bbb = b_nuv(x)
        elif ((x > 8.0) & (x <= 11.0)):
            aaa = a_fuv(x)
            bbb = b_fuv(x)
        extinct_out = av * extinct_calc(x, aaa, bbb, rv)

    return extinct_out



def lsrvel(long, lat, radec=False, mihalas=False, silent=False):
    """delta_v = lsrvel(long, lat, mihalas=False, SILENT=False):

       This program calculates the projection of the velocity vector of the local standard of 
       rest on the sightline specified by (l,b)=(long,lat) used to calculate the shift from 
       heliocentric to LSR velocities: v(LSR) = v(helio) + LSR

          Assumes v(LSR) = 20   km/sec to (l,b)=(56, 22) or
                  v(LSR) = 16.5 km/sec to (l,b)=(53, 25)
                                  from Mihalas & Binney

          Created by JCH 9/27/99

    :param long: Longitude [Galactic unless radec=True]
    :param lat: Latitude [Galactic unless radec=True]
    :param radec: input coordinates are RA/Dec? [default: False]
    :param mihalas: use the Mihalas & Binney definition of LSR [default: False]
    :param SILENT: suppress printing (default: False)
    :return: LSR shift, where v(LSR) = v(helio) + LSR 
    """
    import numpy as np
    from astropy.coordinates import SkyCoord

    # Radio definition coordinates
    llsr = 56.
    blsr = 22.
    lsr_coords = SkyCoord(llsr, blsr, frame='galactic', unit='deg')
    # Radio definition velocity
    vlsr = 20.0

    # M&B coordinates
    lmb = 53.
    bmb = 25.
    mb_coords = SkyCoord(lmb, bmb, frame='galactic', unit='deg')
    # M&B velocity
    vmb = 16.5

    if radec == False:
        input_coords = SkyCoord(long, lat, frame='galactic', unit='deg')
    else:
        input_coords = SkyCoord(long, lat, frame='icrs', unit='deg')

    # Calculate the separations on the sky [given in degrees]:
    dlsr = input_coords.separation(lsr_coords)
    dmb = input_coords.separation(mb_coords)

    # Calculate the projected velocities:
    vlsr_out = vlsr * np.cos(dlsr)
    vmb_out = vmb * np.cos(dmb)

    if silent == False:
        ## Print output...
        print("\n LSR Correction: ")
        print("     LSR     = {0:0.2f} km/s".format(vlsr_out))
        print("     LSR(MB) = {0:0.2f} km/s.".format(vmb_out))
        print("\n v(LSR) = v(helio) + LSR \n")

    def _ret():
        if mihalas == True:
            deltav_out = vmb_out
        else:
            deltav_out = vlsr_out
        return deltav_out

    return _ret()


def rotcurve(long, lat, distance_input=[-99], do_plot=False, tenkpc=False, constant=False, radec=False):
    """
    out = rotcurve(long, lat, distance_input=[...],do_plot=False, tenkpc=False, constant=False,radec=False):

      Returns a Clements (1985) rotation curve for a given direction.  Optionally plots the result. 

    :param long: Longitude coordinate [assumed Galactic w/o radec=True]
    :param lat: Latitude coordinate [assumed Galactic w/o radec=True] 
    :param distance_input: [opt] An np.array holding distances at which to calculate the rotation curve. 
    :param do_plot: Show the results.
    :param tenkpc: Use parameters for a 10 kpc solar circle.
    :param constant: Assume a perfectly flat rotation curve after ~3 kpc.
    :param radec: Input coordinates are RA/Dec. [default: False]
    :return: rotation curve velocity [, distance np.array] 
    """

    import numpy as np
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import matplotlib.pyplot as plt

    # Test for distance np.array:
    if distance_input[0] == -99:
        # If no input np.array, define the np.array to 15 kpc
        distance_array = np.arange(0, 15, 0.1, dtype=np.float64)
        distance_output = True
    else:
        # If there is an input, make sure it doesn't cause troubles with type mismatches.
        distance_array = np.array(distance_input, dtype=np.float64)
        distance_output = False

    # Deal with the coordinates:
    if radec:
        coords = SkyCoord(long, lat, frame='icrs', unit='deg')
    else:
        coords = SkyCoord(long, lat, frame='galactic', unit='deg')

    if tenkpc:
        solarcircle = 10.
        solarvel = 250.00

        ##From Clemens 1985
        Dcoeff = np.array([264.76], dtype=np.float64)
        Ccoeff = np.array(
            [-1931.3363, 1768.47176, -606.461977, 112.102138,
             -11.9698505, +0.7367828, -0.02423453, +0.00032952], dtype=np.float64)
        Bcoeff = np.array([319.8354, -194.3591, +155.30181, -63.00332, 12.142556, -0.869573], dtype=np.float64)
    else:
        solarcircle = 8.5
        solarvel = 220.00

        ##From Clemens 1985
        Dcoeff = np.array([234.88], dtype=np.float64)
        Ccoeff = np.array(
            [-2342.6564, +2507.6039, -1024.06876, +224.562732,
             -28.4080026, +2.0697271, -0.08050808, +0.00129348], dtype=np.float64)
        Bcoeff = np.array([325.0912, -248.1467, 231.87099, -110.73531, +25.073006, -2.110625], dtype=np.float64)

    galcenter_dist = solarcircle * np.sqrt(
        (distance_array / solarcircle) ** 2.-2. * (distance_array / solarcircle) * np.cos(coords.galactic.l) + 1.)
    vel_rot = np.zeros_like(galcenter_dist)

    goodB = np.where(galcenter_dist / solarcircle < 0.45)
    goodC = np.where((galcenter_dist / solarcircle >= 0.45) & (galcenter_dist / solarcircle <= 1.6))
    goodD = np.where(galcenter_dist / solarcircle > 1.6)

    # B coefficient for R/Rsun < 0.45
    if np.size(goodB) != 0:
        for j in np.arange(np.size(Bcoeff)):
            vel_rot[goodB] = vel_rot[goodB] + Bcoeff[j] * galcenter_dist[goodB] ** j

    # C Coefficients for 0.45 < R/Rsun < 1.6
    for i in np.arange(np.size(Ccoeff)):
        vel_rot[goodC] = vel_rot[goodC] + Ccoeff[i]*galcenter_dist[goodC]**i

    # D coefficient for R/Rsun > 1.6
    if np.size(goodD) != 0:
        vel_rot[goodD] = Dcoeff[0]

    if constant:
        vel_rot = solarvel

    vel_radial = (vel_rot / galcenter_dist - solarvel / solarcircle) * solarcircle * np.sin(coords.galactic.l) * np.cos(
        np.absolute(coords.galactic.b))

    # Apply units. This could be a mistake.
    # vel_radial = vel_radial * u.km / u.s
    # distance_array = distance_array * u.kpc

    # Do the plots?
    if do_plot:
        plt.plot(distance_array, vel_radial)
        plt.xlabel('Distance [kpc]')
        plt.ylabel('Velocity [km/s]')
        plt.minorticks_on()

    def _ret():
        if distance_output == True:
            return vel_radial, distance_array
        else:
            return vel_radial

    return _ret()


def kinedist(long, lat, vel, rotcurve=False, do_plot=False, tenkpc=False, constant=False, radec=False):


    """DO NOT USE...not finished.
    """

    import numpy as np
    import matplotlib.pyplot as plt

    # Create a fine grid of distances.
    distance = np.arange(0.01, 15.01, 0.01)
    vel_radial = rotcurve(long, lat, distance_input=distance,
                              tenkpc=tenkpc, constant=constant, radec=radec)


    # Function to extract indeces for velocities "close" to our target velocity.
    def find_velmatches(vel_array, vel_target):
        return np.where(np.abs(vel_array - vel_target) <= 1.)


    # Here are the velocity indeces that match.
    close_velocities = find_velmatches(vel_radial, vel)

    #  TODO Need to decide how we find the "best" distance.
    # Our rotation curve could be double-valued, so we need to find first section that's consistent, then think
    # about the range of acceptable distances.


    bestdist = dist[index] / np.cos(np.absolute(coords.galactic.b))

    print('Distance: {0:f0.2} kpc'.format(bestdist))

    if do_plot:
        if bestdist > 5.:
            upperx = 1.5 * bestdist
        else:
            upperx = 5.

            plt.plot(dist / np.cos(np.absolute(coords.galactic.b)), vel_radial)
            plt.xlim([0, upperx])
            plt.xlabel('Distance (kpc)')
            plt.ylabel('Velocity (km/sec)')

            plt.plot([bestdist], [vel], 'ro')

            rotdist = dist / np.cos(np.absolute(coords.galactic.b))
            rotcurve = vel_radial


    def _ret():
        if rotcurve:
            return rotdist, rotcurve
        else:
            return rotdist


    return _ret()
