from numpy import *
from astropy.coordinates import SkyCoord

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

    # Radio definition coordinates
    llsr = 56.
    blsr = 22.
    lsr_coords = SkyCoord(llsr,blsr,frame='galactic',unit='deg')
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
    vlsr_out = vlsr * cos(dlsr)
    vmb_out = vmb * cos(dmb)
    
    
    if silent == False:
        ## Print output...
        print( "\n LSR Correction: ")
        print( "     LSR     = {0:0.2f} km/s".format(vlsr_out))
        print( "     LSR(MB) = {0:0.2f} km/s.".format(vmb_out))
        print( "\n v(LSR) = v(helio) + LSR \n")


    def _ret():
        if mihalas == True:
            deltav_out = vmb_out
        else:
            deltav_out = vlsr_out
        return deltav_out

    return _ret()

