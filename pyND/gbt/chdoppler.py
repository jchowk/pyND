"""
 computes the projected velocity of the telescope wrt
 six coordinate systems: geo, helio, bary, lsrk, lsrd, gal
 negative velocities mean approach

 The standard LSR is defined as follows: the sun moves at 20.0 km/s
 toward ra=18.0h, dec=30.0 deg in 1900 epoch coords

 Fully vectorized.  All three arguments must have the same dimensions.

 <p>
 This code came via e-mail from Carl Heiles via Tom Bania on 11/04/2004.
 Updated using code found via google from same source on 11/30/2009.
 Local changes:
 <UL>
 <LI> modify this documentation for use by idldoc.
 <LI> removed path argument and related code, replaced by obspos
 argument.
 <LI> default position is the GBT.
 <LI> Observatory longitude was not being passed to juldaytolmst.
 <LI> LSRD added
 <LI> Galactocentric added
 <LI> Checked against aips++ Measures.  Differences were less then 20 m/s in
 one test case (less then 10m/s for geo, bary, and lsrk).
 <LI> Double precision throughout.
 <LI> Relativistic addition of velocities.
 </UL>

 Previous revision history: carlh 29oct04
 <UL>
 <LI> from idpppler_chl; changed calculation epoch to 2000.
 <LI> 19nov04: correct bad earth spin calculation
 <LI> 7 jun 2005: vectorize to make faster for quantity calculations
 <LI> 20 Mar 2007: CH updated documentation
 <LI> 08 Feb 2015: CH fixed doppler additions for array inputs. See
       annotated statements at end of program.
 </UL>

 @param ra {in}{required} The source ra in decimal hours, equinox 2000
 @param dec {in}{required} The source dec in decimal hours, equinox 2000
 @param julday {in}{required} The julian day

 @keyword obspos {in}{optional}{type=double [2]} observatory position
 [East longitude, latitude] in degrees.
 Uses the GBT position if not specified.

 @keyword light {in}{optional}{type=boolean} When set, returns the
 velocity as a fraction of c

 @returns The velocity in km/s, or as a faction of c if
 the keyword /light is specified. the result is a 6-element
 vector whose elements are [geo, helio, bary, lsrk, lsrd, gal].

 @uses <a href="http://idlastro.gsfc.nasa.gov/ftp/pro/astro/baryvel.pro">baryvel</a>
 @uses <a href="http://idlastro.gsfc.nasa.gov/ftp/pro/astro/precess.pro">precess</a>
 @uses <a href="juldaytolmst.html">juldaytolmst</a>

 @version $Id: chdoppler.pro,v 1.7 2009/12/01 17:33:37 bgarwood Exp $
"""

from __future__ import print_function

from numpy import *

def chdoppler(ra, dec, julday, obspos=None, light=None):
#
# Default to GBT if obspos not provided.
    n_params = 3
    def _ret():  return None
    
    if (bitwise_not((obspos is not None))):    
        obspos = zeros([2], "float64")#
        obspos[0] = array([[-(79.e0 + 50.e0 / 60.e0 + 23.3988e0 / 3600.e0)]])
        obspos[1] = array([[(38.e0 + 25.e0 / 60.e0 + 59.2284e0 / 3600.e0)]])
    
    #------------------ORBITAL SECTION-------------------------
    nin = ra.size
    
    #GET THE COMPONENTS OF RA AND DEC, 2000u EPOCH
    rasource = ra * 15.e0 * _sys_dtor
    decsource = dec * _sys_dtor
    
    xxsource = zeros([nin, 3], "float64")
    xxsource[:,0] = cos(decsource) * cos(rasource)
    xxsource[:,1] = cos(decsource) * sin(rasource)
    xxsource[:,2] = sin(decsource)
    pvorbit_helio = zeros([nin], "float64")
    pvorbit_bary = zeros([nin], "float64")
    pvlsrk = zeros([nin], "float64")
    pvlsrd = zeros([nin], "float64")
    pvgal = zeros([nin], "float64")
    
    #GET THE EARTH VELOCITY WRT THE SUN CENTER
    #THEM MULTIPLY BY SSSOURCE TO GET $
    #       PROJECTED VELOCITY OF EARTH CENTER WRT SUN TO THE SOURCE
    for NR in arange(0, (NIN - 1)+(1)):
        baryvel(julday[nr], 2000., vvorbit, velb)
        pvorbit_helio[nr] = total(vvorbit * xxsource[nr,:])
        pvorbit_bary[nr] = total(velb * xxsource[nr,:])
    
    #stop
    
    #-----------------------LSRK SECTION-------------------------
    #THE STANDARD LSRK IS DEFINED AS FOLLOWS: THE SUN MOVES AT 20.0 KM/S
    #TOWARD RA=18.0H, DEC=30.0 DEG IN 1900 EPOCH COORDS
    #using PRECESS, this works out to ra=18.063955 dec=30.004661 in 2000
    #coords.
    ralsrk_rad = 2.e0 * _sys_pi * 18.e0 / 24.e0
    declsrk_rad = _sys_dtor * 30.e0
    precess(ralsrk_rad, declsrk_rad, 1900.e0, 2000.e0, radian=True)
    
    #FIND THE COMPONENTS OF THE VELOCITY OF THE SUN WRT THE LSRK FRAME
    xxlsrk = zeros([nin, 3], "float64")
    xxlsrk[:,0] = cos(declsrk_rad) * cos(ralsrk_rad)
    xxlsrk[:,1] = cos(declsrk_rad) * sin(ralsrk_rad)
    xxlsrk[:,2] = sin(declsrk_rad)
    vvlsrk = 20.e0 * xxlsrk
    
    #PROJECTED VELOCITY OF THE SUN WRT LSRK TO THE SOURCE
    for nr in arange(0, (nin - 1)+(1)):
        pvlsrk[nr] = total(vvlsrk * xxsource[nr,:])
    
    
    #-----------------------LSRD SECTION-------------------------
    #THE LSRD IS DEFINED AS FOLLOWS: THE SUN MOVES AT 16.6 KM/S
    #TOWARD RA=17:49:58.7 hours, DEC=28.07.04 DEG IN 2000 EPOCH COORDS
    
    ralsrd_rad = 2.e0 * pi * (17.e0 + 49.e0 / 60.e0 + 58.7e0 / 3600.e0) / 24.e0
    declsrd_rad = _sys_dtor * (28.e0 + 07.e0 / 60.e0 + 04.0e0 / 3600.e0)
    
    #FIND THE COMPONENTS OF THE VELOCITY OF THE SUN WRT THE LSRD FRAME
    xxlsrd = zeros([nin, 3], "float64")
    xxlsrd[:,0] = cos(declsrd_rad) * cos(ralsrd_rad)
    xxlsrd[:,1] = cos(declsrd_rad) * sin(ralsrd_rad)
    xxlsrd[:,2] = sin(declsrd_rad)
    vvlsrd = 16.6e0 * xxlsrd
    
    #PROJECTED VELOCITY OF THE SUN WRT LSRD TO THE SOURCE
    for nr in arange(0, (nin - 1)+(1)):
        pvlsrd[nr] = total(vvlsrd * xxsource[nr,:])
    
    #-----------------------GALACTOCENTRIC SECTION------------------
    # LSRD + 220 km/s towards RA=21:12:01.1 DEC=48.19.47 in J2000 Epoch
    
    ragal_rad = 2.e0 * pi * (21.e0 + 12.e0 / 60.e0 + 01.1e0 / 3600.e0) / 24.e0
    decgal_rad = _sys_dtor * (48.e0 + 19.e0 / 60.e0 + 47.e0 / 3600.e0)
    
    # Find the components of the velocity of the sun wrt this frame
    xxgal = zeros([nin, 3], "float64")
    xxgal[:,0] = cos(decgal_rad) * cos(ragal_rad)
    xxgal[:,1] = cos(decgal_rad) * sin(ragal_rad)
    xxgal[:,2] = sin(decgal_rad)
    vvgal = 220.e0 * xxgal
    
    #PROJECTED VELOCITY OF THE SUN WRT GAL TO THE SOURCE
    for nr in arange(0, (nin - 1)+(1)):
        pvgal[nr] = total(vvgal * xxsource[nr,:])
    
    #---------------------EARTH SPIN SECTION------------------------
    lat = obspos[1]
    
    lst_mean = 24.e0 / (2.e0 * _sys_pi) * juldaytolmst(julday, obslong=obspos[0])
    
    #MODIFIED EARTH SPIN FROM GREEN PAGE 270
    pvspin = -0.465 * cos(_sys_dtor * lat) * cos(decsource) * sin((lst_mean - ra) * 15. * _sys_dtor)
    
    #---------------------NOW PUT IT ALL TOGETHER------------------
    
    vtotal = zeros([nin, 6], "float64")
    vtotal[:,0] = -pvspin
    vtotal[:,1] = shiftvel(-pvspin, -pvorbit_helio)
    vtotal[:,2] = shiftvel(-pvspin, -pvorbit_bary)
    vtotal[:,3] = shiftvel(vtotal[:,2], -pvlsrk)
    vtotal[:,4] = shiftvel(vtotal[:,2], -pvlsrd)
    vtotal[:,5] = shiftvel(vtotal[:,4], -pvgal)
    
    #the statements below are wrong (CH 8 feb 2015)
    #vtotal[ 3,*]= -shiftvel(-vtotal[2],pvlsrk)
    #vtotal[ 4,*]= -shiftvel(-vtotal[2],pvlsrd)
    #vtotal[ 5,*]= -shiftvel(-vtotal[4],pvgal)
    
    if (light is not None):    
        vtotal = vtotal / (_sys_gc.light_speed * 1.e3)
    
    #stop
    
    return vtotal

