import numpy as np
from astropy.io import fits,ascii
import matplotlib.pyplot as plt

from astropy.table import Table
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.constants as c
import astropy.units as u
import os

# TODO: Work masks into resample method.

# TODO: Change OPTICAL to RADIO LSR
# TODO: Heliocentric to RADIO
# TODO: Estimate RMS / stats
# TODO: fit Gaussian?

class GBTspec(object):

    @classmethod
    def from_ascii(cls,input_filename,
                    format='no_header',
                    header_start=1,
                    data_start=3,
                    **kwargs):

        # Read the data as an astropy table
        a = ascii.read(input_filename,
            format=format,
            header_start=header_start,
            data_start=data_start,
            **kwargs)

        # Temperature system: antenna or brightness?
        if a.colnames[1] != 'Tb':
            efficiency_correction = 1./0.88 # Main beam efficiency
        else:
            efficiency_correction = 1.

        # Fill the spectral information
        velocity = np.array(a[a.colnames[0]])
        Tb = np.array(a[a.colnames[1]])*efficiency_correction

        # Initiate
        slf = cls(velocity, Tb)

        # META DATA:
        # Fill the information about the object:
        slf.meta['filename'] = os.path.abspath(input_filename)

        # Read the file as a whole to grab information.
        #  [Assumes standard GBTIDL ASCII output.]
        rd = open (slf.meta['filename'], "r")

        # Read list of lines
        out = rd.readlines()
        rd.close()

        # Parse the object name.
        object_name = out[0].split()[2]

        # Set the velocity frame.
        frm = out[1].split()[0]
        vsys = out[2].split()[0].split('-')[1]
        slf.meta['veldef'] = frm+'-'+vsys

        # Apply the object name
        slf.meta['object'] = object_name

        simbad_result = Simbad.query_object(slf.meta['object'])
        if isinstance(simbad_result,Table):
            coords = SkyCoord(simbad_result['RA'],simbad_result['DEC'],
                  unit=(u.hourangle, u.deg))
            slf.meta['RA'] = coords.ra.deg[0]
            slf.meta['DEC'] = coords.dec.deg[0]

            slf.meta['l'] = coords.galactic.l.deg
            slf.meta['b'] = coords.galactic.b.deg

        return slf

    @classmethod
    def from_GBTIDL(cls,input_filename,object_name=None):

        # Load the GBTIDL data:
        a = fits.open(input_filename)

        if object_name == None:
            object_name = input("Object to extract: ")

        for j in np.arange(1,np.size(a)):
            xxx = a[j].data
            gd=(xxx['OBJECT'] == object_name)
            if gd.sum() > 0:
                b = xxx[gd]

        if ('b' in locals()) == False:
            print('No matching object found.')
            return

        # Define the spectrum:
        nu0=b['RESTFREQ']
        nu = ((np.arange(np.size(b['DATA']))+1)-b['CRPIX1'])*b['CDELT1'] + b['CRVAL1']

        velocity = (nu0-nu)/nu0 * c.c.to('km/s').value
        Tb = b['DATA'][0]
        if b['TUNIT7'] == 'Ta*':
            Tb /= 0.88 # Main beam efficiency correction

        # Initiate
        slf = cls(velocity, Tb)

        # META DATA
        # Fill the information about the object:
        slf.meta['filename'] = os.path.abspath(input_filename)
        slf.meta['object'] = object_name

        # Fill in some data/information from the GBTIDL format:
        slf.meta['RA'] = b['TRGTLONG'][0]
        slf.meta['DEC'] = b['TRGTLAT'][0]
        # slf.meta['RA'] = b['TRGTLONG']
        # slf.meta['DEC'] = b['TRGTLAT']

        # Define the Galactic coordinates
        slf._fill_Galactic_coords()

        # Details of the GBT data
        slf.meta['veldef'] = b['VELDEF'][0]
        slf.meta['restfreq'] = b['RESTFREQ'][0]
        slf.Tsys = b['TSYS'][0]

        return slf

    @classmethod
    def from_GBTIDLindex(cls,input_filename,
                        object_indx=None):

        # Load the GBTIDL data:
        a = fits.open(input_filename)

        if object_indx == None:
            tbl=GBTspec.index_GBTIDL(input_filename,
                            return_list=True)

            valid = False
            while valid == False:
                try:
                    object_indx = np.int(input("Index of object to load: "))
                except:
                    # Not a numerical entry.
                    print("Enter numerical value from list of indeces.")
                if (object_indx <= np.max(tbl[tbl.colnames[0]])) & (object_indx >= 0):
                    valid = True
                    break
                else:
                    print("Enter an index within the range listed.")



            # What array indeces
            array_indx = tbl['ARRAY_INDECES'][object_indx]
        else:
            tbl=GBTspec.index_GBTIDL(input_filename,
                                        silent=True,
                                        return_list=True)
            array_indx = tbl['ARRAY_INDECES'][object_indx]


        # Set the output object:
        b = a[array_indx[0]].data[array_indx[1]]

        # Define the spectrum:
        nu0=b['RESTFREQ']
        nu = ((np.arange(np.size(b['DATA']))+1)-b['CRPIX1'])*b['CDELT1'] + b['CRVAL1']

        velocity = (nu0-nu)/nu0 * c.c.to('km/s').value
        Tb = b['DATA']
        if b['TUNIT7'] == 'Ta*':
            Tb /= 0.88 # Main beam efficiency correction

        # Initiate
        slf = cls(velocity, Tb)

        # META DATA
        slf.meta['filename'] = os.path.abspath(input_filename)

        # Fill in some data/information from the GBTIDL format:
        slf.meta['object'] = b['OBJECT']
        slf.meta['RA'] = b['TRGTLONG']
        slf.meta['DEC'] = b['TRGTLAT']

        # Define the Galactic coordinates
        slf._fill_Galactic_coords()

        # Details of the GBT data
        slf.meta['veldef'] = b['VELDEF']
        slf.meta['restfreq'] = b['RESTFREQ']

        slf.Tsys = b['TSYS']

        return slf

    def __init__(self, velocity, Tb, mask=None, meta=None):

        # The velocity and brightness temperatures
        self.velocity = velocity.copy()
        self.Tb = Tb.copy()

        if mask is None:
            self.mask = np.repeat(False,len(Tb))
        else:
            self.mask = mask.copy()

        if meta is None:
            self.meta = {}
            # Where are the data?
            self.meta['filename'] = None

            # Information from the GBTIDL structure
            self.meta['object'] = None
            self.meta['RA'] = None
            self.meta['DEC'] = None

            self.meta['veldef'] = None
            self.meta['restfreq'] = None

            self.meta['Tsys'] = None

        else:
            self.meta = meta.copy()



    def _fill_Galactic_coords(self):

            coords = SkyCoord(self.meta['RA'],self.meta['DEC'],unit=u.deg)

            self.meta['l'] = coords.galactic.l.deg
            self.meta['b'] = coords.galactic.b.deg




    def plotspectrum(self,**kwargs):
        try:
            kwargs['color']
        except:
            kwargs['color']='green'

        try:
            ylim=kwargs.pop('ylim')
        except:
            ylim=None

        try:
            xlim=kwargs.pop('xlim')
        except:
            xlim=None

        plt.figure(figsize=(10,5))

        plt.plot(self.velocity,self.Tb,drawstyle='steps-mid',
                    **kwargs)
        plt.axhline(0.,linestyle='--',linewidth=1,color='k',zorder=0)

        # if self.meta['veldef'].find('LSR') > 0:
        #     xlbl_text = 'LSR Velocity [km/s]'
        # elif self.meta['veldef'].find('HELI') > 0:
        #     xlbl_text = 'Heliocentric Velocity [km/s]'
        # else:
        #     xlbl_text = '{0} Velocity [km/s]'.format(self.meta['veldef'])

        xlbl_text = '{0} Velocity [km/s]'.format(self.meta['veldef'])

        plt.xlabel(xlbl_text)
        plt.ylabel('$T_b$ [K]')
        plt.title(self.meta['object'])

        ylim=plt.ylim(ylim)
        xlim=plt.xlim(xlim)


    def index_GBTIDL(input_filename, silent=False,
                        return_list=False):
        """Return an index of information in the GBTIDL FITS file."""
        # Load the GBTIDL data:
        a = fits.open(input_filename)

        # Create the index
        # Counter variable
        i=0
        # Outputs
        indx=[]
        object_names = []
        array_indeces = []

        for j in np.arange(1,np.size(a)):
            xxx = a[j].data
            for k in np.arange(np.size(xxx)):
                indx.append(i)
                object_names.append(xxx['OBJECT'][k])
                array_indeces.append((j,k))

                if silent == False:
                    # print("{0}: {1:20s}   \t{2}".format(i,
                    #     object_names[i],array_indeces[i]))
                    print("{0}: {1:20s}".format(i,
                        object_names[i]))

                # Advance the counter
                i+=1

        if return_list:
            tbl = Table([indx,object_names,array_indeces],
                        names=['INDX','OBJECT','ARRAY_INDECES'])
            return tbl

    def change_veldef(self):
        # Equivalent to GBTDL code `veltovel`, but limited to RADIO/OPTICAL LSR frames.

        light_speed = np.double(c.c.to('m/s').value)
        nu0 = np.double(1420405800.0000000000000000000)

        if not self.meta['restfreq']:
            self.meta['restfreq'] = nu0

        if self.meta['veldef'] == 'OPTICAL-LSR':
            # Change OPTICAL-->RADIO definition of LSR

            # First calculate "TRUE" velocity
            result = self.velocity*1000. / light_speed
            result = (2.e0 * result + result * result) / (2.e0 + 2.e0 * result + result * result)

            # Calculate RADIO definition
            result = 1.e0 - np.sqrt((1.e0 - result) / (1.e0 + result))

            # Update:
            self.velocity = result*light_speed/1000.
            self.meta['veldef'] = 'RADI-LSR'

        elif self.meta['veldef'] == 'RADI-LSR':
            # Change RADIO-->OPTICAL definition of LSR

            # First calculate "TRUE" velocity
            result = self.velocity*1000. / light_speed
            result = (2.e0 * result - result * result) / (2.e0 - 2.e0 * result + result * result)

            # Calculate OPTICAL definition
            result = np.sqrt((1.e0 + result) / (1.e0 - result)) - 1.e0

            # Update veldef:
            self.velocity = result*light_speed/1000.
            self.meta['veldef'] = 'OPTICAL-LSR'
        else:
            print('Inappropriate velocity definition.')


    def columndensity(self,vel_range=[-100,100.]):
        scalefactor = 1.823e18

        npix = len(self.velocity)
        vlh = (self.velocity + np.roll(self.velocity, -1)) / 2.
        vlh[npix - 1] = self.velocity[npix - 1] + \
                        (self.velocity[npix - 1] - self.velocity[npix - 2]) / 2.
        dvl = vlh - np.roll(vlh, 1)
        dvl[0] = 2 * (vlh[0] - self.velocity[0])
        med_dvl = np.median(dvl)

        # Find nearest pixels.
        # For now the integration is done only over the nearest central velocities.
        vl_limits = np.searchsorted(self.velocity,vel_range)

        # Cumulative Sum
        sum = np.sum(self.Tb[vl_limits[0]:vl_limits[1]] * \
                                dvl[vl_limits[0]:vl_limits[1]])

        HIcolumn = scalefactor * sum

        return HIcolumn

    def resample(self, new_velocity,
                fill_value=None, **kwargs):
        """ ADAPTED FROM LINETOOLS REBIN CODE. [https://github.com/linetools/linetools]

        Resample a single spectrum to a new velocity array.

        Uses simple linear interpolation.  The default (and only)
        option conserves counts (and flambda).

        WARNING: Do not trust either edge pixel of the new array.

        Also be aware that neighboring pixels are likely to be
        correlated in a manner that is not described by the error
        array.

        Parameters
        ----------
        new_velocity : array
          New velocity array
        fill_value : float, optional
          Fill value at the edges
          Default = None [filling with NAN]
        """
        from scipy.interpolate import interp1d
        import warnings

        flux = self.Tb
        velocity = self.velocity
        mask = self.mask

        # Deal with nan
        badf = np.any([np.isnan(flux), np.isinf(flux)], axis=0)
        if np.sum(badf) > 0:
            warnings.warn("Ignoring pixels with NAN or INF in flux")
        # Select the good data
        gdf = ~badf
        flux = flux[gdf]

        # Endpoints of original pixels
        npix = len(velocity)

        # Average velocity positions
        vlh = (velocity + np.roll(velocity, -1)) / 2.
        vlh[npix - 1] = velocity[npix - 1] + \
                        (velocity[npix - 1] - velocity[npix - 2]) / 2.
        # Delta velocity
        dvl = vlh - np.roll(vlh, 1)
        dvl[0] = 2 * (vlh[0] - velocity[0])

        # Select "good" data points – those not NAN or INF
        vlh = vlh[gdf]
        dvl = dvl[gdf]

        # To conserve flux, use the cumulative sum as a function of velocity as
        # the basis for the interpolation.

        # Cumulative sum of the brightness temperatures
        cumsum = np.cumsum(flux * dvl)

        # Interpolate the cumulative sum FILL_VALUE should probably be 0.
        fcum = interp1d(vlh, cumsum, fill_value=0., bounds_error=False)

        # Create a reference interpolation to fill/flag pixels outside the range
        # of the original data.
        fcum_ref = interp1d(vlh, cumsum, fill_value=fill_value, bounds_error=False)

        # Endpoints of new pixels
        nnew = len(new_velocity)
        nvlh = (new_velocity + np.roll(new_velocity, -1)) / 2.
        nvlh[nnew - 1] = new_velocity[nnew - 1] + \
                         (new_velocity[nnew - 1] - new_velocity[nnew - 2]) / 2.
        # Pad starting point
        bvl = np.zeros(nnew + 1)
        bvl[0] = new_velocity[0] - (new_velocity[1] - new_velocity[0]) / 2.
        bvl[1:] = nvlh

        # Evaluate
        newcum = fcum(bvl)

        # Rebinned flux
        new_fx = (np.roll(newcum, -1) - newcum)[:-1]

        # Normalize (preserve counts and flambda)
        new_dvl = bvl - np.roll(bvl, 1)
        new_fx = new_fx / new_dvl[1:]

        # Deal with the regions beyond the original data
        if fill_value == None:
            bd_vel = np.isnan(fcum_ref(bvl)[:-1])
            new_fx[bd_vel] = np.nan
        else:
            bd_vel = (fcum_ref(bvl) == fill_value)[:-1]
            new_fx[bd_vel] = fill_value


        # Return new spectrum
        self.velocity = new_velocity

        self.Tb = new_fx
        # Reset the mask
        self.mask = np.repeat(False, len(new_velocity))


    def copy(self):
        new = GBTspec(self.velocity.copy(),self.Tb.copy(),
                        mask=self.mask.copy(), meta = self.meta.copy())
        return new
