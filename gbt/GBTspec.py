import numpy as np
from astropy.io import fits,ascii
import matplotlib.pyplot as plt

from astropy.table import Table
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.constants as c
import astropy.units as u

# TODO: Add copy method.
# TODO: Change OPTICAL to RADIO LSR
# TODO: Heliocentric to RADIO
# TODO: Rebin / congrid spectra
# TODO: Calculate HI column over velocity range
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

        # Temperature system
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
        slf.filename = input_filename

        # Read the file as a whole to grab information.
        #  [Assumes standard GBTIDL ASCII output.]
        rd = open (slf.filename, "r")

        # Read list of lines
        out = rd.readlines()
        rd.close()

        # Parse the object name.
        object_name = out[0].split()[2]

        # Set the velocity frame.
        frm = out[1].split()[0]
        vsys = out[2].split()[0].split('-')[1]
        slf.veldef = frm+'-'+vsys

        # Apply the object name
        slf.object = object_name

        simbad_result = Simbad.query_object(slf.object)
        if isinstance(simbad_result,Table):
            coords = SkyCoord(simbad_result['RA'],simbad_result['DEC'],
                  unit=(u.hourangle, u.deg))
            slf.RA = coords.ra.deg[0]
            slf.DEC = coords.dec.deg[0]

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
        slf.filename = input_filename
        slf.object = object_name

        # Fill in some data/information from the GBTIDL format:
        slf.RA = b['TRGTLONG'][0]
        slf.DEC = b['TRGTLAT'][0]
        slf.veldef = b['VELDEF'][0]

        slf.restfreq = b['RESTFREQ'][0]

        return slf

    @classmethod
    def from_GBTindex(cls,input_filename,object_name=None):

        # Load the GBTIDL data:
        a = fits.open(input_filename)
        tbl=GBTspec.index_GBTIDL(input_filename,
                        return_list=True)

        valid = False
        while valid == False:
            try:
                indx = np.int(input("Index of object to load: "))
            except:
                # Not a numerical entry.
                print("Enter numerical value from list of indeces.")
            if (indx <= np.max(tbl[tbl.colnames[0]])) & \
                    (indx >= 0):
                valid = True
                break
            else:
                print("Enter an index within the range listed.")



        # What array indeces
        array_indx = tbl['ARRAY_INDECES'][indx]

        # from IPython import embed ; embed()

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
        # Fill the information about the object:
        slf.filename = input_filename
        slf.object = object_name

        # Fill in some data/information from the GBTIDL format:
        slf.RA = b['TRGTLONG']
        slf.DEC = b['TRGTLAT']
        slf.veldef = b['VELDEF']

        slf.restfreq = b['RESTFREQ']

        return slf

    def __init__(self, velocity, Tb):
        self.velocity = velocity
        self.Tb = Tb

        self.filename = None

        # Information from the GBTIDL structure
        self.object = None
        self.RA = None
        self.DEC = None
        self.veldef = None
        self.restfreq = None


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

        # if self.veldef.find('LSR') > 0:
        #     xlbl_text = 'LSR Velocity [km/s]'
        # elif self.veldef.find('HELI') > 0:
        #     xlbl_text = 'Heliocentric Velocity [km/s]'
        # else:
        #     xlbl_text = '{0} Velocity [km/s]'.format(self.veldef)

        xlbl_text = '{0} Velocity [km/s]'.format(self.veldef)

        plt.xlabel(xlbl_text)
        plt.ylabel('$T_b$ [K]')
        plt.title(self.object)

        ylim=plt.ylim(ylim);
        xlim=plt.xlim(xlim);

    def index_GBTIDL(input_filename, silent=False,
                        return_list=False):
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

                if ~silent:
                    print("{0}: {1:20s}   \t{2}".format(i,
                        object_names[i],array_indeces[i]))

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

        if not self.restfreq:
            self.restfreq = nu0

        if self.veldef == 'OPTICAL-LSR':
            # Change OPTICAL-->RADIO definition of LSR

            # First calculate "TRUE" velocity
            result = self.velocity*1000. / light_speed
            result = (2.e0 * result + result * result) / (2.e0 + 2.e0 * result + result * result)

            # Calculate RADIO definition
            result = 1.e0 - np.sqrt((1.e0 - result) / (1.e0 + result))

            # Update:
            self.velocity = result*light_speed/1000.
            self.veldef = 'RADI-LSR'

        elif self.veldef == 'RADI-LSR':
            # Change RADIO-->OPTICAL definition of LSR

            # First calculate "TRUE" velocity
            result = self.velocity*1000. / light_speed
            result = (2.e0 * result - result * result) / (2.e0 - 2.e0 * result + result * result)

            # Calculate OPTICAL definition
            result = np.sqrt((1.e0 + result) / (1.e0 - result)) - 1.e0

            # Update veldef:
            self.velocity = result*light_speed/1000.
            self.veldef = 'OPTICAL-LSR'
        else:
            print('Inappropriate velocity definition.')
