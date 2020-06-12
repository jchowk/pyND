import numpy as np
from astropy.io import fits,ascii
import matplotlib.pyplot as plt

from astropy.table import Table
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.constants as c
import astropy.units as u



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
        slf.VELDEF = frm+'-'+vsys

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

        velocity = (nu0-nu)/nu0 * c.c.to('km/s')
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
        slf.VELDEF = b['VELDEF'][0]

        return slf

    def __init__(self, velocity, Tb):
        self.velocity = velocity
        self.Tb = Tb

        self.filename = None

        self.object = None
        self.RA = None
        self.DEC = None
        self.VELDEF = None


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

        # if self.VELDEF.find('LSR') > 0:
        #     xlbl_text = 'LSR Velocity [km/s]'
        # elif self.VELDEF.find('HELI') > 0:
        #     xlbl_text = 'Heliocentric Velocity [km/s]'
        # else:
        #     xlbl_text = '{0} Velocity [km/s]'.format(self.VELDEF)

        xlbl_text = '{0} Velocity [km/s]'.format(self.VELDEF)

        plt.xlabel(xlbl_text)
        plt.ylabel('$T_b$ [K]')
        plt.title(self.object)

        ylim=plt.ylim(ylim);
        xlim=plt.xlim(xlim);
