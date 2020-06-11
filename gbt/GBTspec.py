import numpy as np
from astropy.io import fits,ascii
import matplotlib.pyplot as plt

import astropy.constants as c

class GBTspec(object):

    def __init__(self):
        self.velocity = None
        self.Tb = None

        self.filename = None
        self.velocity_frame = None

        self.object = None

        self.RA = None
        self.DEC = None
        self.VELDEF = None


    def from_ascii(self,input_filename,
                    format='no_header',
                    header_start=1,
                    data_start=3,
                    **kwargs):

        # Fill the information about the object:
        self.filename = input_filename

        # Read the file to grab the object name.
        #  [Assumes standard GBTIDL ASCII output.]
        rd = open (self.filename, "r")

        # Read list of lines
        out = rd.readlines()
        object_name = out[0].split()[2]
        rd.close()

        self.object = object_name


        a = ascii.read(self.filename,
            format=format,
            header_start=header_start,
            data_start=data_start,
            **kwargs)

        self.velocity = a[a.colnames[0]]
        self.Tb = a[a.colnames[1]]

    def from_GBTIDL(self,input_filename,object=None):
        # Fill the information about the object:
        self.filename = input_filename

        if object:
            self.object = object

        # Load the GBTIDL data:
        a = fits.open(self.filename)

        for j in np.arange(1,np.size(a)):
            xxx = a[j].data
            gd=(xxx['OBJECT'] == self.object)
            if gd.sum() > 0:
                b = xxx[gd]

        if ('b' in locals()) == False:
            print('No matching object found.')
            return

        # Define the spectrum:
        nu0=b['RESTFREQ']
        nu = ((np.arange(np.size(b['DATA']))+1)-b['CRPIX1'])*b['CDELT1'] + b['CRVAL1']

        self.velocity = (nu0-nu)/nu0 * c.c.to('km/s')
        self.Tb = b['DATA'][0]
        if b['TUNIT7'] == 'Ta*':
            self.Tb /= 0.88 # Main beam efficiency correction

        # Fill in some data/information from the GBTIDL format:
        self.RA = b['TRGTLONG']
        self.DEC = b['TRGTLAT']
        self.VELDEF = b['VELDEF'][0]


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
