import numpy as np
from astropy.io import fits,ascii
import matplotlib.pyplot as plt

class GBTspec(object):

    def __init__(self):
        self.velocity = None
        self.Tb = None

        self.filename = None
        self.velocity_frame = None

        self.object = None

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

    def plotspectrum(self,**kwargs):
        try:
            kwargs['color']
        except:
            kwargs['color']='green'

        plt.figure(figsize=(10,5))

        plt.plot(self.velocity,self.Tb,drawstyle='steps-mid',
                    **kwargs)
        plt.axhline(0.,linestyle='--',linewidth=1,color='k',zorder=0)
        plt.xlabel('Velocity [km/s]')
        plt.ylabel('$T_b$ [K]')
        plt.title(self.object)

        ylim=plt.ylim(kwargs.get("ylim"));
        xlim=plt.xlim(kwargs.get("xlim"));
