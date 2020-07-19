import numpy as np
from astropy.io import fits
from astropy.table import Table
import pickle

from .GBTspec import *

class GBTsurvey:

    @classmethod
    def from_GBTIDL(cls,survey_file):
    # Load an entire GBT survey.

        # Load the GBTIDL data:
        a = fits.open(survey_file)

        survey_data = []
        survey_meta = {'object':[],
                        'idx':[],
                        'l':[],'b':[],
                        'RA':[], 'DEC':[]}
        # survey_sightline = []
        # survey_RA = []
        # survey_Dec = []
        # survey_l = []
        # survey_b = []

        counter=0
        for j in np.arange(1,len(a)):
            xxx = a[j].data
            for k in np.arange(len(xxx)):
                b = xxx[k]
                # Append the actual object to the list
                bb=GBTsurvey._load_object(b)
                survey_data.append(bb)

                # Append metadata
                survey_meta['object'].append(bb.meta['object'])
                survey_meta['RA'].append(bb.meta['RA'])
                survey_meta['DEC'].append(bb.meta['DEC'])
                survey_meta['l'].append(bb.meta['l'])
                survey_meta['b'].append(bb.meta['b'])
                survey_meta['idx'].append(counter)

                counter+=1

        # Initiate the class
        slf = cls(survey_data)

        # META DATA
        # Fill the information about the survey:
        slf.GBTFITS_file = os.path.abspath(survey_file)
        slf.index = Table(survey_meta)

        return slf

    @classmethod
    def from_pickle(cls,pickle_file):
    # Load an entire GBT survey.

        # Load the GBTIDL data:
        slf = pickle.load(open(survey_file,'rb'))
        slf.data_file=pickle_file

        return slf

    def __init__(self,survey_data):
        self.GBTFITS_file = None
        self.data_file = None
        self.index = None
        self.survey = survey_data


    def _load_object(b):
    # Populate a GBTspec object

        # Define the spectrum:
        nu0=b['RESTFREQ']
        nu = ((np.arange(np.size(b['DATA']))+1)-b['CRPIX1'])*b['CDELT1'] + b['CRVAL1']

        velocity = (nu0-nu)/nu0 * c.c.to('km/s').value
        Tb = b['DATA']
        if b['TUNIT7'] == 'Ta*':
            Tb /= 0.88 # Main beam efficiency correction

        # Initiate object:
        gbt_object = GBTspec(velocity, Tb)

        # Fill in some data/information from the GBTIDL format:
        gbt_object.meta['object'] = b['OBJECT']
        gbt_object.meta['RA'] = b['TRGTLONG']
        gbt_object.meta['DEC'] = b['TRGTLAT']

        # Define the Galactic coordinates
        gbt_object._fill_Galactic_coords()

        # Details of the GBT data
        gbt_object.meta['veldef'] = b['VELDEF']
        gbt_object.meta['restfreq'] = b['RESTFREQ']

        gbt_object.meta['Tsys'] = b['TSYS']

        return gbt_object

    def save2pickle(self,pickle_file):
        """Save a GBTsurvey object to a pickle file."""

        # Load the GBTIDL data:
        self.data_file=pickle_file
        pickle.dump(self,open(pickle_file,'wb'),protocol=2)



    def sort(self,column_name,reverse=False):
        """Sort the GBT data on the basis of a column in the index."""
        yyy=self.index[column_name].data
        ss=np.argsort(yyy)

        if reverse is True:
            ss = ss[::-1]

        # Sort the index table
        self.index = self.index[ss]
        self.index['idx'] = np.arange(len(yyy))

        # Sort the datasets
        self.survey = [self.survey[i] for i in ss]


    def copy(self):
        new_survey = []
        for i in np.arange(len(self.survey)):
            old_GBTspec = self.survey[i]
            new_GBTspec = old_GBTspec.copy()
            new_survey.append(new_GBTspec)


        new = GBTsurvey(new_survey.copy())
        new.index = self.index.copy()
        new.GBTFITS_file = self.GBTFITS_file
        new.data_file = self.data_file

        return new

    def resample_survey(self, new_velocity, fill_value=None):
        survey_size = len(self.index)
        for j in np.arange(survey_size):
            self.survey[j].resample(new_velocity,fill_value)
