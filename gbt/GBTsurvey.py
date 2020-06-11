import numpy as np
from astropy.io import fits

class GBTsurvey:

    def __init__(self,survey_file):
        self.survey_file = survey_file
        self.survey_FITStable = None


    def load_survey(self):
        a = fits.open(self.survey_file)
        self.survey_FITStable = a
