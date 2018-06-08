""" Contains the HIIspec class
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.units.quantity import Quantity
from astropy.table import Table, vstack, Column, MaskedColumn

from pyND.spec import plotzero, plotaxes
import pdb

CACHE = {'full_table': {}, 'data': {}}


class HIIspec(object):
    """Class for holding H II region spectrum and its properties.
    """

   # Init
    def __init__(self, llst_key, verbose=False, set_lines=True,
                 use_MODS_table=True, sort_by='wrest', redo_extra_columns=False):

        # Save
        self.list = llst_key

        if set_lines:
            # Set lines for use (from defined LineList)
            # This sets self._data
            self.set_lines(use_MODS_table=use_MODS_table, verbose=verbose,
                           use_cache=use_cache)
        else:
            self._data = None

    def load_data(self, use_MODS_table=True, tol=1e-3, use_cache=True):
        """Grab the data for the lines of interest
        Also load into CACHE

        Parameters
        ----------
        use_MODS_table : bool, optional
        tol : float, optional
          Tolerance for matching wavelength in AA
        use_cache : bool, optional

        Returns
        -------

        """

        # Define datasets: In order of Priority
        dataset = {
            'ism': [lilp.parse_morton03, lilp.parse_morton00, lilp.parse_verner96,
                    lilp.read_verner94, lilp.read_euv],  # Morton 2003, Morton 00, Verner 96, Verner 94
            'hi': [lilp.parse_morton03],
            # H2 (Abrigail), CO (JXP)
            'molecules': [lilp.read_H2, lilp.read_CO],
            'euv': [lilp.read_euv],  # EUV lines (by hand for now; soon to be Verner96)
            'galaxy': [lilp.read_forbidden, lilp.read_recomb, lilp.read_galabs],
        }

        sets = []
        flag_fval = False  # Update f-values?
        flag_wrest = False  # Update wavelengths?
        flag_gamma = True  # Update gamma values (recommended)








@property
def name(self):
    """ Return the transition names
    """
    return np.array(self._data['name'])


@property
def wrest(self):
    """ Return the rest wavelengths as a Quantity array
    """
    return Quantity(self._data['wrest'])


@property
def Z(self):
    """ Return the Z of the transitions
    """
    return self._data['Z']


@property
def ion(self):
    """ Return the ionization state of the transitions
    """
    return self._data['ion']


def load_data(self, use_ISM_table=True, tol=1e-3, use_cache=True):
    def fit_flux():

    def integrate_flux():
