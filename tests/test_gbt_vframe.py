"""
Tests that GBTIDL SDFITS readers return velocities in the VELDEF frame (e.g.
LSR) rather than the topocentric FREQ-OBS frame, i.e. that VFRAME is applied,
and that the Ta* -> T_B correction does not modify the FITS data in place.
"""
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from pyND.gbt.GBTspec import GBTspec
from pyND.gbt.GBTsurvey import GBTsurvey

GBT_DATA = Path(__file__).parent.parent / 'pyND' / 'gbt' / 'docs' / 'data'
GBTIDL_FILE = GBT_DATA / 'GBTdata.fits'
ASCII_FILE = GBT_DATA / 'HS0033+4300_GBT.dat'
CKMS = 299792.458

needs_data = pytest.mark.skipif(not (GBTIDL_FILE.exists() and ASCII_FILE.exists()),
                                reason='bundled GBT test data not present')


def _ascii_radio_lsr():
    """GBTIDL ASCII export of HS0033+4300 (OPTICAL, LSR) converted to radio velocities."""
    a = np.loadtxt(ASCII_FILE, skiprows=3)
    v_opt, tb = a[:, 0], a[:, 1]
    return CKMS * v_opt / (CKMS + v_opt), tb


def _centroid(vel, tb, vmin=-40., vmax=20.):
    """T_B-weighted centroid of the Galactic HI emission."""
    s = (vel > vmin) & (vel < vmax) & np.isfinite(tb)
    return np.sum(vel[s] * tb[s]) / np.sum(tb[s])


def _hs0033_row():
    for ext in fits.open(GBTIDL_FILE)[1:]:
        for row in ext.data:
            if row['OBJECT'].strip() == 'HS0033+4300':
                return row
    raise LookupError('HS0033+4300 not in test file')


@needs_data
def test_load_object_matches_gbtidl_lsr():
    """GBTsurvey._load_object velocities agree with GBTIDL's LSR export."""
    row = _hs0033_row()
    assert abs(row['VFRAME']) > 3.e3        # test is only meaningful if VFRAME is sizable
    spec = GBTsurvey._load_object(row)
    v_ref, tb_ref = _ascii_radio_lsr()
    assert _centroid(spec.velocity, spec.Tb) == pytest.approx(_centroid(v_ref, tb_ref), abs=0.5)


@needs_data
def test_from_GBTIDL_matches_gbtidl_lsr():
    spec = GBTspec.from_GBTIDL(str(GBTIDL_FILE), object_name='HS0033+4300')
    v_ref, tb_ref = _ascii_radio_lsr()
    assert _centroid(spec.velocity, spec.Tb) == pytest.approx(_centroid(v_ref, tb_ref), abs=0.5)


@needs_data
def test_from_GBTIDLindex_matches_from_GBTIDL():
    tbl = GBTspec.index_GBTIDL(str(GBTIDL_FILE), silent=True, return_list=True)
    names = [str(n).strip() for n in tbl[tbl.colnames[1]]]
    spec_i = GBTspec.from_GBTIDLindex(str(GBTIDL_FILE), object_indx=names.index('HS0033+4300'))
    spec_n = GBTspec.from_GBTIDL(str(GBTIDL_FILE), object_name='HS0033+4300')
    np.testing.assert_allclose(spec_i.velocity, spec_n.velocity)
    np.testing.assert_allclose(spec_i.Tb, spec_n.Tb)


@needs_data
def test_load_object_does_not_modify_fits_data():
    """Calling the loader twice on the same row must give the same T_B."""
    row = _hs0033_row()
    raw = np.array(row['DATA'], dtype=float).copy()
    t1 = GBTsurvey._load_object(row).Tb.copy()
    t2 = GBTsurvey._load_object(row).Tb.copy()
    np.testing.assert_array_equal(np.array(row['DATA'], dtype=float), raw)
    np.testing.assert_allclose(t1, t2)
    np.testing.assert_allclose(t1, raw / 0.88)


def test_frame_velocity_exact_shift():
    """Synthetic row: a line at the rest frequency observed from a frame moving
    at VFRAME appears at radio velocity ~ -VFRAME in the topocentric frame and
    must come back at 0 km/s in the VELDEF frame."""
    from pyND.gbt.GBTspec import _sdfits_velocity

    nu0 = 1.420405752e9
    vframe = -20.e3                                    # m/s
    beta = vframe / 1.e3 / CKMS
    nu_obs_line = nu0 * np.sqrt((1. - beta) / (1. + beta))  # observed freq of a v_frame=0 line
    row = {'RESTFREQ': nu0, 'CRPIX1': 101., 'CDELT1': -5722., 'CRVAL1': nu_obs_line,
           'VFRAME': vframe, 'CTYPE1': 'FREQ-OBS', 'DATA': np.zeros(201)}
    vel = _sdfits_velocity(row)
    assert vel[100] == pytest.approx(0., abs=1.e-6)
    # Without a topocentric axis no shift is applied
    row['CTYPE1'] = 'FREQ-LSR'
    row['CRVAL1'] = nu0
    assert _sdfits_velocity(row)[100] == pytest.approx(0., abs=1.e-6)
