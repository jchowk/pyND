# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

pyND is the Notre Dame CGM group's shared Python library for absorption line spectroscopy and related astrophysics analysis. It is managed with Poetry (`pyproject.toml`) and installed editable via pip.

## Installation

```bash
# Editable install (workaround for Poetry + editable mode):
touch setup.cfg
pip install -e .
rm setup.cfg

# Or directly from GitHub:
pip install git+https://github.com/jchowk/pyND.git
```

## Running Tests

There is no formal test suite. Development is validated through Jupyter notebooks in `docs/Examples/`.

## Module Architecture

Each subpackage lives under `pyND/` and exposes its public API via `__init__.py` using `from .module import *` or named imports. The top-level `pyND/__init__.py` is minimal — submodules are imported explicitly by users (e.g., `from pyND.absorption import get_line`).

| Subpackage | Purpose |
|---|---|
| `absorption/` | Absorption line EW/column density tools; saturation corrections; sensitivity calculations. Uses `linetools` for line data. |
| `analysis/` | Covering factor (Bayesian), asymmetric error distributions, solar abundances (`solarabundance.py` + `data/`). |
| `gal/` | Stellar-to-halo mass relations (Behroozi, Rodriguez-Puebla, Shan, Tinker), virial radius calculations, Schechter function. |
| `spec/` | Wavelength-space utilities: rebinning (`rebin`, `congrid`, `resample`), vacuum/air conversions. |
| `ism/` | Extinction (CCM), helio→LSR velocity, MW rotation curve. `wham.py` fetches WHAM HII-region survey spectra. |
| `hii/` | Emission line fitting (plain and Sherpa), pyMCZ export. |
| `plotting/` | Axis helpers (`plotaxes`, `plotzero`), cosmic-time/lookback-time secondary axes, sky plots. |
| `gbt/` | GBT survey and spectrum classes (`GBTsurvey`, `GBTspec`), Doppler/velocity utilities. |
| `hst/` | HST/COS data tools. |
| `mods/` | LBT/MODS 1D spectrum I/O. |
| `uves/` | VLT/UVES log creation and coaddition. |
| `stis/` | HST/STIS utilities. |
| `cloudy/` | Cloudy photoionization model fitting (`cloudy_fits.py`). |

## Key Dependencies

- `astropy` ≥ 5.0 — units, coordinates, FITS, constants throughout
- `linetools` ≥ 0.3 — atomic line data used heavily in `absorption/`
- `numpy`, `scipy` — numerical core
- `matplotlib` — all plotting

## Code Conventions

- Functions import their own dependencies internally (local imports at function scope are common).
- NumPy-style docstrings throughout.
- `astropy.units` quantities are used where physical units matter; do not silently strip units.
- Submodule `__init__.py` files control what is re-exported — keep them updated when adding public functions.
- No linting or formatting config is present; match the existing style of the file being edited.
