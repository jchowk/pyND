# pyND: Code for the ND CGM group

This repository contains a broad range of code used by the ND CGM group. It iss focused on absorption line studies. It may not be of interest to anyone else, but it helps us keep our ducks in a row.

## Installation:

See the [Installation](Installation.md) documentation.


## Dependencies:

* `numpy`
* `scipy`
* `matplotlib`
* `astropy`
* [`linetools`](https://github.com/linetools/linetools)

## Contributing

Please contact the owner before contributing.

## Structure
A rough outline of the basic structure of the module. This may be incomplete as code is added and not integrated into the list.

* **pyND.absorption** – Absorption line spectroscopy routines.
  * `get_fvals()`
  * `get_line()`
  * `logmean()` – mean of log values
  * `sensitivity()` – calculate limiting EW, log N given SNR, instrument, transition.
  * `sum_components()` – Input: component columns; Returns: total column, err.
  * `correct_saturation()` – Calculate the Savage & Sembach (1991) saturation correction. 


* **pyND.analysis** – General analysis support
  * `solarabundance()` - Retrieve default solar abundances for a list of elements. [note the availability of abundances by mass with Xinit, Yinit, Zinit, Xphoto, Yphoto, Zphoto]
  * `coveringfactor()` - Calculate a covering factor assuming Bayesian statistics using the inverse-beta function. 


* **pyND.gal.halos** – Galaxy DM halo relationships
    * `smhm()` – Calculate M_halo = M200 for a given stellar mass. Uses default relationship [currently `smhm_rodriguez()`.]
    * `smhm_behroozi()` - Calculate M_halo = M200 given stellar mass using Behroozi relationships.
    * `smhm_rodriguez()` – Calculate M_halo = M200 given stellar mass using Rodriguez-Puebla+ (2017) relationships.
    * `smhm_shan()` - Calculate M_halo = M200 given stellar mass using Shan+ (2017) relationships.
    * `smhm_tinker()` - Calulate M_halo for a given stellar mass using Tinker+ (2017) relationships. *Defaults to `smhm_shan` for log M < 11.04.*
    * `virial_radius()` - Calculate virial radius given M_halo (several assumptions available).
    * `calc_r200()` - Convenience method for calculating R200c from M200c [uses `virial_radius`].
    * `log_schechter_function()` - Calculates a 2-component Schechter function given mass function parameters. 

* **pyND.hii** – HII region specific code
  * `fit_lines()` – Fit emission line spectrum, return intensities
  * `fit_lines_sherpa()` – As above, but w/Sherpa fitting code [this is preferred]
  * `get_linelist()` – Load HII region linelist for fitting.
  * `integrate_line_flux()` – Integration of lines given input centers.
  * `prep_pyMCZ()` – Export line intensities for use with pyMCZ code.


* **pyND.ism** – ISM code.
  * `ccm_filters()` – Return A_lambda for input A_V, R_V.
  * `ccm_extinct()` – Return extinction at specified wavelength(s).
  * `lsrvel()` – Heliocentric to LSR conversion for input coordinates.
  * `rotcurve()` – Calculate MW rotation curve for input coordinates.
  * `kinedist()` – Calculate MW kinematic distance for input coords, vel. **[Not complete]**


* **pyND.ism.wham** – Access WHAM data
  * `get_spectrum()` – Load WHAM spectrum for input coords.
  * `get_intensity()` – Calculate WHAM integrated intens. for input coords.


* **pyND.mods** – Tools for MODS data
  * `read_mods1d()` – Load wave, flux, err from a single MODs spectrum.
  * `join_mods1d()` – Join MODS red, blue channels into a single spectrum.

* **pyND.plotting** – Tools for plotting.
    * `plotaxes()`         – Plot lines at 0, 1. Useful for absorption plots.
    * `plotzero()`         – Plot line 0.
    * `cosmictimeaxis()`   – Plot cosmic time since BigBang on top x-axis.
    * `lookbacktimeaxis()` – Plot lookback time on top x-axis.
    * `skyplot()`          – Make sky plots with optional z-scale colors. [Notebook in `docs/Examples/skyplot_usage.ipynb`]

* **pyND.spec** – Tools for spectroscopy.
    * `rebin()`      – Rebin a spectrum by an integer factor. [~IDL rebin]
    * `congrid()`    – Rebin a spectrum by an arbitrary factor. [~IDL congrid]
    * `resample()`   - Resample a spectrum to an _arbitrary_ vel/wave vector, flux conserving.
    * `vac_to_air()` – Convert vacuum to air wavelengths.
    * `air_to_vac()` – Convert air to vacuum wavelengths.

* **pyND.uves** – Tools for UVES data
  * `uves_log()` – Create log file from raw UVES data.
  * `simple_coadd()` – Coadd only spectra with same cen waves.
  * `full_coadd()` – Coadd all spectra for a single object.
