# pyND: Code for the ND CGM group

This is a bit of generic code used by the ND CGM group. It's focused on absorption line studies. It may not be of interest to anyone else, but it helps us keep our ducks in a row.

## Dependencies:

* numpy
* matplotlib
* astropy
* [linetools](https://github.com/linetools/linetools)

## Contributing

Please contact the owner before contributing.


## Structure
A rough outline of the basic structure of the module. This may be incomplete as code is added and not integrated into the list.

* **abs** – Absorption line spectroscopy routines.
  * get_fvals()
  * get_line()
  * logmean() – mean of log values
  * plotaxes() – plot 0, 1 for abs plots.
  * plotzero() – plot 0.
  * sensitivity() – calculate limiting EW, log N given SNR, instrument, transition.
  * sum_components() – Input: component columns; Returns: total column, err.


* **pyND.analysis** – General analysis support
  * solarabundance


* **pyND.hii** – HII region specific code
  * air_to_vac()
  * fit_lines() – Fit emission line spectrum, return intensities
  * fit_lines_sherpa() – As above, but w/Sherpa fitting code [this is preferred]
  * get_linelist() – Load HII region linelist for fitting.
  * integrate_line_flux() – Integration of lines given input centers.
  * prep_pyMCZ() – Export line intensities for use with pyMCZ code.


* **pyND.ism** – ISM research code.
  * ccm_filters() – Return A_lambda for input A_V, R_V.
  * ccm_extinct() – Return extinction at specified wavelength(s).
  * lsrvel() – Heliocentric to LSR conversion for input coordinates.
  * rotcurve() – Calculate MW rotation curve for input coordinates.
  * kinedist() – Calculate MW kinematic distance for input coords, vel. **[Not complete]**


* **pyND.ism.wham** – Access WHAM data
  * get_spectrum() – Load WHAM spectrum for input coords.
  * get_intensity() – Calculate WHAM integrated intens. for input coords.


* **pyND.lbt.mods** – Tools for MODS data
  * read_mods1d() – Load wave, flux, err from a single MODs spectrum.
  * join_mods1d() – Join MODS red, blue channels into a single spectrum.


* **pyND.uves** – Tools for UVES data
  * uves_log() – Create log file from raw UVES data.
  * simple_coadd() – Coadd only spectra with same cen waves.
  * full_coadd() – Coadd all spectra for a single object.
