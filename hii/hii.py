""" Extract HII region line fluxes
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import numpy as np
from astropy.io import ascii
import astropy.constants as c

from astropy.table import Table, Column

from astropy.modeling.parameters import Parameter
from astropy.modeling.functional_models import *

import pdb

class GaussianEmission(BaseGaussian1D):
    """
    One dimensional Gaussian model.

    Parameters
    ----------
    amplitude : float
        Amplitude of the Gaussian.
    redshift : float
        Mean redshift of the Gaussian.
    stddev : float
        Standard deviation of the Gaussian.
    wave0 : float
        Central wavelength of the line.
    """

    amplitude = Parameter(default=1)
    redshift = Parameter(default=0)
    stddev = Parameter(default=1)
    wave0 = Parameter(default = 5008.240)

    @staticmethod
    def evaluate(x, amplitude, redshift, stddev, wave0):
        """
        Gaussian1D model function.
        """
        return amplitude * np.exp(- 0.5 * (x - (1+redshift)*wave0) ** 2 / stddev ** 2)

    @staticmethod
    def fit_deriv(x, amplitude, redshift, stddev, wave0):
        """
        Gaussian1D model function derivatives.
        """

        d_amplitude = np.exp(-0.5 / stddev ** 2 * (x - (1+redshift)*wave0) ** 2)
        d_mean = amplitude * d_amplitude * (x - (1+redshift)*wave0) / stddev ** 2
        d_stddev = amplitude * d_amplitude * (x - (1+redshift)*wave0) ** 2 / stddev ** 3
        d_wave0 = 0
        return [d_amplitude, d_mean, d_stddev, d_wave0]

def get_linelist(use_mods_table=True):
    """Read emission line list data."""
    if use_mods_table:
        data_dir = os.path.join(os.path.dirname(__file__),'data/')
        data_file = 'mods_linelist.dat'

    # Read in the emission line table
    line_data=ascii.read(data_dir+data_file,comment='#',format='no_header',
                 names=['indx','ion','lambda','action','line','A_i',
                        'v_g','sig_g','mode'])
    line_data.keep_columns(['indx','ion','lambda','action','line','A_i','mode'])

    # Create name for identification
    transition_name = []
    for j in np.arange(np.size(line_data)):
        lambda_name = np.str(np.int(np.floor(line_data['lambda'][j])))
        ion_name = line_data['ion'][j]
        transition_name.append(ion_name+lambda_name)

    # Make a column of this, add it to the line_data
    transition_column = Column(transition_name,name='name')
    line_data.add_column(transition_column,index=0)

    # Extract those lines to be part of the measurement group
    use_for_fit = np.where(line_data['action'] == 'f')

    return line_data[use_for_fit]

def fit_lines_sherpa(spec_file, z_init=0., do_plot=True, monte_carlo=False):
    """Fit an HII region spectrum.    """

    from astropy.modeling.fitting import SherpaFitter
    from astropy.modeling import models, fitting
    import matplotlib.pyplot as plt
    from linetools.spectra.xspectrum1d import XSpectrum1D

    def tie_sigma_4862(model):
        # Tie the dispersions to Hbeta 4862
        return model['33'].stddev

    def tie_redshift_4862(model):
        # Tie the redshift to Hbeta 4862
        return model['33'].redshift

    def tie_sigma_5008(model):
        # Tie the dispersions to O III 5008
        return model['35'].stddev

    def tie_redshift_5008(model):
        # Tie the redshift to O III 5008
        return model['35'].redshift

    def tie_ampl_5008(model):
        #Tie 4959 flux to that of 5008
        return model['35'].amplitude*0.350

    def tie_sigma_6585(model):
        # Tie the dispersions to N II 6585
        return model['45'].stddev

    def tie_redshift_6585(model):
        # Tie the redshift to N II 6585
        return model['45'].redshift

    def tie_ampl_6585(model):
        #Tie 6549 flux to that of 6585
        return model['45'].amplitude*0.340

    # Redshift scale:
    scale_factor = (1.+z_init)

    # Read in the spectrum. **ASSUME VACUUM WAVELENGTHS?**
    mods_spec = XSpectrum1D.from_file(spec_file)

    # Set up a convenient wavelength, flux, error arrays
    wave = mods_spec.wavelength.value
    flux = mods_spec.flux.value
    err = mods_spec.sig.value

    # Grab the data for the lines to be fit.
    line_data = get_linelist()
    # Exclude lines outside of the wavelength coverage.
    keep_lines = np.where((line_data['lambda'] >= mods_spec.wvmin.value/scale_factor) &
                          (line_data['lambda'] <= mods_spec.wvmax.value/scale_factor))
    line_data = line_data[keep_lines]


    ##### MODEL DEFINITIONS
    # Define a joint model as the sums of Gaussians for each line

    # Define initial parameters
    amplitude_init = 0.1*np.max(mods_spec.flux)
    stddev_init = 2.

    # Set parameters constraints
    amplitude_bounds = (0.,1.e-12)  # Amplitudes must be >=0.
    stddev_bounds = (1.0,3.00)       # stddev between 0.5 and 5 Ang

    # Redshift constraints
    velocity_range = 500.   # Velocity range to explore about center
    # Calculate the redshift delta
    z_bounds_scale = (velocity_range/c.c.to('km/s').value)*scale_factor
    z_bounds = (z_init-z_bounds_scale, z_init+z_bounds_scale)

    #mean_bounds_scale = (velocity_range / c.c.to('km/s').value) * np.array([-1., 1.]) + 1.
    #mean_bounds = []

    #  Initial Gaussian:
    j=0
    wave0 = line_data['lambda'][j]
    line_center = wave0*scale_factor
    model_name = np.str(line_data['indx'][j])
    # Here we use a class defined above to allow fixing redshifts
    joint_model = GaussianEmission(amplitude=amplitude_init,redshift=z_init,
                                   stddev=stddev_init,wave0=wave0,
                                   name=model_name)
    joint_model.wave0.fixed = True

    # Traditional astropy 1D Gaussian fit:
    # joint_model = Gaussian1D(amplitude=amplitude_init,mean=line_center,
    #                                stddev=stddev_init,
    #                                name=model_name)
    # Set constraints on how much the central value can vary
    # mean_bounds.append([line_center*mean_bounds_scale[0],line_center*mean_bounds_scale[1]])

    #  Loop through the remaining lines:
    for j in np.arange(1,np.size(line_data)):
        wave0 = line_data['lambda'][j]
        line_center = wave0 * scale_factor
        model_name = np.str(line_data['indx'][j])

        joint_model += GaussianEmission(amplitude=amplitude_init,redshift=z_init,
                                        stddev=stddev_init, wave0=wave0,
                                        name=model_name)

    # Now we have to loop through the same models, applying the bounds:
    for k in np.arange(0, np.size(line_data)):
        joint_model[k].bounds['amplitude'] = amplitude_bounds
        joint_model[k].bounds['redshift'] = z_bounds
        joint_model[k].bounds['stddev'] = stddev_bounds
        joint_model[k].wave0.fixed = True

    # TODO Get tied parameters to work.
    #  Loop through lines to define tied parameters together
    #for k in np.arange(0, np.size(line_data)):
        # if line_data['mode'][k] == 't33':
        #     joint_model[k].stddev.tied = tie_sigma_4862
        #     joint_model[k].redshift.tied = tie_redshift_4862
        # elif line_data['mode'][k] == 't35':
        #     joint_model[k].stddev.tied = tie_sigma_5008
        #     joint_model[k].redshift.tied = tie_redshift_5008
        # elif line_data['mode'][k] == 't45':
        #     joint_model[k].stddev.tied = tie_sigma_6585
        #     joint_model[k].redshift.tied = tie_redshift_6585
        #
        # # Tie fluxes for doublets:
        #if line_data['line'][k] == 'd35':
        #    joint_model[k].amplitude.tied = tie_ampl_5008
        #    joint_model[k].amplitude.value = joint_model[k].amplitude.tied(joint_model)
        # if line_data['line'][k] == 'd45':
        #     joint_model[k].amplitude.tied = tie_ampl_6585
        #     joint_model[k].amplitude.value = joint_model[k].amplitude.tied(joint_model)

    # Constrain the wavelengths over which the fits are calculated (don't need the continuum)
    # Initialize the boolean indeces:
    # fit_me = (wave < 0.)
    # fit_delta = 25.  # A guess...
    # for j in np.arange(0,np.size(line_data)):
    #     fit_me += (np.abs(wave/scale_factor - line_data['lambda'][j]) <= fit_delta)

    ##### FITTING
    # Sherpa model fitting from SABA package
    sfit = SherpaFitter(statistic='chi2', optimizer='levmar', estmethod='confidence')
    sfit_lm = SherpaFitter(statistic='chi2', optimizer='neldermead', estmethod='confidence')
    sfit_mc = SherpaFitter(statistic='chi2', optimizer='moncar', estmethod='confidence')
    # Do the fit
    sfitted_model = sfit(joint_model, wave, flux, err = err)
    # Refine with different optimizer
    sfitted_model = sfit_lm(sfitted_model.copy(), wave, flux, err = err)

    if monte_carlo:
        # If requested, do a second fit with the very slow Monte Carlo approach
        sfitted_model = sfit_mc(sfitted_model.copy(), wave, flux, err = err)

    # Create the fitted flux array
    sfitted_flux = sfitted_model(wave)

    # TODO Get error estimates from Sherpa
    # Work out the errors...
    #sfit.est_config['maxiters']=200
    #sfitted_err = sfit.est_errors(sigma=3)

    # Plot the results
    if do_plot:
        plt.clf()
        plt.plot(wave,flux,drawstyle='steps-mid',linewidth=2)
        plt.plot(wave,sfitted_flux,color='orange',linewidth=2)

    ##### Create integrated fluxes and errors
    # The integration range is over +/-stddev * int_delta_factor
    int_delta_factor = _define_integration_delta()

    for j in np.arange(np.size(line_data)):
        # Calculate integrated fluxes, errors;
        #    deal with blended O II 3727/3729 doublet

        mean_lambda = line_data[j]['lambda']*(1.+sfitted_model[j].redshift)

        if line_data[j]['name'] == '[OII]3727':
            # Calculate the integrated fluxes and errors
            iflux, ierr = integrate_line_flux(wave, flux, err,
                                              mean_lambda,
                                              sfitted_model[j].stddev * int_delta_factor,
                                              line3727=True)
        elif line_data[j]['name'] == '[OII]3729':
            # pdb.set_trace()
            # For 3729, use the flux derived for 3726
            iflux = output_table['fit_flux'][j - 1]
            # For 3729, use its own error. This is appropriate for the fitted errors of both liness
            crap, ierr = integrate_line_flux(wave, flux, err, mean_lambda,
                                             sfitted_model[j].stddev * int_delta_factor)
        else:
            # Calculate the integrated fluxes and errors
            iflux, ierr = integrate_line_flux(wave, flux, err, mean_lambda,
                                              sfitted_model[j].stddev * int_delta_factor)

        redshift_out = sfitted_model[j].redshift
        sfitted_flux_out = np.sqrt(2.*np.pi)*sfitted_model[j].amplitude*sfitted_model[j].stddev

        if j == 0:
            # Define and construct the initial table to hold the results
            output_col_names, output_format, output_dtype = _define_output_table()
            output_data = [[line_data[j]['name']],
                           [line_data[j]['ion']],
                           [line_data[j]['lambda']],
                           [line_data[j]['indx']],
                           [line_data[j]['mode']],
                           [mean_lambda],
                           [sfitted_model[j].amplitude.value],
                           [sfitted_model[j].stddev.value],
                           [redshift_out], [sfitted_flux_out],
                           [iflux], [ierr]]
            output_table = Table(output_data, names=output_col_names,dtype=output_dtype)
        else:
            output_table.add_row([line_data[j]['name'],line_data[j]['ion'],
                              line_data[j]['lambda'], line_data[j]['indx'],line_data[j]['mode'],
                              mean_lambda,sfitted_model[j].amplitude.value,
                              sfitted_model[j].stddev.value,
                              redshift_out, sfitted_flux_out, iflux, ierr])

    # Set the output format of the results table:
    colnames = output_table.colnames
    for j in np.arange(np.size(colnames)):
        output_table[colnames[j]].format = output_format[j]

    # Set up the spectral table:
    spec_table = Table([wave,flux,err,sfitted_flux],
                       names=['wave','flux','err','spec_fit'])

    # Write summary FITS files
    file_base = spec_file.rsplit('.')[0]
    table_file = file_base+'.HIIFitTable.fits'
    fit_file = file_base+'.HIIFitSpec.fits'

    output_table.write(table_file,overwrite=True)
    spec_table.write(fit_file,overwrite=True)

    return output_table


def fit_lines(spec_file, z_init=0., do_plot=True):
    from astropy.modeling import models, fitting
    from linetools.spectra.xspectrum1d import XSpectrum1D
    import matplotlib.pyplot as plt

    def tie_sigma_3729(model):
        # Tie the dispersions to O II 3729
        return model['17'].stddev

    def tie_sigma_4862(model):
        # Tie the dispersions to Hbeta 4862
        return model['33'].stddev

    def tie_sigma_5008(model):
        # Tie the dispersions to O III 5008
        return model['35'].stddev

    def tie_ampl_5008(model):
        #Tie 4959 flux to that of 5008
        return model['35'].amplitude*0.350

    def tie_mean_3727_5008(model):
        # Tie the redshift to O III 5008
        return model['35'].mean * 3727.092 / 5008.240

    def tie_mean_3729_5008(model):
        # Tie the redshift to O III 5008
        return model['35'].mean * 3729.875 / 5008.240

    def tie_mean_4364_5008(model):
        # Tie the mean of 4363 auroral line to O III 5008
        return model['35'].mean *  4364.435 / 5008.240
    def tie_mean_4960_5008(model):
        # Tie the mean of 4363 auroral line to O III 5008
        return model['35'].mean *  4960.295 / 5008.240

    def tie_sigma_6585(model):
        # Tie the dispersions to N II 6585
        return model['45'].stddev

    def tie_ampl_6585(model):
        #Tie 6549 flux to that of 6585
        return model['45'].amplitude*0.340



    # Redshift scale:
    scale_factor = (1.+z_init)

    # Read in the spectrum. **ASSUME VACUUM WAVELENGTHS?**
    mods_spec = XSpectrum1D.from_file(spec_file)

    # Set up a convenient wavelength, flux, error arrays
    wave = mods_spec.wavelength.value
    flux = mods_spec.flux.value
    err = mods_spec.sig.value

    # Grab the data for the lines to be fit.
    line_data = get_linelist()
    # Exclude lines outside of the wavelength coverage.
    keep_lines = np.where((line_data['lambda'] >= np.min(wave)/scale_factor) &
                          (line_data['lambda'] <= np.max(wave)/scale_factor))
    line_data = line_data[keep_lines]

    ##########
    # MODEL DEFINITION
    # Define a joint model as the sums of Gaussians for each line

    # Initial parameters
    amplitude_init = 0.1*np.max(mods_spec.flux)
    stddev_init = 2.

    # Constraint parameters:
    amplitude_bounds = (0.,None)
    stddev_bounds = (1.0,3.0)
    velocity_range = 500.   # Velocity range to explore about center
    mean_bounds_scale = (velocity_range / c.c.to('km/s').value) * np.array([-1., 1.]) + 1.
    mean_bounds = []

    #  Initial Gaussian:
    j=0
    wave0 = line_data['lambda'][j]
    line_center = wave0*scale_factor
    model_name = np.str(line_data['indx'][j])

    # Traditional astropy 1D Gaussian fit:
    joint_model = models.Gaussian1D(amplitude=amplitude_init,mean=line_center,
                                    stddev=stddev_init,
                                    name=model_name)
    # Set constraints on how much the central value can vary
    mean_bounds.append([line_center*mean_bounds_scale[0],
                         line_center*mean_bounds_scale[1]])

    #  Loop through the remaining lines:
    for j in np.arange(1,np.size(line_data)):
        wave0 = line_data['lambda'][j]
        line_center = wave0 * scale_factor
        model_name = np.str(line_data['indx'][j])

        # Traditional astropy 1D Gaussian fit:
        joint_model += models.Gaussian1D(amplitude=amplitude_init,mean=line_center,
                                  stddev=stddev_init,
                                  name=model_name)
        # Set constraints on how much the central value can vary
        mean_bounds.append([line_center * mean_bounds_scale[0],
                            line_center * mean_bounds_scale[1]])

    # Now we have to loop through the same models, applying the remaining bounds:
    for k in np.arange(0, np.size(line_data)):
        joint_model[k].bounds['amplitude'] = amplitude_bounds
        joint_model[k].bounds['mean'] = (mean_bounds[k][0],mean_bounds[k][1])
        joint_model[k].bounds['stddev'] = stddev_bounds

        # Tie some parameters together:
        if line_data['mode'][k] == 't33':
            joint_model[k].stddev.tied = tie_sigma_4862
        elif line_data['mode'][k] == 't35':
            joint_model[k].stddev.tied = tie_sigma_4862
            # joint_model[k].stddev.tied = tie_sigma_5008
        elif line_data['mode'][k] == 't45':
            joint_model[k].stddev.tied = tie_sigma_6585

        # Tie amplitudes of doublets
        if line_data['line'][k] == 'd35':
            joint_model[k].amplitude.tied = tie_ampl_5008   # 4959/5008
        if line_data['line'][k] == 'd45':
            joint_model[k].amplitude.tied = tie_ampl_6585   # 6549/6585

        # Finally, tie wavelengths of OII 3727, 3729, OIII 4364, 4960 to 5008
        if line_data['indx'][k] == '16':
            joint_model[k].mean.tied = tie_mean_3727_5008
        if line_data['indx'][k] == '17':
            joint_model[k].mean.tied = tie_mean_3729_5008
        if line_data['indx'][k] == '29':
            joint_model[k].mean.tied = tie_mean_4364_5008
        if line_data['indx'][k] == '34':
            joint_model[k].mean.tied = tie_mean_4960_5008

    # Standard astropy fit:
    fitter = fitting.LevMarLSQFitter()
    fitted_model = fitter(joint_model,wave, flux,
                          weights=1. / err**2,
                          maxiter=500)
    fitted_flux = fitted_model(wave)
    # Print the reason the fitter stops:
    fitter.fit_info['message']

    # Plot the results
    if do_plot:
        plt.clf()
        plt.plot(wave,flux,drawstyle='steps-mid',linewidth=2)
        plt.plot(wave,fitted_flux,color='orange',linewidth=2)

    ##### Create integrated fluxes and errors
    # The integration range is over +/-stddev * int_delta_factor
    int_delta_factor = _define_integration_delta()

    for j in np.arange(np.size(line_data)):
        # Calculate integrated fluxes, errors; deal with the blended O II 3727/3729 doublet
        if line_data[j]['name'] == '[OII]3727':
            # Calculate the integrated fluxes and errors
            iflux, ierr = integrate_line_flux(wave, flux, err, fitted_model[j].mean.value,
                                              fitted_model[j].stddev * int_delta_factor,
                                              line3727=True)
        elif line_data[j]['name'] == '[OII]3729':
            # pdb.set_trace()
            # For 3729, use the flux derived for 3726
            iflux = output_table['fit_flux'][j - 1]
            # For 3729, use its own error. This is appropriate for the fitted errors of both liness
            crap, ierr = integrate_line_flux(wave, flux, err, fitted_model[j].mean.value,
                                             fitted_model[j].stddev * int_delta_factor)
        else:
            # Calculate the integrated fluxes and errors
            iflux, ierr = integrate_line_flux(wave, flux, err, fitted_model[j].mean.value,
                                              fitted_model[j].stddev * int_delta_factor)


        redshift_out = (fitted_model[j].mean/line_data[j]['lambda'] - 1.)
        fitted_flux_out = np.sqrt(2.*np.pi)*fitted_model[j].amplitude*fitted_model[j].stddev

        if j == 0:
            # Define and construct the initial table to hold the results
            output_col_names, output_format, output_dtype = _define_output_table()

            output_data = [[line_data[j]['name']],
                           [line_data[j]['ion']],
                           [line_data[j]['lambda']],
                           [line_data[j]['indx']],
                           [line_data[j]['mode']],
                           [fitted_model[j].mean.value],
                           [fitted_model[j].amplitude.value],
                           [fitted_model[j].stddev.value],
                           [redshift_out], [fitted_flux_out],
                           [iflux], [ierr]]
            output_table = Table(output_data, names=output_col_names,dtype=output_dtype)
        else:
            output_table.add_row([line_data[j]['name'],line_data[j]['ion'],
                              line_data[j]['lambda'], line_data[j]['indx'],line_data[j]['mode'],
                              fitted_model[j].mean.value,fitted_model[j].amplitude.value,
                              fitted_model[j].stddev.value,
                              redshift_out, fitted_flux_out, iflux, ierr])

    # Set the output format of the results table:
    colnames = output_table.colnames
    for j in np.arange(np.size(colnames)):
             output_table[colnames[j]].format = output_format[j]

    # Set up the spectral table:
    spec_table = Table([wave,flux,err,fitted_flux],
                       names=['wave','flux','err','spec_fit'])

    # Write summary FITS files
    file_base = spec_file.rsplit('.')[0]
    table_file = file_base+'.HIIFitTable.fits'
    fit_file = file_base+'.HIIFitSpec.fits'

    output_table.write(table_file,overwrite=True)
    spec_table.write(fit_file,overwrite=True)

    return output_table

def _define_integration_delta(delta=3.):
    return delta

def _define_output_table():
    output_col_names = ['name', 'ion', 'lambda0', 'line_index',
                        'fit_mode', 'fit_lambda', 'fit_amplitude',
                        'fit_stddev', 'fit_redshift', 'fit_flux',
                        'int_flux', 'int_err']
    output_format = ['', '', '0.3f', '', '', '0.3f', '0.3g',
                     '0.2f', '0.7f', '0.3g', '0.3g', '0.3g']
    output_dtype = ['S11', 'S7', '<f8', 'S3', 'S3', '<f8', '<f8', '<f8', '<f8',
                    '<f8', '<f8', '<f8']

    return output_col_names,output_format,output_dtype


def integrate_line_flux(wave, flux, err, center, delta, line3727=False):
    from scipy.integrate import simps

    if line3727:
        integrate_me = ((wave >= center - delta) &
                        (wave <= center + delta + 3.))
    else:
        integrate_me = (np.abs(wave - center) <= delta)

    pix_size = np.abs(np.median(np.roll(wave[integrate_me],-1)-wave[integrate_me]))

    # Integrate the flux
    int_flux = simps(flux[integrate_me])*pix_size

    # Integrate the variance. Center must be subtracted to get the squared delta lam values right.
    int_err =  simps(err[integrate_me]**2)*pix_size**2
    int_err = np.sqrt(int_err)

    # Include 2% flux uncertainty
    int_err += 0.02*int_flux

    return int_flux, int_err

def prep_pyMCZ():
    # TODO Create pyMCZ input file from fitting / integrated results
    print('crap')