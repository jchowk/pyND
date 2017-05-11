""" Extract HII region line fluxes
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import numpy as np
from astropy.io import fits, ascii
import astropy.constants as c
import astropy.units as u

from astropy.table import Table, vstack, Column, MaskedColumn

from astropy.modeling.parameters import Parameter
from astropy.modeling.core import Fittable1DModel

import pdb


class GaussianEmission(Fittable1DModel):
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


def get_lines(use_mods_table=True):
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


def fit_flux(spec_file,z_init):
    from saba import SherpaFitter
    from astropy.modeling import models, fitting
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
    line_data = get_lines()
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
    amplitude_bounds = (0.,None)  # Amplitudes must be >=0.
    stddev_bounds = (0.5,5.0)       # stddev between 0.5 and 5 Ang

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

        # Traditional astropy 1D Gaussian fit:
        # joint_model += Gaussian1D(amplitude=amplitude_init,mean=line_center,
        #                                stddev=stddev_init,
        #                                name=model_name)
        # Set constraints on how much the central value can vary
        #mean_bounds.append([line_center * mean_bounds_scale[0], line_center * mean_bounds_scale[1]])

    # Now we have to loop through the same models, applying the bounds:
    for k in np.arange(0, np.size(line_data)):
        joint_model[k].bounds['amplitude'] = amplitude_bounds
        joint_model[k].bounds['redshift'] = z_bounds
        joint_model[k].bounds['stddev'] = stddev_bounds
        joint_model[k].wave0.fixed = True

        # Tie some parameters together: only works using GaussianEmission
        if line_data['mode'][k] == 't33':
            joint_model[k].stddev.tied = tie_sigma_4862
            joint_model[k].redshift.tied = tie_redshift_4862
        elif line_data['mode'][k] == 't35':
            joint_model[k].stddev.tied = tie_sigma_5008
            joint_model[k].redshift.tied = tie_redshift_5008
        elif line_data['mode'][k] == 't45':
            joint_model[k].stddev.tied = tie_sigma_6585
            joint_model[k].redshift.tied = tie_redshift_6585

        # Tie fluxes for doublets:
        if line_data['line'][k] == 'd35':
            joint_model[k].amplitude.tied = tie_ampl_5008
        if line_data['line'][k] == 'd45':
            joint_model[k].amplitude.tied = tie_ampl_6585

    # Constrain the wavelengths over which the fits are calculated (don't need the continuum)
    # Initialize the boolean indeces:
    fit_me = (wave < 0.)
    fit_delta = 25.  # A guess...
    for j in np.arange(0,np.size(line_data)):
        fit_me += (np.abs(wave/scale_factor - line_data['lambda'][j]) <= fit_delta)

    # Standard astropy fit:
    #fitter = fitting.LevMarLSQFitter()
    # fitted_model = fitter(joint_model,wave, flux,
    #                       weights=1. / err**2,
    #                       maxiter=500)
    # fitted_lines = fitted_model(wave)
    # fitter.fit_info['message']

    ##### FITTING
    # Sherpa model fitting from SABA package
    sfit = SherpaFitter(statistic='chi2gehrels', optimizer='levmar', estmethod='confidence')
    sfitted_model = sfit(joint_model,wave, flux,
                          err = err*1.1)
    # Work out the errors...
    #sfit.est_config['maxiters']=200
    #sfitted_err = sfit.est_errors(sigma=3)

    # Store the output parameters
    out_parameters = (sfitted_model.parameters).reshape(np.size(line_data),3)

    return

def fit_flux_free(spec_file,z_init):
    from astropy.modeling import models, fitting
    from linetools.spectra.xspectrum1d import XSpectrum1D

    def tie_sigma_4862(model):
        # Tie the dispersions to Hbeta 4862
        return model['33'].stddev

    def tie_sigma_5008(model):
        # Tie the dispersions to O III 5008
        return model['35'].stddev

    def tie_ampl_5008(model):
        #Tie 4959 flux to that of 5008
        return model['35'].amplitude*0.350

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
    line_data = get_lines()
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
    stddev_bounds = (1.,5.)
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
            joint_model[k].stddev.tied = tie_sigma_5008
        elif line_data['mode'][k] == 't45':
            joint_model[k].stddev.tied = tie_sigma_6585

        if line_data['line'][k] == 'd35':
            joint_model[k].amplitude.tied = tie_ampl_5008   # 4959/5008
        if line_data['line'][k] == 'd45':
            joint_model[k].amplitude.tied = tie_ampl_6585   # 6549/6585

    # Constrain the wavelengths over which the fits are calculated (don't need the continuum)
    # Initialize the boolean indeces:
    fit_me = (wave < 0.)
    fit_delta = 25.  # A guess...
    for j in np.arange(0,np.size(line_data)):
        fit_me += (np.abs(wave/scale_factor - line_data['lambda'][j]) <= fit_delta)

    # Standard astropy fit:
    fitter = fitting.LevMarLSQFitter()
    fitted_model = fitter(joint_model,wave, flux,
                          weights=1. / err**2,
                          maxiter=500)
    fitted_flux = fitted_model(wave)
    fitter.fit_info['message']

    # Store the output parameters
    out_parameters = (fitted_model.parameters).reshape(np.size(line_data),3)
    amplitude_out = []
    mean_out = []
    stddev_out = []
    redshift_out = []
    fitted_flux_out = []
    integrated_flux_out = []
    integrated_err_out = []

    for j in np.arange(np.size(line_data)):
        amplitude_out.append(fitted_model[j].amplitude.value)
        mean_out.append(fitted_model[j].mean.value)
        stddev_out.append(fitted_model[j].stddev.value)
        redshift_out.append(fitted_model[j].mean/line_data[j]['lambda'] - 1.)
        fitted_flux_out.append(np.sqrt(2.*np.pi)*fitted_model[j].amplitude*fitted_model[j].stddev)

        iflux,ierr = integrate_flux(wave,flux,err,mean_out[j],stddev_out[j]*4.)
        integrated_flux_out.append(iflux)
        integrated_err_out.append(ierr)
        



    return

def integrate_flux(wave,flux,err,center,delta):
    from scipy.integrate import simps

    integrate_me = (np.abs(wave - center) <= delta)
    # Integrate the flux
    int_flux = simps(flux[integrate_me],x=wave[integrate_me]-center)

    # Integrate the variance. Center must be subtracted to get the squared delta lam values right.
    int_err =  simps(err[integrate_me]**2,x=(wave[integrate_me]-center)**2)
    int_err = np.sqrt(int_err)

    return int_flux, int_err
