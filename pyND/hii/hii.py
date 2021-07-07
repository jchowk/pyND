""" Extract HII region line fluxes
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import numpy as np
from astropy.io import ascii, fits
import astropy.constants as c

from astropy.table import Table, Column

from astropy.modeling.parameters import Parameter
from astropy.modeling.functional_models import Fittable1DModel
from astropy.modeling import models, fitting

from pyND.lbt.mods import air_to_vac, vac_to_air

import pdb

class GaussianEmission(Fittable1DModel):
    """One dimensional Gaussian model.

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
        x0 = (1.+redshift)*wave0
        return amplitude * np.exp(- 0.5 * (x -  x0)**2 / stddev**2)

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

def fit_lines_sherpa(spec_file, z_init=0., file_out=None,
                        do_plot=True, monte_carlo=False):
    """Fit an HII region spectrum using Sherpa package.    """

    from astropy.modeling.fitting import SherpaFitter
    import matplotlib.pyplot as plt
    from linetools.spectra.xspectrum1d import XSpectrum1D

    # Redshift scale:
    scale_factor = (1.+z_init)

    # Read in the spectrum. **ASSUME VACUUM WAVELENGTHS?**
    mods_spec = XSpectrum1D.from_file(spec_file)

    # Set up a convenient wavelength, flux, error arrays
    wave = mods_spec.wavelength.value
    flux = mods_spec.flux.value
    err = mods_spec.sig.value

    ###### ------ FOR TESTING!! ------
    ### To test this, let's constrain ourselves to only the wavelengths between ~Hbeta, OIII
    # g = np.where((wave >= 4000) & (wave <= 5400.))
    # wave = wave[g]
    # flux = flux[g]
    # err = err[g]

    # Load the data for the lines to be fit. Starts with MANGA line list, modified for MODS.
    line_data = get_linelist()
    # Exclude lines outside of the wavelength coverage.
    keep_lines = np.where((line_data['lambda'] <= np.max(wave)/scale_factor) &
                          (line_data['lambda'] >= np.min(wave)/scale_factor))[0]
    keep_line_index = line_data['indx'][keep_lines]

    ##### MODEL DEFINITIONS
    # Define initial parameters
    amplitude_init = 0.1*np.max(mods_spec.flux)
    stddev_init = 1.5

    amplitude_bounds, stddev_bounds, velocity_range = _define_bounds()

    # Calculate the redshift delta
    z_bounds_scale = (velocity_range/c.c.to('km/s').value)*scale_factor
    z_bounds = (z_init-z_bounds_scale, z_init+z_bounds_scale)

    # Define a joint model as the sums of Gaussians for each line
    #  Gaussian for first line:
    j=0
    wave0 = line_data['lambda'][j]
    line_center = wave0*scale_factor
    model_name = np.str(line_data['indx'][j])
    # Here we use a custom Gaussian class to  fix redshifts together
    joint_model = GaussianEmission(amplitude=amplitude_init,redshift=z_init,
                                   stddev=stddev_init,wave0=wave0,
                                   name=model_name)
    # The rest wavelength is not a free parameter:
    joint_model.wave0.fixed = True

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
        # The rest wavelength is not a free parameter:
        joint_model[k].wave0.fixed = True

    # TODO Get tied parameters to work.
    # Tie some parameters together, checking that reference lines
    #  are actually covered by the spectrum:
    for k in np.arange(0, np.size(line_data)):
        if (line_data['mode'][k] == 't33') & (np.in1d(33,keep_line_index)):
            joint_model[k].stddev.tied = _tie_sigma_4862
            joint_model[k].redshift.tied = _tie_redshift_4862
        elif (line_data['mode'][k] == 't35') & (np.in1d(35,keep_line_index)):
            joint_model[k].stddev.tied = _tie_sigma_5008
            joint_model[k].redshift.tied = _tie_redshift_5008
        elif (line_data['mode'][k] == 't45') & (np.in1d(45,keep_line_index)):
            joint_model[k].stddev.tied = _tie_sigma_6585
            joint_model[k].redshift.tied = _tie_redshift_6585
        elif (line_data['mode'][k] == 't46') & (np.in1d(46,keep_line_index)):
            joint_model[k].stddev.tied = _tie_sigma_6718
            joint_model[k].redshift.tied = _tie_redshift_6718

        # 3727/3729 lines:
        if line_data['name'][k] == '[OII]3727':
            joint_model[k].stddev.tied = _tie_sigma_3729
            joint_model[k].redshift.tied = _tie_redshift_3729

        # Tie amplitudes of doublets
        if (line_data['line'][k] == 'd35') & (np.in1d(35,keep_line_index)):
            joint_model[k].amplitude.tied = _tie_ampl_5008   # 4959/5008
        if (line_data['line'][k] == 'd45') & (np.in1d(45,keep_line_index)):
            joint_model[k].amplitude.tied = _tie_ampl_6585   # 6549/6585

    ##### FITTING
    # Sherpa model fitting from SABA package
    sfit = SherpaFitter(statistic='chi2', optimizer='levmar',
                        estmethod='confidence')
    sfit_lm = SherpaFitter(statistic='chi2', optimizer='neldermead',
                        estmethod='confidence')
    sfit_mc = SherpaFitter(statistic='chi2', optimizer='moncar',
                        estmethod='confidence')
    # Do the fit
    sfitted_model = sfit(joint_model, wave, flux, err = err)
    # Refine with different optimizer
    temp_model = sfitted_model.copy()
    sfitted_model = sfit_lm(temp_model, wave, flux, err = err)

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
    output_construct = 0

    for j in np.arange(np.size(line_data)):
        # Calculate integrated fluxes, errors;
        #  -- First test that the lines are in the range covered by data
        if np.in1d(j,keep_lines):
            mean_lambda = line_data[j]['lambda']*(1.+sfitted_model[j].redshift)

            #    deal with blended O II 3727/3729 doublet
            if line_data[j]['name'] == '[OII]3727':
                # Calculate the integrated fluxes and errors
                iflux, ierr = integrate_line_flux(wave, flux, err,
                                                  mean_lambda,
                                                  sfitted_model[j].stddev * int_delta_factor,
                                                  line3727=True)
            elif line_data[j]['name'] == '[OII]3729':
                # pdb.set_trace()
                # For 3729, use the flux derived for 3726
                iflux = output_table['int_flux'][j - 1]
                #iflux = 0.
                # For 3729, use its own error. This is appropriate for the fitted errors of both lines
                crap, ierr = integrate_line_flux(wave, flux, err, mean_lambda,
                                                 sfitted_model[j].stddev * int_delta_factor)
            else:
                # Calculate the integrated fluxes and errors
                iflux, ierr = integrate_line_flux(wave, flux, err, mean_lambda,
                                                  sfitted_model[j].stddev * int_delta_factor)

            redshift_out = (sfitted_model[j].redshift)[0]
            sfitted_flux_out = np.sqrt(2.*np.pi)*sfitted_model[j].amplitude*sfitted_model[j].stddev

            if output_construct == 0:
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
                               [iflux], [ierr], [iflux/ierr]]
                output_table = Table(output_data, names=output_col_names,dtype=output_dtype)

                output_construct = 1
            else:
                output_table.add_row([line_data[j]['name'],line_data[j]['ion'],
                                  line_data[j]['lambda'], line_data[j]['indx'],line_data[j]['mode'],
                                  mean_lambda,sfitted_model[j].amplitude.value,
                                  sfitted_model[j].stddev.value,
                                  redshift_out, sfitted_flux_out,
                                  iflux, ierr, iflux/ierr])

    # Set the output format of the results table:
    colnames = output_table.colnames
    for j in np.arange(np.size(colnames)):
        output_table[colnames[j]].format = output_format[j]

    # Set up the spectral table:
    spec_table = Table([wave,flux,err,sfitted_flux],
                       names=['wave','flux','err','spec_fit'])

    # Write summary FITS files
    if file_out is None:
        file_base = spec_file.strip('.fits')
    else:
        file_base = file_out

    table_file = file_base+'.HIIFitTable.fits'
    fit_file = file_base+'.HIIFitSpec.fits'

    output_table.write(table_file,overwrite=True)
    spec_table.write(fit_file,overwrite=True)

    return output_table


def fit_lines(spec_file, z_init=0., do_plot=True, file_out = None,
                no_tie_5008 = False):
    # from astropy.modeling import models, fitting
    from linetools.spectra.xspectrum1d import XSpectrum1D
    import matplotlib.pyplot as plt

    # Redshift scale:
    scale_factor = (1.+z_init)

    # Read in the spectrum. **ASSUME VACUUM WAVELENGTHS?**
    mods_spec = XSpectrum1D.from_file(spec_file)
    # mods_spec = Table.read(spec_file)

    # Set up a convenient wavelength, flux, error arrays
    wave = mods_spec.wavelength.value
    flux = mods_spec.flux.value
    err = mods_spec.sig.value
    # wave = mods_spec['wave']
    # flux = mods_spec['flux']
    # err = mods_spec['err']

    # Load the data for the lines to be fit. Starts with MANGA line list, modified for MODS.
    line_data = get_linelist()
    # Exclude lines outside of the wavelength coverage.
    keep_lines = np.where((line_data['lambda'] <= np.max(wave)/scale_factor) &
                          (line_data['lambda'] >= np.min(wave)/scale_factor))[0]
    keep_line_index = line_data['indx'][keep_lines]

    ##########
    # MODEL DEFINITION
    # Define a joint model as the sums of Gaussians for each line

    # Initial parameters
    amplitude_init = 0.1*np.max(flux)
    stddev_init = 1.5

    # Constraint parameters:
    amplitude_bounds, stddev_bounds, velocity_range = _define_bounds()

    # Start up the constraints on the mean wavelength  [separate bounds for each lambda0]
    mean_bounds_scale = (velocity_range / c.c.to('km/s').value) * np.array([-1., 1.]) + 1.
    mean_bounds = []

    ##### DEFINE THE MODEL:
    #  Initial Gaussian:
    j=0
    wave0 = line_data['lambda'][j]
    line_center = wave0*scale_factor
    model_name = np.str(line_data['indx'][j])

    # Traditional astropy 1D Gaussian model:
    joint_model = models.Gaussian1D(amplitude=amplitude_init,mean=line_center,
                                    stddev=stddev_init,
                                    name=model_name)
    # Set constraints on how much the central value can vary
    mean_bounds.append([line_center*mean_bounds_scale[0],
                         line_center*mean_bounds_scale[1]])

    #  Loop through the remaining lines to create their Gaussian model:
    for j in np.arange(1,np.size(line_data)):
        wave0 = line_data['lambda'][j]
        line_center = wave0 * scale_factor
        model_name = np.str(line_data['indx'][j])

        # Traditional astropy 1D Gaussian fit:
        joint_model += models.Gaussian1D(amplitude=amplitude_init,
                                             mean=line_center,
                                             stddev=stddev_init,
                                             name=model_name)
        # Set constraints on how much the central value can vary
        mean_bounds.append([line_center * mean_bounds_scale[0],
                            line_center * mean_bounds_scale[1]])

    # Now we have to loop through the same models, applying the
    #  remaining bounds. This includes tying parameters:
    for k in np.arange(0, np.size(line_data)):
        joint_model[k].bounds['amplitude'] = amplitude_bounds
        joint_model[k].bounds['mean'] = (mean_bounds[k][0],mean_bounds[k][1])
        joint_model[k].bounds['stddev'] = stddev_bounds

        # Tie some parameters together, checking that reference lines
        #  are actually covered by the spectrum:
        if (line_data['mode'][k] == 't33') & (np.in1d(33,keep_line_index)):
            joint_model[k].stddev.tied = _tie_sigma_4862
        elif (line_data['mode'][k] == 't35') & (np.in1d(35,keep_line_index)):
            #joint_model[k].stddev.tied = _tie_sigma_4862
            joint_model[k].stddev.tied = _tie_sigma_5008
        elif (line_data['mode'][k] == 't45') & (np.in1d(45,keep_line_index)):
            joint_model[k].stddev.tied = _tie_sigma_6585

        # Tie amplitudes of doublets
        if (line_data['line'][k] == 'd35') & (np.in1d(35,keep_line_index)):
            joint_model[k].amplitude.tied = _tie_ampl_5008   # 4959/5008
        if (line_data['line'][k] == 'd45') & (np.in1d(45,keep_line_index)):
            joint_model[k].amplitude.tied = _tie_ampl_6585   # 6549/6585

        # Finally, tie wavelengths of OII 3727, 3729, OIII 4364, 4960 to 5008
        if not no_tie_5008:
            if (line_data['indx'][k] == '16') & (np.in1d(35,keep_line_index)):
                joint_model[k].mean.tied = _tie_mean_3727_5008
            if (line_data['indx'][k] == '17') & (np.in1d(35,keep_line_index)):
                joint_model[k].mean.tied = _tie_mean_3729_5008
            if (line_data['indx'][k] == '29') & (np.in1d(35,keep_line_index)):
                joint_model[k].mean.tied = _tie_mean_4364_5008
            if (line_data['indx'][k] == '34') & (np.in1d(35,keep_line_index)):
                joint_model[k].mean.tied = _tie_mean_4960_5008

    # TODO Assess quality of emission line fits.
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


    # Test whether we've constructed output or not:
    output_construct = 0

    # Loop through the list of lines:
    #  --One could imagine only looping over lines covered by the spectrum,
    #     but this screws up the way we tie parameters.
    for j in np.arange(np.size(line_data)):
        # Only fit those lines that are covered by the spectrum
        if np.in1d(j,keep_lines):
            # Calculate integrated fluxes, errors; deal with the blended O II 3727/3729 doublet
            if line_data[j]['name'] == '[OII]3727':
                # Calculate the integrated fluxes and errors
                iflux, ierr = integrate_line_flux(wave, flux, err, fitted_model[j].mean.value,
                                                  fitted_model[j].stddev * int_delta_factor,
                                                  line3727=True)
            elif line_data[j]['name'] == '[OII]3729':
                # pdb.set_trace()
                # For 3729, use the flux derived for 3726
                iflux = output_table['int_flux'][j - 1]
                #iflux = 0.
                # For 3729, use its own error. This is appropriate for the fitted errors of both liness
                crap, ierr = integrate_line_flux(wave, flux, err, fitted_model[j].mean.value,
                                                 fitted_model[j].stddev * int_delta_factor)
            else:
                # Calculate the integrated fluxes and errors
                iflux, ierr = integrate_line_flux(wave, flux, err, fitted_model[j].mean.value,
                                                  fitted_model[j].stddev * int_delta_factor)


            redshift_out = (fitted_model[j].mean/line_data[j]['lambda'] - 1.)
            fitted_flux_out = np.sqrt(2.*np.pi)*fitted_model[j].amplitude*fitted_model[j].stddev

            if output_construct == 0:
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
                               [iflux], [ierr], [iflux/ierr]]
                output_table = Table(output_data, names=output_col_names,dtype=output_dtype)

                output_construct = 1
            else:
                output_table.add_row([line_data[j]['name'],line_data[j]['ion'],
                                  line_data[j]['lambda'], line_data[j]['indx'],line_data[j]['mode'],
                                  fitted_model[j].mean.value,fitted_model[j].amplitude.value,
                                  fitted_model[j].stddev.value,
                                  redshift_out, fitted_flux_out,
                                  iflux, ierr, iflux/ierr])

    # Set the output format of the results table:
    colnames = output_table.colnames
    for j in np.arange(np.size(colnames)):
             output_table[colnames[j]].format = output_format[j]

    # Set up the spectral table:
    spec_table = Table([wave,flux,err,fitted_flux],
                       names=['wave','flux','err','spec_fit'])

    # # Excise lines that weren't part of the fitting process:
    # output_table = output_table[keep_lines]

    # Write summary FITS files
    if file_out is None:
        file_base = spec_file.strip('.fits')
    else:
        file_base = file_out

    table_file = file_base+'.HIIFitTable.fits'
    fit_file = file_base+'.HIIFitSpec.fits'

    output_table.write(table_file,overwrite=True)
    spec_table.write(fit_file,overwrite=True)

    return output_table


##
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
    # Include 2% flux variance
    int_err += (0.02*int_flux)**2
    # Return the err (=sqrt(variance))
    int_err = np.sqrt(int_err)

    return int_flux, int_err

def prep_pyMCZ(filebase, file_input=None, file_output=None,
                data=None, fitted_flux=False):
    from astropy.table import Table,Column
    import numpy as np

    # TODO Create a path to output both integrated, fitted fluxes in different "regions."
    # User chooses fitted or integrated flux
    if fitted_flux:
        flux_label = 'fit_flux'
    else:
        flux_label = 'int_flux'

    if data is None and file_input is not None:
        num_files = np.size(file_input)
        if num_files > 1:
            data = Table.read(file_input[0])
            for j in np.arange(1,num_files):
                temp = Table.read(file_input[j])
                data = [data,temp]
            num_regions = num_files
        else:
            data = [Table.read(file_input)]
            num_regions = 1
    elif data is None:
        filename=filebase+'.HIIFitTable.fits'
        data = [Table.read(filename)]
        num_regions = 1
    else:
        data = [data]
        if np.shape(np.shape(data))[0] == 1:
            num_regions = 1
        else:
            num_regions = (np.shape(data))[0]

    # Define some information for pyMCZ input files:
    header_text = ';# galnum,[OII]3727,Hg,Hb,[OIII]4959,[OIII]5007,[OI]6300,Ha,[NII]6584,' \
                  '[SII]6717,[SII]6731,[SIII]9069,[SIII]9532'
    mcz_lines = header_text.split(',')[1:]

    hii_labels = '[OII]3727,Hg4341,Hb4862,[OIII]4960,[OIII]5008,[OI]6302,' \
                 'Ha6564,[NII]6585,' \
                  '[SII]6718,[SII]6732,[SIII]9071,[SIII]9533'
    hii_labels = hii_labels.split(',')
    num_lines = np.size(hii_labels)

    mcz_flux = (np.full(num_regions*num_lines,np.nan,
                        dtype=np.float32)).reshape(num_regions,num_lines)
    mcz_err = (np.full(num_regions*num_lines,np.nan,
                       dtype=np.float32)).reshape(num_regions, num_lines)
    flux_scalefactor = 1.e16
    region_label = np.full(num_regions,'_.',dtype='S2')

    for k in np.arange(num_regions):
        region = data[k]
        region_label[k] = region_label[k].replace('_',np.str(k+1))
        for j in np.arange(num_lines):
            gg = np.where(region['name'] == hii_labels[j])
            # pdb.set_trace()
            if np.size(gg) == 1:
                if region['significance'][gg] >= 3.:
                    mcz_flux[k,j] = region[flux_label][gg]*flux_scalefactor
                    mcz_err[k,j] = region['int_err'][gg]*flux_scalefactor

            # Fix OII line flux to include both members of the doublet
            if hii_labels[j] == '[OII]3727' and fitted_flux is True:
                # The OII lines are summed together already for the integrated data.
                # For the _fitted_ values, the total error is that of the 3727 line, but
                # the fluxes need to be summed to be appropriate for the SEL methods.
                oo = np.where(region['ion'] == '[OII]')[0]
                oxy2flux = (region[flux_label][oo]).sum() * flux_scalefactor
                oxy2err = (region['int_err'][oo[0]] * flux_scalefactor)

                mcz_flux[k, j] = oxy2flux
                mcz_err[k, j] = np.sqrt(oxy2err)

    flux_table = Table(mcz_flux,names=mcz_lines)
    err_table = Table(mcz_err,names=mcz_lines)

    colnames = flux_table.colnames
    for j in np.arange(np.size(colnames)):
        flux_table[colnames[j]].format = '0.3f'
        err_table[colnames[j]].format = '0.3f'

    region_column = Column(region_label,name=' galnum')
    flux_table.add_column(region_column,index=0)
    err_table.add_column(region_column,index=0)

    # Set the output filenames:
    if file_output is not None:
        file_meas = file_output.strip('.fits') + '_meas.txt'
        file_err = file_output.strip('.fits') + '_err.txt'
    else:
        file_meas = filebase + '_meas.txt'
        file_err = filebase + '_err.txt'

    # Write output files
    flux_table.write(file_meas, format='ascii.commented_header',
                     delimiter='\t', comment=';# ', overwrite=True)
    err_table.write(file_err, format='ascii.commented_header',
                     delimiter='\t', comment=';# ', overwrite=True)


########## HELPER FUNCTIONS ##########
def _define_integration_delta(delta=3.):
    return delta

def _define_bounds():
    # Set parameters constraints
    amplitude_bounds = (0.,1.e-12)  # Amplitudes must be >=0.
    stddev_bounds = (1.0,3.00)       # stddev between 0.5 and 5 Ang

    # Redshift constraints
    velocity_range = 750.   # Velocity range to explore about center

    return amplitude_bounds, stddev_bounds, velocity_range

##
def _define_output_table():
    output_col_names = ['name', 'ion', 'lambda0', 'line_index',
                        'fit_mode', 'fit_lambda', 'fit_amplitude',
                        'fit_stddev', 'fit_redshift', 'fit_flux',
                        'int_flux', 'int_err','significance']
    output_format = ['', '', '0.3f', '', '', '0.3f', '0.3g',
                     '0.2f', '0.7f', '0.3g', '0.3g', '0.3g','0.1f']
    output_dtype = ['S11', 'S7', '<f8', 'S3', 'S3', '<f8', '<f8', '<f8', '<f8',
                    '<f8', '<f8', '<f8','<f8']

    return output_col_names,output_format,output_dtype

#####
## Create the functional ties for parameters.
# (0, '[OII]3727', 16)
# (1, '[OII]3729', 17)
# (10, 'Hb4862', 33)
# (12, '[OIII]5008', 35)
# (15, 'Ha6564', 44)
# (16, '[NII]6585', 45)
# (17, '[SII]6718', 46)
# (18, '[SII]6732', 47)
# (22, '[SIII]9071', 60)
# (23, '[SIII]9533', 61)

def _tie_sigma_3729(model):
    # Tie the dispersions to O II 3729
    return model.stddev_2

def _tie_redshift_3729(model):
    # Tie the dispersions to O II 3729
    return model.redshift_2

def _tie_sigma_4862(model):
    # Tie the dispersions to Hbeta 4862
    return model.stddev_10

def _tie_redshift_4862(model):
    # Tie the redshift to Hbeta 4862
    return model.redshift_10

def _tie_sigma_5008(model):
    # Tie the dispersions to O III 5008
    return model.stddev_12

def _tie_redshift_5008(model):
    # Tie the redshift to O III 5008
    return model.redshift_12

def _tie_ampl_5008(model):
    #Tie 4959 flux to that of 5008
    return model.amplitude_12*0.350

def _tie_sigma_6585(model):
    # Tie the dispersions to N II 6585
    return model.stddev_16

def _tie_redshift_6585(model):
    # Tie the redshift to N II 6585
    return model.redshift_16

def _tie_ampl_6585(model):
    #Tie 6549 flux to that of 6585
    return model.amplitude_16*0.340

def _tie_sigma_6718(model):
    # Tie the dispersions to S II 6718
    return model.stddev_17

def _tie_redshift_6718(model):
    # Tie the redshift to S II 6718
    return model.redshift_17

def _tie_mean_3727_5008(model):
    # Tie the redshift to O III 5008
    return model.mean_12 * 3727.092 / 5008.240

def _tie_mean_3729_5008(model):
    # Tie the redshift to O III 5008
    return model.mean_12 * 3729.875 / 5008.240

def _tie_mean_4364_5008(model):
    # Tie the mean of 4363 auroral line to O III 5008
    return model.mean_12 *  4364.435 / 5008.240

def _tie_mean_4960_5008(model):
    # Tie the mean of 4363 auroral line to O III 5008
    return model.mean_12 *  4960.295 / 5008.240
