import numpy as np
from numpy import mean
from pyND.plotting import plotzero,plotaxes

def get_line(input_line, verbose=False, closest=True, line_list='ISM'):
    """Extract spectral line information from linetools lists.

    :param input_line -  linetools-style input, either symbol ('CIV 1548') or wavelength (1548.1)
    :param verbose(=False) - Print everything?
    :param closest(=True) - Allow an imprecise match [recommended]
    :param line_list(='ISM') - Which linetools line list to use [e.g., Strong, Gal, H2, etc.]
    :return: linetools line information.
    """
    from linetools.lists.linelist import LineList

    # Load in the linetools line lists.
    line_list = LineList(line_list, verbose=verbose, closest=closest)
    user_line = line_list[input_line]

    # Check that we've actually got a line:
    if not user_line:
        ## BAIL!!
        print('No line information.')

        return False
    else:
        return user_line

def get_fvals(input_line, log_lambdaf=False, wavelength=False,line_list='ISM'):
    """Extract f-values from linetools lists.

    :param input_line: linetools-style input, either symbol ('CIV 1548') or wavelength (1548.1)
    :param log_lambdaf: If True, return log lam*f instead of f
    :param wavelength:  If True, return also the wavelength.
    :return:
    """

    # Load in the linetools line lists.
    user_line = get_line(input_line,
            verbose=False, closest=True, line_list=line_list)

    if not user_line:
        ## BAIL!!
        print('No line information.')
        return False
    else:
        # Determine f-value for transition:
        fval = user_line['f']
        # Determine wavelength for transition:
        wave_out = user_line['wrest'].value
        # Log lambda*f
        lf = np.log10(wave_out*fval)

    def _ret():
        if log_lambdaf:
            # _rv = [round(lf,3)]
            _rv = round(lf, 3)
        else:
            _rv = round(fval, 4)

        if wavelength:
            _rv = _rv, wave_out
            # _rv.append(wave_out)
        # return tuple(_rv)
        return _rv

    return _ret()


def sensitivity(input_line, snr, fwhm, bvalue=False, instrument='COS',
                return_results=False):
    """Calculate EW, column densities limits achievable for assumed observational parameters/results.
    	Following Wakker et al. (1996, ApJ, 473, 834), which is based on (Kaper et al. 1966, Bull.
    	Astron. Inst. Netherlands, 18, 465).

    :param input_line: linetools-style line ID [e.g., 'CIV 1548']
    :param snr: SNR at line center.
    :param fwhm: FWHM in km/s
    :param bvalue: Input is b-value rather than FWHM: True / False [Default: False]
    :param instrument: cos,fuse,e140h,e140m,e230h,e230m,hires,mike,uves,uves03,mods03,mods06,mods,other
    :param return_results: If True, return a dict of results.     [Default: False]
    :return: wavelength, SNR, FWHM
    """
    # import
    import numpy as np
    from linetools.lists.linelist import LineList

    import astropy.units as u
    import astropy.constants as c

    _opt = (bvalue, instrument)

    # TODO Flag for S/N per resel (snr_per_resel = False)


    # Old code; doesn't take advantage of units (?)
    lightspeed = c.c.to('km / s').value

    # Load in the linetools line lists.
    line_list = LineList('ISM', verbose=False)
    user_line = line_list[input_line]

    # Check that we've actually got a line:
    if not user_line:
        ## BAIL!!
        print('No line information.')
        bad = True
    else:
        bad=False

    # Determine f-value for transition:
    fval = user_line['f']
    # Determine wavelength for transition:
    wave_out = user_line['wrest'].value
    # Log lambda*f
    lf = np.log10(wave_out*fval)

    # bad_instrument:
    ##Was an instrument specified?
    # if instrument is None:
    #     instrument = ' '
    #     print('No (or inappropriate) instrument specified.')
    #     read('Which instrument (cos,fuse,e140h,e140m,e230h,e230m,hires,mods03,mods06,mike,other)? ', instrument)
    #     print()

    ##Look up Instrument parameters:
    _expr = instrument.lower()
    if _expr == 'fuse':
        delta = 1.9 * 3.
        lsfwidth = 20.
    elif _expr == 'e140h':
        delta = lightspeed / 228000.
        lsfwidth = lightspeed / 114000.
    elif _expr == 'e140m':
        delta = lightspeed / 91700.
        lsfwidth = lightspeed / 45800.
    elif _expr == 'e230h':
        delta = lightspeed / 228000.
        lsfwidth = lightspeed / 114000.
    elif _expr == 'e230m':
        delta = lightspeed / 60000.
        lsfwidth = lightspeed / 30000.
    elif _expr == 'hires':
        delta = 2.
        lsfwidth = lightspeed / 45000.
    elif _expr == 'uves':
        delta = 0.78
        lsfwidth = lightspeed / 67000.
    elif _expr == 'uves03':
        delta = 0.78
        lsfwidth = lightspeed / 120000.
    elif _expr == 'mods03':
        delta = 30.
        lsfwidth = lightspeed / 4000.
    elif _expr == 'mods06':
        delta = 30.
        lsfwidth = lightspeed / 2000.
    elif _expr == 'mods':
        delta = 30.
        lsfwidth = lightspeed / 2000.
    elif _expr == 'cos':
        delta = 2.2
        lsfwidth = lightspeed / 18000.
    elif _expr == 'mike':
        if wave_out < 5000:
            delta = 0.02 / wave_out * lightspeed
            lsfwidth = lightspeed / 30000.
        else:
            delta = 0.05 / wave_out * lightspeed
            lsfwidth = lightspeed / 30000.
    elif _expr == 'other':
        ##    'other' :  begin
        delta = input('Enter instrument pixel spacing [km/s]: ')
        lsfwidth = input('Enter instrumental line spread function [fwhm, in km/s]: ')
        ##
    else:
        print('Assuming COS G130M/G160M')
        delta = 2.2
        lsfwidth = lightspeed / 18000.

    if (bvalue is True):
        fwhm = np.sqrt(fwhm ** 2. + (lsfwidth / 1.67) ** 2.)
        eqwidth = 6.5e-3 * wave_out * (delta * fwhm) ** (0.5) / snr
    else:
        fwhm = np.sqrt(fwhm ** 2. + lsfwidth ** 2.)
        eqwidth = 5.0e-3 * wave_out * (delta * fwhm) ** (0.5) / snr

    ColDens = 3. * eqwidth * (1.13e17 / (wave_out * 10 ** (lf)))

    # Allow for the case where we can't find a line value
    if bad:
        print("-------------------------------------")
        print('Eq.Width <  {0:0.2f} mA (1-sigma)'.format(eqwidth))
        print("-------------------------------------")
        print("    No f-value found")
    else:
        print("---------------------------------------")
        print(' Eq.Width <  {0:0.2f}     mA   (1-sigma)'.format(eqwidth))
        print('    N     <  {0:0.2g} cm**-2 (3-sigma)'.format(ColDens))
        print('  log N   <  {0:0.2f}         (3-sigma)'.format(np.log10(ColDens)))
        print("---------------------------------------")
        if bvalue is True:
            print('--> For f({0:0.3f}) = {1:0.3f}; b = {2:0.1f}'.format(wave_out, fval, fwhm))
        else:
            print('--> For f({0:0.3f}) = {1:0.3f}; FWHM = {2:0.1f} km/s'.format(wave_out, fval, fwhm))

    def _ret():
        _rv = [input_line, wave_out * u.angstrom,
               np.round(fval, 3), np.round(eqwidth * 3.e-3, 2) * u.angstrom, ColDens / u.cm ** 2]
        return tuple(_rv)

    if return_results:
        return _ret()


def sum_components(Col, Err, return_output=True, print_output=True):
    """
    output = sum_components(ColArray, ErrArray, return_output=True, print_output=True)

    :param Col: Array of input column density.
    :param Err: Array of input errors.
    :return: Output column, output error

    	11/10/97 -- Created by jch.
    	04/09/17 -- Translated to python.

    """

    import numpy as np

    n_params = 2
    sz1 = np.size(Col)

    err_scale = np.log10(np.exp(1))
    sum_col = 0
    sum_err = 0

    for i in np.arange(sz1):
        sum_col += 10.** (Col[i])
        sum_err += (10.**Col[i] * Err[i]/err_scale)**2

    sum_err = (sum_err) ** (0.5)

    log_sumCol = np.round(np.log10(sum_col),3)
    #log_sumErr = np.log10(1 + (sum_err / sum_col))
    log_sumErr = np.round((sum_err / sum_col)*err_scale,3)

    ## log_sumCol = alog10(total(10.0D**Col))
    ## log_sumErr = sqrt( total((10.0D**Col*2.3*Err)**2) )/total(10D**Col)/2.3


    if print_output:
        print("\t log N = {0} +/- {1}".format(log_sumCol, log_sumErr))

    def _ret():  return (log_sumCol, log_sumErr)

    if return_output:
        return _ret()


def logmean(log_data, log_err, verbose=True, return_straight=False):
    """Calculate the weighted average (and errors) of a series of values given in log-space.
    	log_data   = Array holding the log-space data.
    	log_err = Array holding the log-space errors.

    	11/01/97 -- Created: jch
    	04/10/17 -- Transferred from IDL to Python
    """

    import numpy as np

    def _ret():
        if return_straight == True:
            return (np.log10(best), (best_err/best)/np.log(10.))
        else:
            return (log_mean, log_mean_err)

    # Convert inputs to numpy arrays
    log_data = np.array(log_data)
    log_err = np.array(log_err)

    # Convert to linear data, errors
    linear_data = 10. ** log_data
    linear_err = linear_data*(log_err * np.log(10))

    sum_data = (linear_data/linear_err**2).sum()
    sum_err = (1./linear_err**2).sum()

    mean = sum_data / sum_err
    mean_err = (1 / sum_err) ** (0.5)

    log_mean = np.log10(mean)
    log_mean_err = (mean_err / mean) / np.log(10.)


    best = linear_data.mean()
    best_err = np.sqrt((linear_err**2).sum())

    if verbose:
        print("Weighted Mean:")
        print("{0:0.3g} +/- {1:0.4g}".format(mean,mean_err))
        print("{0:0.3f} +/- {1:0.3f}".format(log_mean,log_mean_err))

        print("Straight Mean:")
        print("{0:0.4g} +/- {1:0.4g}".format(best, best_err))
        print("{0:0.3f} +/- {1:0.3f}".format(np.log10(best), (best_err/best)/np.log(10.)))

    return _ret()



def correct_saturation(logN, err_logN, corr_max = 0.15, printresults=True):
    """Correct mildly-saturated Na(v) results for doublets following the recommendation in Savage & Sembach (1991). **This assumes a factor of 2x difference in f-values.

    CALLING SEQUENCE:
    sscorrect, logn_array, err(logn)_array

    INPUTS:
    logN     -- a two-element array holding the Na(v) values [strong, weak].
    err_logN -- a two-element array holding the Na(v) errors [strong, weak].

    OPTIONAL INPUTS:
    corr_max = 0.15     -- The assumed maximum value of valid corrections.
    printresults = True -- print results to the console.

    OUTPUTS:
    logNf     -- Final corrected column density.
    err_logNf -- Error in final column density (-2 == saturated).
    """

    # The difference in columns is
    diff=logN[1]-logN[0]

    # Check that the column densities aren't reversed.
    # Assume they are if the difference in columns is negative.
    if diff < 0:
        print, "The logN difference is negative! Did you put the strong line first? If so, there may be a problem, BUT ASSUMING THERE IS NONE...."

        logN = logN[::-1]
        err_logN = err_logN[::-1]
        diff=logN[1]-logN[0]

    # Calculate the correction using polynomial fit to SS1991 results.
    ssdiff=diff
    sscorrect=16.026*ssdiff**3 - 0.507*ssdiff**2 \
                + 0.9971*ssdiff + 5.e-5


    err_logNf = np.sqrt(err_logN[1]**2 + \
        ((48.078*ssdiff**2 - 1.014*ssdiff + 0.9971)*np.sqrt(err_logN[1]**2 \
        + err_logN[0]**2) )**2 )

    if sscorrect >= corr_max:
        print("Your correction exceeded the maximum correction, d(logN)_max = {0:0.3f}.".format(corr_max))
        print("Applying the maximum correction and assuming the final result is saturated.")

        logNf = logN[1]+corr_max
        err_logNf = -2
        if printresults:
            print("Final result: logN_final > {0:0.3f}".format(logNf))
    else:
        logNf=logN[1]+sscorrect
        if printresults:
            print("Final result: logN_final = {0:0.3f}+/-{1:0.3f}".format(logNf,err_logNf))

    return logNf, err_logNf
