def get_line(input_line, verbose=False, closest=True):
    """
    output = abs_getline(input_line, verbose=False, closest=True)

    :param input_line: linetools-style input, either symbol ('CIV 1548') or wavelength (1548.1)  
    :return: linetools line information. 
    """
    from linetools.lists.linelist import LineList

    # Load in the linetools line lists.
    # TODO: allow user alteration of line list.
    line_list = LineList('ISM', verbose=verbose, closest=closest)
    user_line = line_list[input_line]

    # Check that we've actually got a line:
    bad = False
    if not user_line:
        ## BAIL!!
        print('No line information.')
        bad = True

    return user_line

def get_fvals(input_line, log_lambdaf=False, wavelength=False):
    """
    output = get_fvals(input_line,log_lambdaf=False,wavelength=False)

    :param input_line: linetools-style input, either symbol ('CIV 1548') or wavelength (1548.1)  
    :param log_lambdaf: If True, return log lam*f instead of f
    :param wavelength:  If True, return also the wavelength.
    :return: 
    """

    import numpy as np
    from linetools.lists.linelist import LineList

    # Load in the linetools line lists.
    user_line = get_line(input_line,verbose=False, closest=True)

    # Determine f-value for transition:
    fval = user_line['f']
    # Determine wavelength for transition:
    wave_out = user_line['wrest'].value
    # Log lambda*f
    lf = user_line['log(w*f)']

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


def sensitivity(input_line, snr, fwhm, bvalue=False, instrument='COS', return_results=False):
    """
        out = abs_sensitivity(input_line, snr, fwhm, bvalue=False, instrument='COS',return_results=False)

    	Calculate EW, column densities limits achievable for assumed observational parameters/results. 
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

    # Old code; doesn't take advantage of units (?)
    lightspeed = c.c.to('km / s').value

    # Load in the linetools line lists.
    line_list = LineList('ISM', verbose=False)
    user_line = line_list[input_line]

    # Check that we've actually got a line:
    bad = False
    if not user_line:
        ## BAIL!!
        print('No line information.')
        bad = True

    # Determine f-value for transition:
    fval = user_line['f']
    # Determine wavelength for transition:
    wave_out = user_line['wrest'].value
    # Log lambda*f
    lf = user_line['log(w*f)']

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
        print('    N     <  {0:0.2g} cm^-2 (3-sigma)'.format(ColDens))
        print('  log N   <  {0:0.2f}         (3-sigma)'.format(np.log10(ColDens)))
        print("---------------------------------------")
        if bvalue is True:
            print('--> For f({0:0.3f}) = {1:0.3f}; b = {2:0.1f}').format(wave_out, fval, fwhm)
        else:
            print('--> For f({0:0.3f}) = {1:0.3f}; FWHM = {2:0.1f} km/s').format(wave_out, fval, fwhm)

    def _ret():
        _rv = [input_line, wave_out * u.angstrom,
               np.round(fval, 3), np.round(eqwidth * 3.e-3, 2) * u.angstrom, ColDens / u.cm ** 2]
        return tuple(_rv)

    if return_results:
        return _ret()
