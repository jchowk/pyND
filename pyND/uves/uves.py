from __future__ import print_function, absolute_import, division, unicode_literals

def uves_log(filespec="ADP*.fits",
            outputfilebase="UVESdatalog",
            rawspec=False,
            noFITS = False,
            browser=False):
    """
    Output a quick overview of *reduced* UVES data in a directory

    :param filespec: optional regex code for finding files (def: *.fits)
    :param outputfilebase: optional output table basename (def: UVESdatalog)
    :param rawspec: Raw data (default: False)? Defult keys for reduced data are ['OBJECT', 'TEXPTIME','WAVELMIN', 'WAVELMAX','SNR', 'SPEC_RES','DATE-OBS']. For Raw data
    :param noFITS: suppress writing FITS table (def: False)
    :param browser: show results in browser (default: False)
    :return: astropy.Table log
    """

    import numpy as np
    from astropy.io import fits
    from astropy.table import Table, Column
    import glob

    # Following http://stackoverflow.com/questions/21583647/reading-headers-from-multiple-fits-file-from-the-same-directory

    # Default
    if rawspec:
        # if (filespec=="ADP*.fits"): filespec='UVES*.fits'
        keys = ['OBJECT', 'EXPTIME','HIERARCH ESO INS PATH',
                'HIERARCH ESO DPR CATG','HIERARCH ESO DPR TYPE',
                'DATE-OBS']
        ktypes = ['<U25','float','<U25',
                  '<U25','<U25',
                  '<U25']
    else:
        keys = ['OBJECT', 'TEXPTIME','WAVELMIN', 'WAVELMAX',
                'SNR', 'SPEC_RES','DATE-OBS']
        ktypes = ['<U25','float','float','float',
                'float','float','<U25']

    # Test for files matching the search term; bail if none available.
    fls = glob.glob(filespec)
    if np.shape(fls)[0] == 0:
        print("ERROR: No files matching the search term.")
        return

    # Setup for getting header keywords
    hdu = 0
    values = []
    fitsNames = []

    # dir = './'
    #for fitsName in glob.glob(dir + '*.fits'):
    for fitsName in fls:
        # Pull the header infor right from the file
        header = fits.getheader(fitsName, hdu)
        values.append([header.get(key) for key in keys])

        # if you want the fits file name only without the full path then
        # fitsNames.append(os.path.split(fitsName)[1])
        fitsNames.append(fitsName)

    #####
    # Create the output table:

    # if noPandas:
    ###############
    # Create a table container.
    # One trick is to use the data types in the first "values" to let
    # astropy guess datatypes. To use this trick, you need to specify the
    # column names in the table

    # if np.isin('OBJECT',keys):
    #     indx = np.where(np.array(keys) == 'OBJECT')[0][0]
    #     objlen = np.zeros(len(values),dtype='int')
    #     for i in np.arange(len(values)):
    #         objlen[i] = len(values[i][indx])

    # Old approach that doesn't keep the right str length for objects:
    # row0 = dict(zip(keys, values[0]))
    # t = Table([row0], names=keys)
    # for i in range(1, len(values)):

    # Create the base table:
    t = Table(names=np.chararray.replace(keys,'HIERARCH ESO ',''),dtype=ktypes)

    # now add all the other rows. again, because dict didn't preserve column
    # order, you have to repeat the dict here.
    for i in range(len(values)):
        t.add_row(values[i])

    # add the filenames column
    new_column = Column(name='fitsName', data=fitsNames)
    t.add_column(new_column)

    # # Prefer Pandas data frame approach, since it keeps the right length
    # #  on strings.
    # noPandas = False
    # try:
    #     import pandas as pd
    # except ModuleNotFoundError:
    #     noPandas = True
    #
    # if ~noPandas:
    #     pdfr = pd.DataFrame(np.transpose(values),index=keys)
    #     if np.isin('OBJECT',keys):
    #         objcol = Table.from_pandas(pdfr.transpose())['OBJECT']
    #         t.replace_column('OBJECT',objcol)



    # np.transpose(pdfr).to_csv('example.txt', sep='|',index=False)
    # save the file
    # http://docs.astropy.org/en/stable/table/io.html
    t.write(outputfilebase+'.txt', format='ascii.fixed_width',
            delimiter='|', overwrite=True)
    if ~noFITS:
        t.write(outputfilebase+'.fits',overwrite=True)


    if browser:
        t.show_in_browser()


# if __name__ == "__main__":
#     main()

    return t

def simple_coadd(uves_table=None, outputbase=None, airtovac=True):
    """simple_coadd(uves_table=None, outputbase=None,airtovac=True):

    Combine multiple UVES exposures obtained with the same set-up into a single spectrum, weighting by the inverse variance of the input data.

    For multiple set-ups, outputs an individual file for each unique wavelength range.
    """

    import numpy as np
    from astropy.table import Table
    from linetools.spectra.xspectrum1d import XSpectrum1D

    if uves_table==None:
    #   There is no input table; let's create one.
        uves_table = uves_log()

    # Select unique wavelengths to identify the setups
    # Table should be the output from UVES log or at least contain wavelengths and filenames.
    uniq_setup, uniq_indeces, uniq_inverse = np.unique(
        np.around(uves_table['WAVELMIN'], 1),
        return_index=True, return_inverse=True)

    num_setups = np.size(uniq_indeces)

    # Loop over the setups.
    for j in np.arange(num_setups):
        # Work out the setups that are to be coadded
        setup_obs = (np.where(uniq_inverse == j))[0]
        num_setup_obs = np.size(setup_obs)

        # Hold object name
        setup_obj = uves_table[uniq_indeces[j]]['OBJECT']

        # Loop over the files in this setup
        for k in np.arange(num_setup_obs):
            # inspec = fits.getdata(uves_table[setup_obs[k]]['fitsName'])
            specfile=uves_table[setup_obs[k]]['fitsName']
            inspec = Table.read(specfile)

            # Load the spectrum with XSpectrum1D:
            inspec = XSpectrum1D.from_file(specfile)
            inspec.meta['airvac'] = 'air'

            # Unless user requests, we transform to vacuum.
            if airtovac == True:
                inspec.airtovac()

            # Check for sig = 0:
            badErrors = (inspec.sig == 0)
            inspec.flux[badErrors] = np.nan
            inspec.sig[badErrors] = np.nan
            inspec.ivar[badErrors] = np.nan

            # Set up the arrays if it's the first spectrum.
            if k==0:
                out_file_list = specfile
                out_wave = inspec.wavelength.value
                # Do weighted average; calculate the inverse variance and the inv var-weighted flux.
                out_inv_variance = inspec.ivar.value
                out_flux_weighted = \
                    inspec.flux.value*inspec.ivar.value
            # For subsequent spectra:
            else:
                #Test for same array length:
                new_length = np.min([np.size(inspec.flux.value),
                                        np.size(out_flux_weighted)])

                out_file_list += ', '+specfile
                out_inv_variance[:new_length] += \
                    inspec.ivar.value[:new_length]
                out_flux_weighted[:new_length] += \
                    inspec.flux.value[:new_length]*inspec.ivar.value[:new_length]

        # Calculate the output weighted mean flux and error.
        out_flux = out_flux_weighted / out_inv_variance
        out_err = np.sqrt(1. / out_inv_variance)

        # Book-keeping: set up filename
        if outputbase == None:
            # If there is no base for the filenames, use the object name.
            outputfilename = uves_table[setup_obs[k]]['OBJECT']
        else:
            outputfilename=outputbase

        outputfilename = "{0}.uves.{1:0.0f}.{2:0.0f}.fits".format(outputfilename,
                   uves_table[setup_obs[k]]['WAVELMIN']*10.,
                   uves_table[setup_obs[k]]['WAVELMAX']*10.)
        # Set up the output table
        outputtable = Table([out_wave,out_flux,out_err],
                            names=['wave','flux','err'])
        # Load the spectrum in the standard way to get the units right:
        inspec = Table.read(uves_table[0]['fitsName'])
        # Assign units
        outputtable['wave'].unit = inspec['WAVE'].unit
        outputtable['flux'].unit = inspec['FLUX'].unit
        outputtable['err'].unit = inspec['ERR'].unit

        outputtable.meta = inspec.meta
        outputtable.meta['COMMENT'] = 'Simple coadd of '+out_file_list

        outputtable.write(outputfilename,overwrite=True)
        print('Wrote '+outputfilename+'.')

def full_coadd(uves_table=None, outputbase=None,
                wavelength_range=None, airtovac=True):
    """full_coadd(uves_table=None, outputbase=None,
                    wavelength_range=None, airtovac=True):

    Combine multiple UVES exposures UVES into a single spectrum, weighting by the inverse variance of the input data. Creates a single spectrum covering the full wavelength range of the individual exposures.

    Note: no guarantee the fluxes will be continuous if there are issues with the flux calibration of individual datasets.
    """

    import numpy as np
    import astropy.units as u
    from astropy.table import Table
    from linetools.spectra.xspectrum1d import XSpectrum1D

    if uves_table==None:
    #   There is no input table; let's create one.
        uves_table = uves_log()

    # Create minmax wavelengths for filename, wavelength array
    wave_minmax = []
    wave_minmax.append(np.min(uves_table['WAVELMIN'])*10.)
    wave_minmax.append(np.max(uves_table['WAVELMAX'])*10.)
    wave_minmax = np.around(wave_minmax)

    if wavelength_range == None:
        wavelength_range = wave_minmax

    # Store the string to make the filename:
    file_wave = '{0:0.0f}.{1:0.0f}'.format(wavelength_range[0],wavelength_range[1])

    # TODO: Allow for linear-space regular array
    # Create the wavelength array. For now uniform in log space
    log_steps = True
    if log_steps:
        wavelength_range = np.log10(wavelength_range)
        wavelength_step = 2.e-6 # log-space step
        out_wave = np.arange(wavelength_range[0],wavelength_range[1],
                             wavelength_step,dtype=np.float64)
        out_wave = 10.**out_wave

    # Create the holder arrays
    out_inv_variance = np.zeros(np.size(out_wave),dtype=np.float32)
    out_weighted_flux = np.zeros(np.size(out_wave),dtype=np.float32)

    # Loop over files:
    for j in np.arange(np.size(uves_table)):
        # Filename:
        specfile=uves_table[j]['fitsName']
        print('{0}: {1}'.format(j,specfile))

        # Load the spectrum with XSpectrum1D:
        inspec = XSpectrum1D.from_file(specfile)
        inspec.meta['airvac'] = 'air'

        # Unless user requests, we transform to vacuum.
        if airtovac == True:
            inspec.airtovac()

        # Check for sig = 0:
        badErrors = (inspec.sig == 0)
        inspec.flux[badErrors] = np.nan
        inspec.sig[badErrors] = np.nan
        inspec.ivar[badErrors] = np.nan

        # Interpolate the UVES fluxes onto the wavelength grid using the XSpectrum1D functionality
        grid_spec = inspec.rebin(out_wave*u.angstrom,do_sig=True,grow_bad_sig=True)

        # Add new fluxes, errors into output weighting, flux arrays
        # The inverse variance is stored in the XSpectrum1D class

        out_inv_variance += grid_spec.ivar.value
        out_weighted_flux += grid_spec.flux.value*grid_spec.ivar.value


    #####
    # Calculate the final output weighted mean flux and error by taking out the weighting
    out_flux = out_weighted_flux / out_inv_variance
    out_err = np.sqrt(1. / out_inv_variance)

    # Book-keeping: set up filename
    if outputbase == None:
        # If there is no base for the filenames, use the object name.
        outputbase = uves_table[0]['OBJECT']
    # Create final filename
    outputfilename = "{0}.uves.{1}.fits".format(outputbase,file_wave)

    # Set up the output table
    outputtable = Table([out_wave,out_flux,out_err],
                        names=['wave','flux','err'])

    # Load the spectrum in the standard way to get the units right:
    inspec = Table.read(uves_table[0]['fitsName'])
    # Assign units
    outputtable['wave'].unit = inspec['WAVE'].unit
    outputtable['flux'].unit = inspec['FLUX'].unit
    outputtable['err'].unit = inspec['ERR'].unit

    # TODO: Get the header information to save
    # The header:
    outputtable.meta = inspec.meta
    outputtable.meta['NCOMBINE'] = np.size(uves_table)
    outputtable.meta['HISTORY'] = 'Coadd of {0}'.format(uves_table[0]['fitsName'])
    for j in np.arange(1,np.size(uves_table)):
        outputtable.meta['HISTORY'] += ', {0}'.format(uves_table[j]['fitsName'])

    # TODO: Switch back to the XSpectrum1D format?
    # # Put the results into the XSpectrum1D format
    #out_spec = XSpectrum1D.from_tuple((out_wave*u.angstrom,out_flux,out_err),
    #                                   verbose=False)
    #out_spec.write(outputfilename)

    # Write the output
    outputtable.write(outputfilename,overwrite=True)
    print('Wrote '+outputfilename+'.')


def setwave(hdr):
    """ Generate wavelength array from a header

    Parameters
    ----------
    hdr : FITS header

    Returns
    -------
    wave : ndarray
      No units yet

    Adopted from linetools code.
    """
    import numpy as np

    # Parse the header
    npix = hdr['NAXIS1']
    crpix1 = hdr['CRPIX1'] if 'CRPIX1' in hdr else 1.
    crval1 = hdr['CRVAL1']

    cdelt1 = hdr['CDELT1']
    dc_flag = hdr['DC-FLAG']

    # Generate wave array
    wave = crval1 + cdelt1 * (np.arange(npix) + 1. - crpix1)
    if dc_flag == 1:
        wave = 10.**wave # Log

    return wave

def read_popler(datfil,savfil=''):
    """ Read in the output from a UVES_popler-reduction.

        Much of this code built from linetools examples.
    """
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table
    import os

    # hdulist = fits.open(os.path.expanduser(datfil), **kwargs)
    hdulist = fits.open(os.path.expanduser(datfil))
    hdr0 = hdulist[0].header

    # Check to see if this is a POPLER output file
    #  [code from lintools]
    poplerFile=False
    if 'history' in hdr0:
        for row in hdr0['history']:
            if 'UVES POst Pipeline Echelle Reduction' in row:
                poplerFile = True
    if poplerFile == False:
        print('File does not to appear to be output of UVES_popler.')
        return

    wave = setwave(hdr0)
    co = hdulist[0].data[3]         #  Continuum
    flux = hdulist[0].data[0] * co  #  Flux
    err = hdulist[0].data[1] * co   #  Error

    # Correct some _really_ bad error pixels
    bd=(err < -flux)
    err[bd] = np.median(err)*2.

    output_table = Table([wave,flux,err,co],names=['wave','flux','err','continuum'])

    if np.str.find(savfil,'.fits') > 0:
        output_table.write(savfil)

    return output_table
