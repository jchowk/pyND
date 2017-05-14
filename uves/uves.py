def uves_log(filespec="ADP*.fits", outputfile="UVESdatalog.txt", browser=False):
    """
    Output a quick overview of *reduced* UVES data in a directory
    
    :param filespec: optional regex code for finding files (def: *.fits) 
    :param outputfile: optional output table name (def: UVESdatalog.txt)
    :param browser: show results in browser (default: False)
    :return: astropy.Table log
    """

    import glob, os
    from astropy.io import fits
    import numpy as np

    from astropy.table import Table, Column, Row

    # Keynames of interest
    keys = ['OBJECT', 'TEXPTIME', 'WAVELMIN', 'WAVELMAX', 'SNR', 'SPEC_RES','DATE-OBS']

    # Following http://stackoverflow.com/questions/21583647/reading-headers-from-multiple-fits-file-from-the-same-directory

    dir = './'
    hdu = 0
    # get header keyword values
    values = []
    fitsNames = []
    #for fitsName in glob.glob(dir + '*.fits'):
    for fitsName in glob.glob(filespec):
        # Pull the header infor right from the file
        header = fits.getheader(fitsName, hdu)
        values.append([header.get(key) for key in keys])

        # if you want the fits file name only without the full path then
        # fitsNames.append(os.path.split(fitsName)[1])
        fitsNames.append(fitsName)

    ###############
    # Create a table container.

    # One trick is to use the data types in the first "values" to let astropy guess
    # datatypes. To use this trick, you need to specify the column names in the
    # table

    row0 = [dict(zip(keys, values[0]))]
    t = Table(row0, names=keys)

    # now add all the other rows. again, because dict didn't preserve column order, you have to repeat
    # the dict here.
    for i in range(1, len(values)):
        t.add_row(values[i])

    # add the filenames column
    # t.add_column
    new_column = Column(name='fitsName', data=fitsNames)
    # t.add_column(new_column, 0)
    t.add_column(new_column)

    # save the file
    # http://docs.astropy.org/en/stable/table/io.html
    t.write(outputfile, format='ascii', delimiter='|', overwrite=True)

    if browser:
        t.show_in_browser()


# if __name__ == "__main__":
#     main()

    return t


def simple_coadd(uves_table=None, outputbase=None):

    import numpy as np
    from astropy.io import fits,ascii
    from astropy.table import Table,Column

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
            # Set up the arrays if it's the first spectrum.
            if k==0:
                out_file_list = specfile
                out_wave = inspec['WAVE'][0,:]
                # Do weighted average; calculate the inverse variance and the inv var-weighted flux.
                out_inv_variance = 1./(inspec['ERR'][0,:]**2)
                out_flux_weighted = inspec['FLUX'][0,:]/(inspec['ERR'][0,:]**2)
            # For subsequent spectra:
            else:
                #Test for same array length:
                if np.size(out_flux_weighted) > np.size(inspec['ERR'][0,:]):
                    new_length = np.size(inspec['ERR'][0,:])
                elif  np.size(out_flux_weighted) < np.size(inspec['ERR'][0,:]):
                    new_length = np.size(inspec['ERR'][0,:])
                else:
                    new_length = np.size(out_flux_weighted)

                out_file_list += ', '+specfile
                out_inv_variance[:new_length] += 1./(inspec['ERR'][0,:new_length]**2)
                out_flux_weighted[:new_length] += inspec['FLUX'][0,:]/(inspec['ERR'][0,:new_length]**2)

        # Calculate the output weighted mean flux and error.
        out_flux = out_flux_weighted / out_inv_variance
        out_err = np.sqrt(1. / out_inv_variance)

        # Book-keeping: set up filename
        if outputbase == None:
            #     If there is no base for the filenames, use the object name.
            outputfilename = uves_table[setup_obs[k]]['OBJECT']
        else:
            outputfilename=outputbase

        outputfilename = "{0}.uves.{1:0.0f}.{2:0.0f}.fits".format(outputfilename,
                                                           uves_table[setup_obs[k]]['WAVELMIN']*10.,
                                                           uves_table[setup_obs[k]]['WAVELMAX'] * 10.)
        # Set up the output table
        outputtable = Table([out_wave,out_flux,out_err],
                            names=['wave','flux','err'])
        outputtable['wave'].unit = inspec['WAVE'].unit
        outputtable['flux'].unit = inspec['FLUX'].unit
        outputtable['err'].unit = inspec['ERR'].unit

        outputtable.meta = inspec.meta
        outputtable.meta['COMMENT'] = 'Simple coadd of '+out_file_list

        outputtable.write(outputfilename,overwrite=True)
        print('Wrote '+outputfilename+'.')

def full_coadd(uves_table=None, outputbase=None, wavelength_range=None, air2vac=True):

    import numpy as np
    from astropy.io import fits,ascii
    import astropy.units as u
    from astropy.table import Table,Column
    from linetools.spectra.xspectrum1d import XSpectrum1D
    import pdb

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
        if air2vac == True:
            inspec.airtovac()

        # Check for sig = 0:
        badErrors = (inspec.sig == 0)
        inspec.flux[badErrors] = np.nan
        inspec.sig[badErrors] = np.nan
        inspec.ivar[badErrors] = np.nan

        # Interpolate the UVES fluxes onto the wavelength grid using the XSpectrum1D functionality
        grid_spec = inspec.rebin(out_wave*u.angstrom,do_sig=True,grow_bad_sig=True)

        # Work out where the newly-gridded spectrum goes in the output arrays
        sect = np.where((out_wave >= grid_spec.wvmin.value) & (out_wave <= grid_spec.wvmax.value))[0]

        # Add new fluxes, errors into output weighting, flux arrays
        # The inverse variance is stored in the XSpectrum1D class
        out_inv_variance[sect] += grid_spec.ivar.value
        out_weighted_flux[sect] += grid_spec.flux.value*grid_spec.ivar.value

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