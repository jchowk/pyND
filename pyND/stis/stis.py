def plotx1d(inputfile,outputfile=None):
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import numpy as np

    # Load the data
    a=fits.open(inputfile)
    ahdr=a[0].header
    b=a[1].data

    # Scale the flux:
    scaleExponent = np.floor(np.log10(np.median(b['flux'])))
    b['flux']*=10.**(-scaleExponent)


    # How many orders to plot?
    num_orders = np.size(b['sporder'])

    # Axis labels:
    xlabelText = r'Wavelength [$\AA$]'
    ylabelText = r'Flux [$10^{{{0:0.0f}}}$ cgs]'.format(scaleExponent)

    ylabel_x = 0.04
    xlabel_y = 0.05


    # Object name:
    objectName = ahdr['targname']

    ##Output filename:
    if outputfile == None:
        outputfile = inputfile.replace('.fits','.pdf')

    num_columns=1
    num_rows = 4
    figsize=(8,10)


    # Set-up the multi-page PDF:
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(outputfile)

    # How many pages?
    num_pages = np.ceil(num_orders / num_rows)


    for pg in np.arange(num_pages):
        # Set up the plot grid
        fig=plt.figure(figsize=figsize)
        figspec = fig.add_gridspec(ncols=num_columns,
                        nrows=num_rows,hspace=0.25)

        for row in np.arange(num_rows):
            for col in np.arange(num_columns):
                spec_indx = (num_orders-1)-np.int(pg*(num_rows*num_columns)+(num_rows*col + row))

                if (spec_indx >= 0):
                    ax = fig.add_subplot(figspec[row,col])

                    wave=b['WAVELENGTH'][spec_indx][2:-2]
                    flux=b['FLUX'][spec_indx][2:-2]

                    ax.plot(wave,flux,drawstyle='steps-mid')
                    ax.axhline(0.,linewidth=1,color='k',linestyle='--')

                    # txt_lbl=objectName+' [index={0}: Order {1}]'.format(spec_indx,b['sporder'][spec_indx])
                    txt_lbl='[index={0}: Order {1}]'.format(spec_indx,b['sporder'][spec_indx])
                    ax.text(1.01, 0.5,txt_lbl,rotation=-90,
                            ha='left',va='center',color='red',fontsize=8,transform=ax.transAxes)

                    ax.set_xlim(np.min(wave),np.max(wave))


        if (pg == num_pages-1):
            ax.set_xlabel(xlabelText, fontsize=16)
            ax.set_ylabel(ylabelText, fontsize=16)
        else:
            # Apply the x,y axis labels
            fig.text(0.5,xlabel_y,xlabelText, fontsize=16,
                        ha='center')
            fig.text(ylabel_x,0.5,ylabelText,fontsize=16,
                        va='center',rotation='vertical')
        pp.savefig()
        plt.close()

    pp.close()
    plt.close()
