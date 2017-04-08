import glob, os
from astropy.io import fits
import numpy as np

from astropy.table import Table, Column, Row

def uves_log(filespec="*.fits", outputfile="UVESdatalog.txt", browser=False):
    """
    Output a quick overview of *reduced* UVES data in a directory
    
    :param filespec: optional regex code for finding files (def: *.fits) 
    :param outputfile: optional output table name (def: UVESdatalog.txt)
    :param browser: show results in browser (default: False)
    :return: astropy.Table log
    """

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