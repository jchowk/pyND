So, I'll also attach one more thing, with no warranty and relatively
little support.  I use IDL a lot, including to analyze Cloudy results.
I've put together a little IDL routine that reads in all of the Cloudy
output files in a given directory (assuming they're the ones ending in
.out) and stores the resulting information in a FITS file.
 
cloudy_fits, output_filename, OUTPUT_STRUCTURE = cloudy, $
	INPUT_FILES = 'inputfilebase' 
 
All of the inputs are optional.  So, to make a FITS file from my
calculations (but using only those files that start with qso), I would
do the following:
 
IDL> cloudy_fits, 'cloudy_output.fits', output=cloud1, input='qso*' 
 
This would make the file cloudy_output.fits and populate the IDL
structure cloud1 with the resulting information.  Another way to do
the last bit after the fact is (if you have the IDL astronomy user's
library):
 
IDL> cloud1=mrdfits('cloudy_output.fits',1) 
 
