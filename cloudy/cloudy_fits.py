"""
#
"""

from numpy import *

def float_after_text(string_input, text_value, number_digits):
##Find that float after the given text marker:
    n_params = 3
    def _ret():  return None
    
    if number_digits.size == 0:    
        number_digits = 7
    flagpos = strpos(string_input, text_value)
    output = array(strmid(string_input, flagpos + strlen(text_value), number_digits), copy=0).astype(float32)
    return output


def float_at_position(string_input, string_position, number_digits):
##Find that float after the given text marker:
    n_params = 3
    def _ret():  return None
    
    if number_digits.size == 0:    
        number_digits = 10
    
    output = array(strmid(string_input, string_position, number_digits), copy=0).astype(float32)
    
    return output


def get_ion_fraction(string_input, ionization_state):
##Find that float after the given text marker:

    n_params = 2
    def _ret():  return None
    
    number_digits = 7.
    string_position = 11 + 7. * (ionization_state - 1)
    
    substring_input = strmid(string_input, string_position, number_digits)
    
    if (bitwise_and((strpos(substring_input, "-") == -1), (strpos(substring_input, '0.00') == -1))):    
        output = -30.0
    else:    
        output = array(substring_input, copy=0).astype(float32)
    
    return output



######################################################################
##Main program
######################################################################

def cloudy_fits(output_filename, OUTPUT_STRUCTURE=None, INPUT_FILES=None, HELP=None):
    """
     NAME: CLOUDY_FITS
    
     PURPOSE:
    
     Read in a directory of CLOUDY output files and save the
     pertinent information into a FITS file.
    
    
     CALLING SEQUENCE:
    
        CLOUDY_FITS, output_filename, OUTPUT_STRUCTURE = cloudy, $
             INPUT_FILES = 'search_string'
    
     INPUTS:
    
       output_filename -- Filename for output FITS file.  Appends
                          '.fits'.  DEFAULT: cloudy_output.fits
    
    
     KEYWORD PARAMETERS:
    
       OUTPUT_STRUCT   -- Set this keyword to a variable to have that
                          variable filled with output structure.
    
       INPUT_FILES     -- Set this keyword to string to constrain files
                          read in by cloudy_fits.  All files assumed to
                          end in '.out' suffix.
    
     MODIFICATION HISTORY:
    
       03/15/02 -- JCH -- Created.
       03/21/02 -- JCH -- Updated for use with different versions of
                           Cloudy.
       03/25/02 -- JCH -- Changed search string for H II column to accommodate
                           multiple versions of the Cloudy output.
       03/27/02 -- JCH -- Added filename and caution tags to the output
                           structure.
                       -- Changed h1column,h2column,h12column tags to
                           hold logarithmic column densities.
                       -- Changed code involving strsplit/str_sep
                           routines to use the appropriate routine for
                           the version of IDL being used.
       04/10/02 -- JCH -- Added stop tag; collect stopping criterion,
                           and warnings as well as cautions.
       04/16/02 -- JCH -- Added log_jnu tag: value of J_nu.
    
       04/26/06 -- JCH -- Adapted to v06.02.
    
       05/18/07 -- JCH -- Included CI, CII, OI, and SiII fine structure columns.
                           New tags: c1star0,c1star1,c1star2,c2star,
                                     o1star0,o1star1,o1star2,si2star
    
       03/17/11 -- JCH -- Updated procedure to be compatible with
                           Cloudy v08.01.
    
    
    """

    n_params = 1
    cloudy = OUTPUT_STRUCTURE
    inputfilebase = INPUT_FILES
    help = HELP
    _opt = (cloudy, inputfilebase, help)
    def _ret():
        _optrv = zip(_opt, [cloudy, inputfilebase, help])
        _rv = [output_filename]
        _rv += [_o[1] for _o in _optrv if _o[0] is not None]
        return tuple(_rv)
    
    if (help is not None):    
        print ''
        print 'Calling Sequence: '
        print '        CLOUDY_FITS, output_filename, ' + 'OUTPUT_STRUCTURE = cloudy, $'
        print '                       INPUT_FILES = inputfilebase, /HELP'
        retall()
    
    
    if output_filename.size == 0:    
        output_filename = 'cloudy_output.fits'
    else:    
        if bitwise_or(strlen(output_filename) <= 5, strmid(output_filename, strlen(output_filename) - 5, 5) != '.fits'):    
            output_filename = output_filename + '.fits'
    
    ##Find CLOUDY output files:
    if bitwise_not((inputfilebase is not None)):    
        search_string = '*.out'
        
        ##Find the input Cloudy files
        spawn('ls ' + search_string, cloudy_files)
        
    else:    
        if inputfilebase.size == 0:    
            suffix = strmid(inputfilebase[0], 0, strlen(inputfilebase[0]) - 4)
            if suffix != '.out':    
                cloudy_files = inputfilebase + '.out'
            else:    
                cloudy_files = inputfilebase
        else:    
            spawn('ls ' + inputfilebase, cloudy_files)
    
    ##Number of Cloudy files to use:
    num_files = cloudy_files.size
    
    ########################################################################
    ##
    ## Structure definitions:
    
    
    ##Structure elem will hold abundances and ionization fractions:
    elem = element_info(h=zeros([3], "float32"), he=zeros([4], "float32"), c=zeros([8], "float32"), n=zeros([9], "float32"), o=zeros([10], "float32"), Neon=zeros([12], "float32"), mg=zeros([14], "float32"), al=zeros([15], "float32"), si=zeros([16], "float32"), p=zeros([17], "float32"), s=zeros([18], "float32"), cl=zeros([19], "float32"), ar=zeros([20], "float32"), ca=zeros([21], "float32"), ti=zeros([23], "float32"), cr=zeros([25], "float32"), mn=zeros([26], "float32"), fe=zeros([27], "float32"), ni=zeros([29], "float32"), zn=zeros([31], "float32"))
    
    ## Save the elemental tags (and number of them) out of the
    ## elem structure:
    elements = strmid(tag_names(elem), 0, 2)
    num_OF_elements = elements.size
    
    ## Structure excite will hold columns of excited states:
    excite = excitation_info(c1star0=-30.0, c1star1=-30.0, c1star2=-30.0, c2star=-30.0, o1star0=-30.0, o1star1=-30.0, o1star2=-30.0, si2star=-30.0)
    
    ## Save the excited state tags (and number of them)
    excited_states = tag_names(excite)
    num_OF_excited = excited_states.size
    
    
    ##Structure cloudy_results will hold results of CLOUDY models...
    ##  includes elem structure defined above.
    cloudy = cloudy(filename='', cloudy_version='', input_block=strarr(50), cautions=strarr(50), stop=strarr(5), exit=strarr(2), h12column=0., h1column=0., h2column=0., hhcolumn=0., log_u=0., log_jnu=0., temperature=0., density=0.)
    
    ##Now make an array with one structure for each CLOUDY model.
    cloudy = (cloudy)*ones([num_files])
    
    ##
    ######################################################################
    
    
    input_string = ''            # Initialize string to hold input.
    
    ##Loop through the CLOUDY output files.
    for ll in arange(0, (num_files - 1)+(1)):
    
    ######################################################################
    ##
    ##(0): Open file for reading, save filename, cloudy version into structure:
    ##
        cloudy[ll].filename = cloudy_files[ll]
        openr(input_unit, cloudy_files[ll], get_lun=True)
        
        ##Identify the version of CLOUDY from the first line.
        ##Separate input_string into an array:
        readf(input_unit, input_string)
        if array(_sys_version.release, copy=0).astype(float32) >= 5.3:    
            input_string_array = strsplit(input_string, ' ', extract=True)
        else:    
            input_string_array = str_sep(input_string, ' ')
        
        cloudy[ll].cloudy_version = input_string_array[1]
        
        ##Save a copy of this for later as a float!
        cloudy_version = array(cloudy[ll].cloudy_version, copy=0).astype(float32)
        
        
        ##
        ##(1): Find  value of J_nu for model:
        ##
        ##
        
        nu912 = 2.998e10 / (911.9e-8) #Frequency of Lyman limit
        
        text_flag = 'U(1.0----):'
        while strmid(strtrim(input_string, 2), 0, strlen(text_flag)) != text_flag:
            readf(input_unit, input_string)
        
        text_flag = 'nuJnu(912A):' #in ergs/cm^2/s
        nu_jnu = float_after_text(input_string, text_flag, 9)
        
        ##Convert value to log of  ergs/cm^2/s/Hz/Sr
        cloudy[ll].log_jnu = log10(nu_jnu / (4. * _sys_pi * nu912))
        
        
        ##
        ##(2): Get abundances of elements...
        ##
        
        ##Look for the 'Gas Phase Chemical Composition' text.
        text_flag = 'Gas Phase Chemical Composition'
        while strpos(input_string, text_flag) == -1:
            readf(input_unit, input_string)
        
        ##stop
        
        ##Elemental abundances are given over the next few lines:
        element_input = ''
        readf(input_unit, input_string)
        while strtrim(input_string, 2) != '':
            element_input = element_input + ' ' + strtrim(input_string, 2)
            readf(input_unit, input_string)
        
        if array(_sys_version.release, copy=0).astype(float32) >= 5.3:    
            element_input_array = strsplit(element_input, ' :', extract=True)
        else:    
            element_input_array = str_sep(strcompress(element_input), ':')
        
        element_input_name = element_input_array[0::2]
        element_input_abun = element_input_array[1::2]
        
        
        for j in arange(0, (element_input_name.size - 1)+(1)):
            neon = where(ravel(strlowcase(element_input_name) == 'ne'))[0]
            if neon != -1:    
                element_input_name[neon] = 'Neon'
            
            elem_index = where(ravel(strlowcase(element_input_name[j]) == strlowcase(tag_names(cloudy[ll]))))[0]
            elem_index = elem_index[0]
            
            
            if elem_index != -1:    
                cloudy[ll][elem_index][0] = element_input_abun[j]
        
        ##stop
        
        ##
        ##(3): Fill string array with information regarding the set-up of
        ##     the CLOUDY calculation.
        ##
        
        ##Look for the beginning of the info. block
        text_flag = '> Cloudy '
        while stregex(input_string, text_flag) == -1:
            readf(input_unit, input_string)
        
        ## OLD version:
        ##    WHILE strmid(strtrim(input_string, 1), 0, 40) NE text_flag DO $
        
        
        ##Fill cloudy.input_block array with information from CLOUDY file
        i = 0
        cloudy[ll].input_block[i] = strtrim(input_string, 2)
        
        ##Look for the end of the information block
        text_flag = '> Log\(U\)'
        while stregex(input_string, text_flag) == -1:
        ## OLD version:
        ##    WHILE strmid(strtrim(input_string, 1), 0, 41) NE  text_flag DO BEGIN
        
            readf(input_unit, input_string)
            i = i + 1
            cloudy[ll].input_block[i] = strtrim(input_string, 2)
            
        i = 0
        
        ##
        ##(4): Find ionization parameter
        ##
        text_flag = ' IONIZE PARMET:'
        while strmid(input_string, 0, strlen(text_flag)) != text_flag:
            readf(input_unit, input_string)
        
        input_string = strtrim(strcompress(input_string), 2)
        
        ##Separate input_string into an array:
        if array(_sys_version.release, copy=0).astype(float32) >= 5.3:    
            input_string_array = strsplit(input_string, ' ', extract=True)
        else:    
            input_string_array = str_sep(input_string, ' ')
        
        ion_param = where(ravel(input_string_array == "U(1-)"))[0]
        cloudy[ll].log_u = array(input_string_array[ion_param + 1], copy=0).astype(float32)
        
        ##This whole block could also be done:
        ##cloudy[ll].log_u = float_after_text(input_string, "U(1-)", 8)
        
        
        
        ##
        ##(5): Find average temperature, density
        ##
        
        text_flag = strcompress('Te      Te(Ne)   Te(NeNp)  Te(NeHe+)')
        while strmid(strcompress(strtrim(input_string, 2)), 0, strlen(text_flag)) != text_flag:
            readf(input_unit, input_string)
        
        ##  Next line is the interesting one...
        readf(input_unit, input_string)
        input_string = strcompress(strtrim(input_string, 2))
        
        ##Positional indexes of temperature and density values
        temp_index = 1
        dens_index = 9  ##<-- This changed from 8 in Cloudy96
        
        ##Extract entries in this line:
        if array(_sys_version.release, copy=0).astype(float32) >= 5.3:    
            input_string_array = strsplit(input_string, ' ', extract=True)
        else:    
            input_string_array = str_sep(strcompress(strtrim(input_string, 2)), ' ')
        
        ##Grab temperature and density:
        cloudy[ll].temperature = array(input_string_array[temp_index], copy=0).astype(float32)
        cloudy[ll].density = array(input_string_array[dens_index], copy=0).astype(float32)
        
        
        ##
        ##(6): Find hydrogen columns
        ##
        
        ##We have to treat Cloudy13 differently.
        
        ##Cloudy10 and earlier:
        if cloudy_version < 13:    
            
            ##
            ##(6a): Total hydrogen column:
            ##
            
            text_flag = 'Log10 Column density (cm^-2)'
            while strpos(input_string, text_flag) == -1:
                readf(input_unit, input_string)
            
            ##  Next line is the interesting one...
            readf(input_unit, input_string)
            input_string = strcompress(strtrim(input_string, 2))
            
            ##Positional indexes of temperature and density values
            htot_index = 2
            
            ##Extract entries in this line:
            if array(_sys_version.release, copy=0).astype(float32) >= 5.3:    
                input_string_array = strsplit(input_string, ' ', extract=True)
            else:    
                input_string_array = str_sep(strcompress(strtrim(input_string, 2)), ' ')
            
            ##Grab Total H column
            cloudy[ll].h12column = array(input_string_array[htot_index], copy=0).astype(float32)
            
            ##
            ##(6b): H I, H II, and H2 columns:
            ##
            ##text_flag='1 2 3 4 5 6 7 8 9'
            ##WHILE strmid(strtrim(strcompress(input_string),2), 0, strlen(text_flag)) $
            ## NE text_flag DO $
            ##   readf, input_unit, input_string
            
            text_flag = 'Log10 Column density (cm^-2)'
            while strpos(input_string, text_flag) == -1:
                readf(input_unit, input_string)
            
            ##Extract the column densities:
            if array(_sys_version.release, copy=0).astype(float32) >= 5.3:    
                input_string_array = strsplit(input_string, ' -', extract=True)
            else:    
                input_string_array = str_sep(strcompress(strtrim(input_string, 2)), ' ')
            
            ##Grab hydrogen columns:
            h1_index = 1
            h2_index = 2
            hh_index = 3
            
            cloudy[ll].h1column = input_string_array[h1_index]
            cloudy[ll].h2column = input_string_array[h2_index]
            cloudy[ll].hhcolumn = input_string_array[hh_index]
            
            ##Correct for non-existent H2:
            if cloudy[ll].hhcolumn == 30.:    
                cloudy[ll].hhcolumn = -30.00
            
            
        else:    
            ##Cloudy 13 version:
            
            ##
            ##(6a): H I, H II, and H2 columns:
            ##
            
            text_flag = 'Log10 Column density (cm^-2)'
            while strpos(input_string, text_flag) == -1:
                readf(input_unit, input_string)
            
            ##Extract the column densities:
            if array(_sys_version.release, copy=0).astype(float32) >= 5.3:    
                input_string_array = strsplit(input_string, ' -', extract=True)
            else:    
                input_string_array = str_sep(strcompress(strtrim(input_string, 2)), ' ')
            
            ##Grab hydrogen columns:
            h1_index = 1
            h2_index = 2
            hh_index = 3
            
            cloudy[ll].h1column = input_string_array[h1_index]
            cloudy[ll].h2column = input_string_array[h2_index]
            cloudy[ll].hhcolumn = input_string_array[hh_index]
            
            ##Correct for non-existent H2:
            if cloudy[ll].hhcolumn == 30.:    
                cloudy[ll].hhcolumn = -30.00
            
            ##
            ##(6b): Total hydrogen column:
            ##
            
            hcolumns = 10 ** array([[cloudy[ll].h1column, cloudy[ll].h2column]])
            
            cloudy[ll].h12column = log10(total(hcolumns))
            
            
        
        
        ##
        ##(7): Find excited state column densities
        ##
        
        ##First, check that they exist:
        spawn('grep "CII\*" ' + cloudy_files[ll], grepout)
        
        if grepout != '':    
            text_flag = 'Exc state'
            while strpos(input_string, text_flag) == -1:
                readf(input_unit, input_string)
            excitation_string = input_string
            
            ## Si II* is on the next line
            readf(input_unit, input_string)
            excitation_string = excitation_string + input_string ##<< Concatenate the 2 lines
            
            for j in arange(0, (num_OF_excited - 1)+(1)):
                exc_index = where(ravel(strupcase(tag_names(cloudy[ll])) == excited_states[j]))[0]
                
                _expr = strlowcase(excited_states[j])
                ##C I:
                if _expr == c1star0:    
                    cloudy[ll][exc_index] = float_after_text(excitation_string, 'C11*', 6)
                elif _expr == c1star1:    
                    cloudy[ll][exc_index] = float_after_text(excitation_string, 'C12*', 6)
                elif _expr == c1star2:    
                    cloudy[ll][exc_index] = float_after_text(excitation_string, 'C13*', 6)
                    
                    ##C II:
                elif _expr == c2star:    
                    cloudy[ll][exc_index] = float_after_text(excitation_string, 'CII*', 6)
                    
                    ##O I:
                elif _expr == o1star0:    
                    cloudy[ll][exc_index] = float_after_text(excitation_string, 'O11*', 6)
                elif _expr == o1star1:    
                    cloudy[ll][exc_index] = float_after_text(excitation_string, 'O12*', 6)
                elif _expr == o1star2:    
                    cloudy[ll][exc_index] = float_after_text(excitation_string, 'O13*', 6)
                    
                    ##Si II:
                elif _expr == si2star:    
                    cloudy[ll][exc_index] = float_after_text(excitation_string, 'Si2*', 6)
                    
                else:
                    raise RuntimeError('no match found for expression')
                
        
        ##
        ##(8): Grab ionization fractions for abundant elements (A=1 to 30)
        ##      [make sure to get averaged over radius...future will get
        ##      both?]
        
        ##Find the table of ionization fractions averaged over radius.
        text_flag = 'Log10 Mean Ionisation (over radius)'
        while strpos(input_string, text_flag) == -1:
            readf(input_unit, input_string)
        
        elem_text = strtrim(strmid(input_string, 1, 10), 2)
        
        while input_string != '':
            elem_text = strtrim(strmid(input_string, 1, 10), 2)
            _expr = strlowcase(elem_text)
            if _expr == hydrogen:    
                elem_abbrev = 'h'
            elif _expr == helium:    
                elem_abbrev = 'he'
            elif _expr == carbon:    
                elem_abbrev = 'c'
            elif _expr == nitrogen:    
                elem_abbrev = 'n'
            elif _expr == oxygen:    
                elem_abbrev = 'o'
            elif _expr == neon:    
                elem_abbrev = 'neon'
            elif _expr == sodium:    
                elem_abbrev = 'na'
            elif _expr == magnesium:    
                elem_abbrev = 'mg'
            elif _expr == aluminium:    
                elem_abbrev = 'al'
            elif _expr == silicon:    
                elem_abbrev = 'si'
            elif _expr == phosphorus:    
                elem_abbrev = 'p'
            elif _expr == sulphur:    
                elem_abbrev = 's'
            elif _expr == chlorine:    
                elem_abbrev = 'cl'
            elif _expr == argon:    
                elem_abbrev = 'ar'
            elif _expr == calcium:    
                elem_abbrev = 'ca'
            elif _expr == titanium:    
                elem_abbrev = 'ti'
            elif _expr == chromium:    
                elem_abbrev = 'cr'
            elif _expr == manganese:    
                elem_abbrev = 'mn'
            elif _expr == iron:    
                elem_abbrev = 'fe'
            elif _expr == cobalt:    
                elem_abbrev = 'co'
            elif _expr == nickel:    
                elem_abbrev = 'ni'
            elif _expr == copper:    
                elem_abbrev = 'cu'
            elif _expr == zinc:    
                elem_abbrev = 'zn'
                ##ELSE : GOTO, bad_element ; There's always a bad element.
            else:
                raise RuntimeError('no match found for expression')
            
            ##Find the tag name that goes with this element:
            elem_index = where(ravel(elem_abbrev == strlowcase(tag_names(cloudy[ll]))))[0]
            elem_index = elem_index[0]
            
            ##Make sure nothing happens if there's no corresponding
            ##elemental tag name.
            ##IF elem_index EQ -1 THEN GOTO, bad_element
            
            num_OF_electrons = cloudy[ll][elem_index].size
            
            for kk in arange(1, (num_OF_electrons - 1)+(1)):
            
                cloudy[ll][elem_index][kk] = get_ion_fraction(input_string, kk)
                
            
            ##print, elem_index, elem_text
            
            # bad_element:
            
            readf(input_unit, input_string)
        
        
        ##
        ## (9) Store stopping criterion, cautions, warnings, and such:
        ##
        
        
        spawn('egrep -- "Calculation stopped" ' + cloudy_files[ll], grep_stop)
        cloudy[ll].stop = grep_stop
        
        spawn('egrep -- "C-|W-|(!)" ' + cloudy_files[ll] + ' | egrep -v -- ">>>>>>>>>> Warning!"', grep_cautions)
        cloudy[ll].cautions = grep_cautions
        
        spawn('tail -2 ' + cloudy_files[ll], grep_exit)
        cloudy[ll].exit = grep_exit
        
        
        ##
        ######################################################################
        
        
        close(input_unit)
        free_lun(input_unit)
    
    
    ## --> This could be uncommented if an IDL SAVE file is preferred:
    ##save, file = 'cloudy.save', cloudy
    
    
    ######################################################################
    ##
    ## Output structure to FITS binary table.
    
    ##Write the structure to a FITS binary table:
    mwrfits(cloudy, output_filename, create=True)
    
    ## Let the user know what happened.
    print 'Wrote ' + output_filename
    
    ##
    ######################################################################
    
    
    
    
    return _ret()





