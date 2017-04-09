;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;	CLOUDY_FITS.PRO:



FUNCTION float_after_text, string_input, text_value, number_digits
  ;;Find that float after the given text marker:
  IF n_elements(number_digits) EQ 0 THEN number_digits = 7
  flagpos = strpos(string_input, text_value)
  output = $
     float(strmid(string_input, $
                  flagpos+strlen(text_value), number_digits))  
  RETURN, output
END

FUNCTION float_at_position, string_input, string_position, number_digits
  ;;Find that float after the given text marker:
  IF n_elements(number_digits) EQ 0 THEN number_digits = 10

  output = $
     float(strmid(string_input, string_position, number_digits))  

  RETURN, output
END

FUNCTION get_ion_fraction, string_input, ionization_state
  ;;Find that float after the given text marker:

  number_digits = 7.
  string_position = 11+7.*(ionization_state-1)

  substring_input = strmid(string_input, string_position, number_digits)

  IF ((strpos(substring_input, "-") EQ -1) AND $
      (strpos(substring_input, '0.00') EQ -1) ) THEN output = -30.0 $
  else output = float(substring_input)

  RETURN, output
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Main program
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;+
; NAME: CLOUDY_COOLING
;
; PURPOSE:
;
; Read in a directory of CLOUDY output files and save the
; pertinent information into a FITS file.  This includes the 
; cooling information for various species.
;
; CALLING SEQUENCE:
;
;    CLOUDY_FITS, output_filename, OUTPUT_STRUCTURE = cloudy, $
;         INPUT_FILES = 'search_string'
;
; INPUTS:
;
;   output_filename -- Filename for output FITS file.  Appends
;                      '.fits'.  DEFAULT: cloudy_cooling.fits
;
;
; KEYWORD PARAMETERS:
;
;   OUTPUT_STRUCT   -- Set this keyword to a variable to have that
;                      variable filled with output structure.
;
;   INPUT_FILES     -- Set this keyword to string to constrain files
;                      read in by cloudy_fits.  All files assumed to 
;                      end in '.out' suffix.
;
; MODIFICATION HISTORY:
;
;   03/15/02 -- JCH -- Created.
;   03/21/02 -- JCH -- Updated for use with different versions of
;                       Cloudy.
;   03/25/02 -- JCH -- Changed search string for H II column to accommodate
;                       multiple versions of the Cloudy output.
;   03/27/02 -- JCH -- Added filename and caution tags to the output 
;                       structure. 
;                   -- Changed h1column,h2column,h12column tags to
;                       hold logarithmic column densities.
;                   -- Changed code involving strsplit/str_sep
;                       routines to use the appropriate routine for
;                       the version of IDL being used.
;   04/10/02 -- JCH -- Added stop tag; collect stopping criterion,
;                       and warnings as well as cautions.
;   
;   04/10/02 -- JCH -- Modified cloudy_fits routine to hold cooling info.
;   04/16/02 -- JCH -- Added log_jnu tag: value of J_nu.
;-

PRO cloudy_cooling, output_filename, $
                 OUTPUT_STRUCTURE = cloudy, INPUT_FILES = inputfilebase


  IF n_elements(output_filename) EQ 0 THEN $
     output_filename = 'cloudy_cooling.fits' $
  ELSE $
     IF strlen(output_filename) LE 5 OR $
     strmid(output_filename, strlen(output_filename)-5, 5) NE '.fits' $
     THEN $
     output_filename = output_filename+'.fits'

  ;;Find CLOUDY output files:
  IF NOT keyword_set(inputfilebase) THEN search_string = '*.out' $
     ELSE search_string = inputfilebase+'.out'

  spawn,'ls '+search_string,cloudy_files

  num_files = n_elements(cloudy_files) ;Number of Cloudy files in the directory

  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; Structure definitions:

  ;;Cooling information:
  cool = {cooling_info, heating: 0.0, cooling: 0.0, $
          c2cool: 0.0, n2cool: 0.0, $
          lyacool: 0.0, htotcool: 0.0, $
          bfh: 0.0, bfhe: 0.0, totm: 0.0, $
          grainheat: 0.0, graincool: 0.0}

  ;;Structure elem will hold abundances and ionization fractions:
  elem = {element_info, h:fltarr(3), he:fltarr(4), c:fltarr(8), $
          n:fltarr(9), o:fltarr(10), Neon:fltarr(12), mg:fltarr(14), $
          al:fltarr(15), si:fltarr(16), p:fltarr(17), s:fltarr(18), $
          ar:fltarr(19), ca:fltarr(21), ti:fltarr(23), cr:fltarr(25), $
          mn:fltarr(26), fe:fltarr(27), ni:fltarr(29), zn:fltarr(31)}

  ;; Save the elemental tags (and number of them) out of the
  ;; elem structure:
  elements = strmid(tag_names(elem), 0, 2)
  num_OF_elements = n_elements(elements) 
  

  ;;Structure cloudy_results will hold results of CLOUDY models...
  ;;  includes elem structure defined above.
  cloudy = $
     {cloudy_cool, filename:'', cloudy_version:'', $
      input_block: strarr(50), cautions: strarr(15), stop: strarr(5), $
      h12column:0., h1column:0., h2column:0., $
      log_u:0., log_jnu:0., temperature:0., density:0., $
      inherits element_info, inherits cooling_info} 

  ;;Now make an array with one structure for each CLOUDY model.
  cloudy = replicate(cloudy, num_files)                    

  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  input_string = ''             ; Initialize string to hold input.
  
  ;;Loop through the CLOUDY output files.
  FOR ll = 0, num_files-1 DO BEGIN 
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;
    ;;(0): Open file for reading and save filename into cloudy structure:
    ;;
    cloudy[ll].filename = cloudy_files[ll]
    openr, input_unit, cloudy_files[ll], /get_lun

    ;;
    ;;(1): Get abundances of elements...
    ;;

    ;;Look for the 'Chemical Composition' text.
    text_flag = 'Chemical Composition'
    WHILE strmid(strtrim(input_string, 2), 0, strlen(text_flag)) $
       NE text_flag DO $
       readf, input_unit, input_string

    ;; Loop through the elements in the list...
    FOR j = 0, num_OF_elements-1 DO BEGIN
      IF strlen(elements[j]) EQ 1 THEN $
         elem_tag = strupcase(elements[j]+' :') $
      ELSE elem_tag = strupcase(elements[j]+':')
      
      ;;Elements are arranged in three rows...damn!
      IF strpos(strupcase(input_string), elem_tag) EQ -1 THEN $
         readf, input_unit, input_string
      
      ;;Find tag index corresponding to current element:
      tag_position = $
         where(strmid(tag_names(cloudy[ll]), 0, 2) EQ $
               strupcase(elements[j]))
      tag_position = tag_position[0]
      
      
      IF tag_position EQ -1 THEN GOTO, bad_craziness $
      ELSE BEGIN 
        
        ;;Here we go!
        cloudy[ll].(tag_position)[0] = $
           float_after_text(strupcase(input_string), elem_tag, 8)

      ENDELSE 

      bad_craziness:

    ENDFOR


    ;;
    ;;(1.8): Find  value of J_nu for model:
    ;;
    ;;

    nu912 = 2.998e10/(911.9e-8) ;Frequency of Lyman limit
    
    text_flag = 'U(1.0----):'
    WHILE strmid(strtrim(input_string, 2), 0, strlen(text_flag)) $
       NE text_flag DO readf, input_unit, input_string
    
    text_flag = 'nuJnu(912A):'  ;in ergs/cm^2/s
    nu_jnu = $
       float_after_text(input_string, text_flag, 9)

    ;;Convert value to log of  ergs/cm^2/s/Hz/Sr
    cloudy[ll].log_jnu = alog10(nu_jnu/(4.*!pi*nu912))


    ;;
    ;;(2): Fill string array with information regarding the set-up of
    ;;     the CLOUDY calculation.
    ;;

    ;;Look for the beginning of the info. block
    text_flag = '*********************************> Cloudy' 
    WHILE strmid(strtrim(input_string, 1), 0, 41) NE text_flag DO $
       readf, input_unit, input_string


    ;;Identify the version of CLOUDY from the first line of the info
    ;; block
    cloudy[ll].cloudy_version = strmid(strtrim(input_string, 1), 43, 5)

    ;;Fill cloudy.input_block array with information from CLOUDY file
    i = 0
    cloudy[ll].input_block[i] = strtrim(input_string, 2)

    ;;Look for the end of the information block
    text_flag = '*********************************> Log(U)'
    WHILE strmid(strtrim(input_string, 1), 0, 41) NE  text_flag DO BEGIN
      
      readf, input_unit, input_string
      i = i+1
      cloudy[ll].input_block[i] = strtrim(input_string, 2)

    ENDWHILE 
    i = 0


    ;;
    ;;(2-cool): Cooling line intensities...
    ;;
    text_flag = 'Intrinsic line intensities'
    WHILE strmid(strtrim(input_string, 2), 0, strlen(text_flag)) NE $
       strtrim(text_flag, 2) DO $
       readf, input_unit, input_string 

    ;;Now find the various heating/cooling specifics:
    text_flag = 'BFH1    0'
    WHILE strmid(strtrim(input_string, 2), 0, strlen(text_flag)) NE text_flag DO $
       readf, input_unit, input_string     
    cloudy[ll].bfh = float_after_text(input_string, text_flag, 7)

    text_flag = 'BFHe    0'
    WHILE strmid(strtrim(input_string, 2), 0, strlen(text_flag)) NE text_flag DO $
       readf, input_unit, input_string     
    cloudy[ll].bfhe = float_after_text(input_string, text_flag, 7)

    text_flag = 'TotM    0'
    WHILE strmid(strtrim(input_string, 2), 0, strlen(text_flag)) NE text_flag DO $
       readf, input_unit, input_string     
    cloudy[ll].totm = float_after_text(input_string, text_flag, 7)

    text_flag = 'GrGH    0'
    WHILE strmid(strtrim(input_string, 2), 0, strlen(text_flag)) NE text_flag DO $
       readf, input_unit, input_string     
    cloudy[ll].grainheat = float_after_text(input_string, text_flag, 7)

    text_flag = 'GrGC    0'
    WHILE strmid(strtrim(input_string, 2), 0, strlen(text_flag)) NE text_flag DO $
       readf, input_unit, input_string     
    cloudy[ll].graincool = float_after_text(input_string, text_flag, 7)

        
    ;;Some you have to really search for...
    ;;condition = -1
    ;;text_flag = 'Clin  912'
    ;;WHILE condition EQ -1 DO BEGIN 
    ;;  readf, input_unit, input_string   
    ;;  condition = strpos(input_string, text_flag)
    ;;ENDWHILE
    ;;cloudy[ll].htotcool = float_after_text(input_string, text_flag, 7)
    ;;condition = -1

    ;;Better to do this with a grep command...
    text_flag = 'Clin  912'
    spawn, 'grep "'+text_flag+'" '+cloudy_files[ll], grep_output
    cloudy[ll].htotcool = float_after_text(grep_output[0], text_flag, 7)

    text_flag = "Cool 1216"
    spawn, 'grep "'+text_flag+'" '+cloudy_files[ll], grep_output
    cloudy[ll].lyacool = float_after_text(grep_output[0], text_flag, 7)

    text_flag = "C  2  157-|C  2  157 -|C  2  157  -"
    spawn, 'egrep "'+text_flag+'" '+cloudy_files[ll], grep_output
    text_flag = "C  2  157"
    cloudy[ll].c2cool = float_after_text(grep_output[0], text_flag , 7)

    text_flag = "N  2  205-|N  2  205 -|N  2  205  -"
    spawn, 'egrep "'+text_flag+'" '+cloudy_files[ll], grep_output
    text_flag = "N  2  205"
    cloudy[ll].n2cool = float_after_text(grep_output[0], text_flag , 7)

    ;;
    ;;(3): Find ionization parameter
    ;;
    text_flag = ' IONIZE PARMET:'
    WHILE strmid(input_string, 0, strlen(text_flag)) NE text_flag DO $
       readf, input_unit, input_string 
    
    input_string = strtrim(strcompress(input_string), 2)

    ;;Separate input_string into an array:
    IF float(!version.release) GE 5.3 THEN $
       input_string_array = strsplit(input_string, ' ', /extract) $
    ELSE input_string_array = str_sep(input_string, ' ')

    ion_param = where(input_string_array EQ "U(1-)")
    cloudy[ll].log_u = float(input_string_array[ion_param+1])

    ;;This whole block could also be done:
    ;;cloudy[ll].log_u = float_after_text(input_string, "U(1-)", 8)
   
   
    ;;
    ;;(3-cool): Find total heating/cooling:
    ;;
    text_flag = ' ENERGY BUDGET:'
    WHILE strmid(input_string, 0, strlen(text_flag)) NE text_flag DO $
       readf, input_unit, input_string 
    
    input_string = strtrim(strcompress(input_string), 2)

    ;;Separate input_string into an array:
    IF float(!version.release) GE 5.3 THEN $
       input_string_array = strsplit(input_string, ' ', /extract) $
    ELSE input_string_array = str_sep(input_string, ' ')

    heat_index = where(input_string_array EQ "Heat:")
    cloudy[ll].heating = float(input_string_array[heat_index+1])
    wind_chill = where(input_string_array EQ "Coolg:")
    cloudy[ll].cooling = float(input_string_array[wind_chill+1])

    ;;
    ;;(4): Find hydrogen columns
    ;;
    text_flag = strmid('Column density  H12:', 0, 15)
    WHILE strmid(strtrim(input_string, 2), 0, 15) NE text_flag DO $
       readf, input_unit, input_string 
    
    cloudy[ll].h12column = float_after_text(input_string, "H12:", 10)
    cloudy[ll].h1column  = float_after_text(input_string, "HI:", 10)
    cloudy[ll].h2column  = float_after_text(input_string, "II:", 10)

    ;;This could also be done such...
    ;;input_string = strtrim(strcompress(input_string), 2)
    ;;input_string_array = strsplit(input_string, ': ',/extract)
    ;;
    ;;h12index=where(input_string_array eq "H12")+1
    ;;h1index=where(input_string_array eq "HI")+1
    ;;h2index = $
    ;;   where((input_string_array eq "H II") or $
    ;;         (input_string_array EQ "HII"))+1
    ;;
    ;;cloudy[ll].h12column=float(input_string_array[h12index])
    ;;cloudy[ll].h1column=float(input_string_array[h1index])
    ;;cloudy[ll].h2column=float(input_string_array[h2index])


    ;;Store the results logarithmically
    cloudy[ll].h12column = alog10(cloudy[ll].h12column)
    cloudy[ll].h1column  = alog10(cloudy[ll].h1column)
    cloudy[ll].h2column  = alog10(cloudy[ll].h2column)


    ;;
    ;;(5): Find average temperature, density
    ;;

    text_flag = strcompress('Te       Te(Ne)    Te(NeNp)   Te(NeHe+)')
    WHILE strmid(strcompress(strtrim(input_string, 2)) $
                 , 0, strlen(text_flag)) $
       NE text_flag DO $
       readf, input_unit, input_string 

    ;;  Next line is the interesting one...
    readf, input_unit, input_string
    input_string = strcompress(strtrim(input_string, 2))
    
    ;;Positional indexes of temperature and density values
    temp_index = 1
    dens_index = 8

    ;;Extract entries in this line:
    IF float(!version.release) GE 5.3 THEN $
       input_string_array = strsplit(input_string,' ', /extract) $
    ELSE $
       input_string_array=str_sep(strcompress(strtrim(input_string,2)),' ')
    
    ;;Grab temperature and density:
    cloudy[ll].temperature = float(input_string_array[temp_index])
    cloudy[ll].density = float(input_string_array[dens_index])


    ;;
    ;;(6): Grab ionization fractions for abundant elements (A=1 to 30)
    ;;      [make sure to get averaged over radius...future will get
    ;;      both?]

    
    ;;First get there...find the next occurence of 'Hydrogen'
    text_flag = ' Hydrogen'
    WHILE strmid(input_string, 0, strlen(text_flag)) NE text_flag DO $
       readf, input_unit, input_string 

    ;;Check to see that this is averaged over radius:
    IF strpos(input_string, "Ionisation (over radius)") EQ -1 THEN BEGIN 
      readf, input_unit, input_string 
      WHILE strmid(input_string, 0, strlen(text_flag)) NE text_flag DO $
         readf, input_unit, input_string 
    ENDIF 

    elem_text = strtrim(strmid(input_string, 1, 10), 2)

    WHILE elem_text NE 'Zinc' DO BEGIN 
      elem_text = strtrim(strmid(input_string, 1, 10), 2)
      CASE strlowcase(elem_text) OF 
        'hydrogen'  : elem_abbrev = 'h'
        'helium'    : elem_abbrev = 'he' 
        'carbon'    : elem_abbrev = 'c' 
        'nitrogen'  : elem_abbrev = 'n' 
        'oxygen'    : elem_abbrev = 'o' 
        'neon'      : elem_abbrev = 'neon' 
        'sodium'    : elem_abbrev = 'na' 
        'magnesium' : elem_abbrev = 'mg' 
        'aluminium' : elem_abbrev = 'al' 
        'silicon'   : elem_abbrev = 'si' 
        'phosphorus': elem_abbrev = 'p' 
        'sulphur'   : elem_abbrev = 's' 
        'chlorine'  : elem_abbrev = 'cl' 
        'argon'     : elem_abbrev = 'ar' 
        'calcium'   : elem_abbrev = 'ca' 
        'titanium'  : elem_abbrev = 'ti' 
        'chromium'  : elem_abbrev = 'cr' 
        'manganese' : elem_abbrev = 'mn' 
        'iron'      : elem_abbrev = 'fe' 
        'cobalt'    : elem_abbrev = 'co' 
        'nickel'    : elem_abbrev = 'ni' 
        'copper'    : elem_abbrev = 'cu' 
        'zinc'      : elem_abbrev = 'zn' 
        ELSE : GOTO, bad_element ; There's always a bad element.
      ENDCASE
      
      ;;Find the tag name that goes with this element:
      elem_index = where(elem_abbrev EQ $
                         strlowcase(tag_names(cloudy[ll])))
      elem_index = elem_index[0]
      
      ;;Make sure nothing happens if there's no corresponding
      ;;elemental tag name.
      IF elem_index EQ -1 THEN GOTO, bad_element

      num_OF_electrons = n_elements(cloudy[ll].(elem_index))

      FOR kk = 1, num_OF_electrons-1 DO BEGIN         

        cloudy[ll].(elem_index)[kk] = $
           get_ion_fraction(input_string, kk)
        
      ENDFOR 

      ;;print, elem_index, elem_text
      
      bad_element:

      readf, input_unit, input_string 
    ENDWHILE 


    ;;
    ;; (7) Store stopping criterion, cautions, warnings, and such:
    ;;

    spawn, 'egrep -- "Calculation stopped" '+cloudy_files[ll], grep_stop
    cloudy[ll].stop = grep_stop

    spawn, 'egrep -- "C-|W-|(!)" '+cloudy_files[ll]+ $
           ' | egrep -v -- ">>>>>>>>>> Warning!"', grep_cautions
    cloudy[ll].cautions = grep_cautions

    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    close, input_unit
    free_lun, input_unit
  ENDFOR 


  ;; --> This could be uncommented if an IDL SAVE file is preferred:
  ;;save, file = 'cloudy.save', cloudy


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 
  ;; Output structure to FITS binary table.

  ;;Write the structure to a FITS binary table:
  mwrfits, /create, cloudy, output_filename

  ;; Let the user know what happened.
  print, 'Wrote '+output_filename

  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



END 




