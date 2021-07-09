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

PRO cloudy_fits, output_filename, OUTPUT_STRUCTURE = cloudy, $
                 INPUT_FILES = inputfilebase, HELP = help
;+
; NAME: CLOUDY_FITS
;
; PURPOSE:
;
; Read in a directory of CLOUDY output files and save the
; pertinent information into a FITS file.
;
;
; CALLING SEQUENCE:
;
;    CLOUDY_FITS, output_filename, OUTPUT_STRUCTURE = cloudy, $
;         INPUT_FILES = 'search_string'
;
; INPUTS:
;
;   output_filename -- Filename for output FITS file.  Appends
;                      '.fits'.  DEFAULT: cloudy_output.fits
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
;   04/16/02 -- JCH -- Added log_jnu tag: value of J_nu.
;
;   04/26/06 -- JCH -- Adapted to v06.02.
;
;   05/18/07 -- JCH -- Included CI, CII, OI, and SiII fine structure columns.
;                       New tags: c1star0,c1star1,c1star2,c2star,
;                                 o1star0,o1star1,o1star2,si2star
;
;   03/17/11 -- JCH -- Updated procedure to be compatible with
;                       Cloudy v08.01.  
;
;
;-



   IF keyword_set(help) THEN BEGIN 
      print, ''
      print, 'Calling Sequence: '
      print, '        CLOUDY_FITS, output_filename, ' + $
             'OUTPUT_STRUCTURE = cloudy, $'
      print, '                       INPUT_FILES = inputfilebase, /HELP'
      retall
   ENDIF 


   IF n_elements(output_filename) EQ 0 THEN $
      output_filename = 'cloudy_output.fits' $
   ELSE $
      IF strlen(output_filename) LE 5 OR $
      strmid(output_filename, strlen(output_filename)-5, 5) NE '.fits' $
      THEN $
         output_filename = output_filename+'.fits'

   ;;Find CLOUDY output files:
   IF NOT keyword_set(inputfilebase) THEN BEGIN 
      search_string = '*.out' 

      ;;Find the input Cloudy files
      spawn,'ls '+search_string,cloudy_files

   ENDIF ELSE IF n_elements(inputfilebase) EQ 0 THEN BEGIN 
      suffix = strmid(inputfilebase[0], 0, strlen(inputfilebase[0])-4)
      IF suffix NE '.out' THEN $
         cloudy_files = inputfilebase+'.out' ELSE $
            cloudy_files = inputfilebase    
   ENDIF ELSE BEGIN 
      spawn, 'ls '+inputfilebase, cloudy_files
   ENDELSE 

   ;;Number of Cloudy files to use:
   num_files = n_elements(cloudy_files) 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;
   ;; Structure definitions:


   ;;Structure elem will hold abundances and ionization fractions:
   elem = {element_info, h:fltarr(3), he:fltarr(4), c:fltarr(8), $
           n:fltarr(9), o:fltarr(10), Neon:fltarr(12), mg:fltarr(14), $
           al:fltarr(15), si:fltarr(16), p:fltarr(17), s:fltarr(18), cl:fltarr(19), $
           ar:fltarr(20), ca:fltarr(21), ti:fltarr(23), cr:fltarr(25), $
           mn:fltarr(26), fe:fltarr(27), ni:fltarr(29), zn:fltarr(31)}

   ;; Save the elemental tags (and number of them) out of the
   ;; elem structure:
   elements = strmid(tag_names(elem), 0, 2)
   num_OF_elements = n_elements(elements) 

   ;; Structure excite will hold columns of excited states:
   excite =  {excitation_info, c1star0:-30.0, c1star1:-30.0, c1star2:-30.0, $
              c2star:-30.0, o1star0:-30.0, o1star1:-30.0, o1star2:-30.0, si2star:-30.0}

   ;; Save the excited state tags (and number of them) 
   excited_states = tag_names(excite)
   num_OF_excited = n_elements(excited_states) 
   

   ;;Structure cloudy_results will hold results of CLOUDY models...
   ;;  includes elem structure defined above.
   cloudy = $
      {cloudy, filename:'', cloudy_version:'', $
       input_block: strarr(50), cautions: strarr(50), $
       stop: strarr(5), exit: strarr(2), $
       h12column:0., h1column:0., h2column:0., hhcolumn:0.,$
       log_u:0., log_jnu:0., temperature:0., density:0., $
       inherits excitation_info, inherits element_info} 

   ;;Now make an array with one structure for each CLOUDY model.
   cloudy = replicate(cloudy, num_files)                    

   ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   input_string = ''            ; Initialize string to hold input.
   
   ;;Loop through the CLOUDY output files.
   FOR ll = 0, num_files-1 DO BEGIN 
      
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;
      ;;(0): Open file for reading, save filename, cloudy version into structure:
      ;;
      cloudy[ll].filename = cloudy_files[ll]
      openr, input_unit, cloudy_files[ll], /get_lun

      ;;Identify the version of CLOUDY from the first line.
      ;;Separate input_string into an array:
      readf, input_unit, input_string
      IF float(!version.release) GE 5.3 THEN $
         input_string_array = strsplit(input_string, ' ', /extract) $
      ELSE input_string_array = str_sep(input_string, ' ')

      cloudy[ll].cloudy_version = input_string_array[1]

      ;;Save a copy of this for later as a float!
      cloudy_version = float(cloudy[ll].cloudy_version)  


      ;;
      ;;(1): Find  value of J_nu for model:
      ;;
      ;;

      nu912 = 2.998e10/(911.9e-8) ;Frequency of Lyman limit
      
      text_flag = 'U(1.0----):'
      WHILE strmid(strtrim(input_string, 2), 0, strlen(text_flag)) $
         NE text_flag DO readf, input_unit, input_string
      
      text_flag = 'nuJnu(912A):' ;in ergs/cm^2/s
      nu_jnu = $
         float_after_text(input_string, text_flag, 9)

      ;;Convert value to log of  ergs/cm^2/s/Hz/Sr
      cloudy[ll].log_jnu = alog10(nu_jnu/(4.*!pi*nu912))


      ;;
      ;;(2): Get abundances of elements...
      ;;

      ;;Look for the 'Gas Phase Chemical Composition' text.
      text_flag = 'Gas Phase Chemical Composition'
      WHILE strpos(input_string,text_flag) eq -1 do $
         readf, input_unit, input_string

      ;;stop

      ;;Elemental abundances are given over the next few lines:
      element_input = ''
      readf,input_unit, input_string
      while strtrim(input_string,2) ne '' do begin 
         element_input = element_input + ' '+strtrim(input_string,2)
         readf,input_unit, input_string
      end 

      IF float(!version.release) GE 5.3 THEN $
         element_input_array = strsplit(element_input, ' :', /extract) $
      ELSE element_input_array = str_sep(strcompress(element_input), ':')

      element_input_name = element_input_array[0:*:2]
      element_input_abun = element_input_array[1:*:2]


      for j=0,n_elements(element_input_name)-1 do begin 
         neon=where(strlowcase(element_input_name) eq 'ne')
         if neon ne -1 then element_input_name[neon]='Neon'

         elem_index = where(strlowcase(element_input_name[j]) EQ $
                            strlowcase(tag_names(cloudy[ll])))
         elem_index = elem_index[0]
         

         if elem_index ne -1 then $
            cloudy[ll].(elem_index)[0] = element_input_abun[j]
      end 
      
      ;;stop

      ;;
      ;;(3): Fill string array with information regarding the set-up of
      ;;     the CLOUDY calculation.
      ;;

      ;;Look for the beginning of the info. block
      text_flag = '> Cloudy '
      WHILE stregex(input_string,text_flag) EQ -1 DO $
         readf, input_unit, input_string

      ;; OLD version:
      ;;    WHILE strmid(strtrim(input_string, 1), 0, 40) NE text_flag DO $


      ;;Fill cloudy.input_block array with information from CLOUDY file
      i = 0
      cloudy[ll].input_block[i] = strtrim(input_string, 2)

      ;;Look for the end of the information block
      text_flag = '> Log\(U\)'
      WHILE stregex(input_string, text_flag) EQ -1 DO BEGIN 
         ;; OLD version:
         ;;    WHILE strmid(strtrim(input_string, 1), 0, 41) NE  text_flag DO BEGIN
         
         readf, input_unit, input_string
         i = i+1
         cloudy[ll].input_block[i] = strtrim(input_string, 2)

      ENDWHILE 
      i = 0

      ;;
      ;;(4): Find ionization parameter
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
      ;;(5): Find average temperature, density
      ;;

      text_flag = strcompress('Te      Te(Ne)   Te(NeNp)  Te(NeHe+)')
      WHILE strmid(strcompress(strtrim(input_string, 2)) $
                   , 0, strlen(text_flag)) $
         NE text_flag DO $
            readf, input_unit, input_string 

      ;;  Next line is the interesting one...
      readf, input_unit, input_string
      input_string = strcompress(strtrim(input_string, 2))
      
      ;;Positional indexes of temperature and density values
      temp_index = 1
      dens_index = 9  ;;<-- This changed from 8 in Cloudy96

      ;;Extract entries in this line:
      IF float(!version.release) GE 5.3 THEN $
         input_string_array = strsplit(input_string,' ', /extract) $
      ELSE $
         input_string_array=str_sep(strcompress(strtrim(input_string,2)),' ')
      
      ;;Grab temperature and density:
      cloudy[ll].temperature = float(input_string_array[temp_index])
      cloudy[ll].density = float(input_string_array[dens_index])


      ;;
      ;;(6): Find hydrogen columns
      ;;

      ;;We have to treat Cloudy13 differently.

      ;;Cloudy10 and earlier:    
      IF cloudy_version LT 13 THEN BEGIN  

         ;;
         ;;(6a): Total hydrogen column:
         ;;

         text_flag = 'Log10 Column density (cm^-2)'
         WHILE strpos(input_string,text_flag) EQ -1 DO $
            readf, input_unit, input_string 

         ;;  Next line is the interesting one...
         readf, input_unit, input_string
         input_string = strcompress(strtrim(input_string, 2))
         
         ;;Positional indexes of temperature and density values
         htot_index = 2

         ;;Extract entries in this line:
         IF float(!version.release) GE 5.3 THEN $
            input_string_array = strsplit(input_string,' ', /extract) $
         ELSE $
            input_string_array=str_sep(strcompress(strtrim(input_string,2)),' ')
         
         ;;Grab Total H column
         cloudy[ll].h12column = float(input_string_array[htot_index])

         ;;
         ;;(6b): H I, H II, and H2 columns:
         ;;
         ;;text_flag='1 2 3 4 5 6 7 8 9'
         ;;WHILE strmid(strtrim(strcompress(input_string),2), 0, strlen(text_flag)) $
         ;; NE text_flag DO $
         ;;   readf, input_unit, input_string

         text_flag='Log10 Column density (cm^-2)'
         WHILE strpos(input_string, text_flag) EQ -1 DO $
            readf, input_unit, input_string

         ;;Extract the column densities:
         IF float(!version.release) GE 5.3 THEN $
            input_string_array = strsplit(input_string,' -', /extract) $
         ELSE $
            input_string_array=str_sep(strcompress(strtrim(input_string,2)),' ')

         ;;Grab hydrogen columns:
         h1_index = 1
         h2_index = 2
         hh_index = 3

         cloudy[ll].h1column  = input_string_array[h1_index]
         cloudy[ll].h2column  = input_string_array[h2_index]
         cloudy[ll].hhcolumn  = input_string_array[hh_index]

         ;;Correct for non-existent H2:
         IF cloudy[ll].hhcolumn EQ 30. THEN cloudy[ll].hhcolumn = -30.00


      ENDIF ELSE BEGIN 
         ;;Cloudy 13 version:

         ;;
         ;;(6a): H I, H II, and H2 columns:
         ;;

         text_flag='Log10 Column density (cm^-2)'
         WHILE strpos(input_string, text_flag) EQ -1 DO $
            readf, input_unit, input_string

         ;;Extract the column densities:
         IF float(!version.release) GE 5.3 THEN $
            input_string_array = strsplit(input_string,' -', /extract) $
         ELSE $
            input_string_array=str_sep(strcompress(strtrim(input_string,2)),' ')

         ;;Grab hydrogen columns:
         h1_index = 1
         h2_index = 2
         hh_index = 3

         cloudy[ll].h1column  = input_string_array[h1_index]
         cloudy[ll].h2column  = input_string_array[h2_index]
         cloudy[ll].hhcolumn  = input_string_array[hh_index]

         ;;Correct for non-existent H2:
         IF cloudy[ll].hhcolumn EQ 30. THEN cloudy[ll].hhcolumn = -30.00

         ;;
         ;;(6b): Total hydrogen column:
         ;;

         hcolumns = 10^[cloudy[ll].h1column, $
                        cloudy[ll].h2column]

         cloudy[ll].h12column = alog10(total(hcolumns))         


      ENDELSE 


      ;;
      ;;(7): Find excited state column densities
      ;;

      ;;First, check that they exist:
      spawn,'grep "CII\*" '+cloudy_files[ll],grepout
      
      IF grepout NE '' THEN BEGIN 
         text_flag = 'Exc state'
         WHILE strpos(input_string, text_flag) EQ -1 DO $
            readf, input_unit, input_string
         excitation_string = input_string
         
         ;; Si II* is on the next line
         readf, input_unit, input_string
         excitation_string = $
            excitation_string+input_string ;;<< Concatenate the 2 lines
         
         FOR j=0, num_OF_excited-1 DO BEGIN 
            exc_index = where(strupcase(tag_names(cloudy[ll])) $
                              EQ excited_states[j])
            
            CASE strlowcase(excited_states[j]) OF 
               ;;C I:
               'c1star0': cloudy[ll].(exc_index) = $
                  float_after_text(excitation_string,'C11*',6)
               'c1star1': cloudy[ll].(exc_index) = $
                  float_after_text(excitation_string,'C12*',6)
               'c1star2': cloudy[ll].(exc_index) = $
                  float_after_text(excitation_string,'C13*',6)

               ;;C II:
               'c2star' : cloudy[ll].(exc_index) = $
                  float_after_text(excitation_string,'CII*',6)

               ;;O I:
               'o1star0': cloudy[ll].(exc_index) = $
                  float_after_text(excitation_string,'O11*',6)
               'o1star1': cloudy[ll].(exc_index) = $
                  float_after_text(excitation_string,'O12*',6)
               'o1star2': cloudy[ll].(exc_index) = $
                  float_after_text(excitation_string,'O13*',6)

               ;;Si II:
               'si2star': cloudy[ll].(exc_index) = $
                  float_after_text(excitation_string,'Si2*',6)

            ENDCASE 
            
         ENDFOR
      ENDIF 

      ;;
      ;;(8): Grab ionization fractions for abundant elements (A=1 to 30)
      ;;      [make sure to get averaged over radius...future will get
      ;;      both?]

      ;;Find the table of ionization fractions averaged over radius.
      text_flag='Log10 Mean Ionisation (over radius)'
      WHILE strpos(input_string, text_flag) EQ -1 DO $
         readf, input_unit, input_string

      elem_text = strtrim(strmid(input_string, 1, 10), 2)

      WHILE input_string NE '' DO BEGIN 
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
      ;; (9) Store stopping criterion, cautions, warnings, and such:
      ;;

      
      spawn, 'egrep -- "Calculation stopped" '+cloudy_files[ll], grep_stop
      cloudy[ll].stop = grep_stop

      spawn, 'egrep -- "C-|W-|(!)" '+cloudy_files[ll]+ $
             ' | egrep -v -- ">>>>>>>>>> Warning!"', grep_cautions
      cloudy[ll].cautions = grep_cautions

      spawn, 'tail -2 '+cloudy_files[ll], grep_exit
      cloudy[ll].exit = grep_exit
      
      
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




