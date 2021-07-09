PRO mkhmrtable, data, redshift, output_filename
   
;---------------------------------------------------------------------
;  Routine to output Haardt, Madau, & Rees spectrum in a format
;  suitable for use with the CLOUDY table read command.
;   
;   ***Near jumps in the energy distribution, where the HMR spectrum
;      gives two flux values for the same energy, the higher energy
;      data point has been shifted to higher energy by 0.0001%.
;   
;   Output: file with CLOUDY interpolation table containing energy points
;           (in Rydbergs) and log F_nu points.
;
;   Created   6/10/99 by jch.   
;   Modified  3/21/02 by jch - Changed for use with data structures.
;---------------------------------------------------------------------   
   
  IF n_params() EQ 0 THEN BEGIN 

    print, 'mkcloudytable, data_structure, redshift[, filename]'
    print, '   ** Default filename is hmr_zZ.ZZZ.cloudy'

    retall 
  ENDIF 

  ;;Find nearest redshift:
  difference = abs(data.z-redshift)
  hold = min(difference, best)
  
  ;;Pull out the closest redshift data:
  z_best = data[best].z
  wave = data[best].wave
  fnu  = data[best].fnu
  energy = data[best].energy

  sort_list = sort(energy)
  wave = wave[sort_list]
  fnu  = fnu[sort_list]
  energy = energy[sort_list]
  
  log_fnu = alog10(fnu)


  ;;Derive output filename:
  IF NOT keyword_set(output_filename) THEN $
     output_filename = 'hmr_z'+ $
     strmid(strtrim(string(redshift), 2), 0, 5)+'.cloudy'
     
   openw, unit, output_filename, /get_lun
  
   ;;First output line needs "interpolate" command:
   printf, unit, float(0.), float(-35.), $
           format='("interpolate (",f11.6,1x,f8.3,")")'
   
   FOR i=0, n_elements(energy)-1 DO BEGIN 
     ;;If the adjacent energy entries are the same, then offset one
     ;; relative to the other...CLOUDY needs monotonically-increasing
     ;; data.
     IF i GT 0. THEN IF (energy[i] EQ energy[i-1]) THEN $
        energy[i] =  energy[i]*1.001

     ;;Output remaining values with "continue" command:
      printf, unit, energy[i], log_fnu[i], $
       format='("continue (",f11.6,1x,f8.3,")")'
   ENDFOR 
   printf, unit, float(9999.999), float(-35.), $
           format='("continue (",f11.6,1x,f8.3,")")'

   close, unit
   free_lun, unit
   
   print, 'Closest redshift: z = '+ $
          strmid(strtrim(string(z_best), 2), 0, 5)
   print, 'Wrote '+output_filename

END 
