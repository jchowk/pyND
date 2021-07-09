PRO mkcloudytable, aaaenergy, aaafnu, output_filename, $
   NUFNU=NUFNU, WAVE=WAVELENGTH, FREQUENCY=FREQUENCY
   
;---------------------------------------------------------------------
;  Routine to create file for use with Cloudy 'table read' commands.
;
;    
;
;  5/29/07 -- jch -- Created from mkhmrtable.pro
;---------------------------------------------------------------------   
   
;; Example Cloudy output and output from this code:
;;
;1234567890123456789012345678901234567890123
;1.052e-08       2.379e-11       1.000e+00   <-- Cloudy
;1.000e-06       1.000e-30       1.000e+00   <-- This code


   IF n_params() EQ 0 THEN BEGIN 

      print, 'mkcloudytable, energy[Ryd], fnu '
      print, '                 [, filename, /nufnu, /wave, /frequency]'
      print, '   ** Default filename is continuum_input.cloudy'

      retall 
   ENDIF 

   hplanck = 6.63e-27

   ;;Fill energy: input may be in Ryd, wavelength, or frequency
   IF keyword_set(wavelength) THEN $
      energy = 6.63e-27*2.998e10/(aaaenergy*1.e-8)/1.602e-12/13.606 $
   ;;(hplanck*2.998e10)/(aaaenergy*1.e-8)/(13.598*1.602e-12) $
   ELSE IF keyword_set(frequency) THEN  $
      energy = (6.63e-27*aaaenergy)/13.598/1.602e-12 $
   ELSE energy = aaaenergy

   ;;Normalize the fluxes
   IF keyword_set(NUFNU) THEN $
      fnu = aaafnu/energy/max(aaafnu) ELSE $
         fnu = aaafnu/max(aaafnu)


   ;;Make sure the energies are ascending
   ss = sort(energy)
   energy = energy[ss]
   fnu = fnu[ss]


   ;;Add low and high frequency points:
   energy = [1.e-6, energy, 7.e6]
   fnu = [1.e-30, fnu, 1.e-30]

   ;;Check for too low fluxes
   ugly = where(finite(fnu) EQ 0.)
   IF ugly[0] NE -1 THEN fnu[ugly] = 1.e-30


   ;;Open output file
   IF NOT keyword_set(output_filename) THEN $
      output_filename = 'continuum_input.cloudy'
   
   openw, unit, output_filename, /get_lun
   
   ;;Print a header...
   printf, unit, "#ener   Tran Contin     trn coef "

   FOR i=0, n_elements(energy)-1 DO BEGIN 

      printf, unit, energy[i], fnu[i], $
              format='(e10.3,6x,e10.3,"       1.000e+00")'

   ENDFOR 

   close, unit
   free_lun, unit
   
   print, 'Wrote '+output_filename

END 
