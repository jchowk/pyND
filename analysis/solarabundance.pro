FUNCTION solarabundance, input, error=besterr, photo=PHOTO, meteor=METEOR
;;
;;
;; 3/7/2010 - Now uses abundances from: Asplund, M., Grevesse, N.,
;;            Sauval, A.J., & Scott, P. (2009), ARAA, 47, 481
;;
;; Calling sequence:
;; result = solarabundance(element_array, [error=err_out, /photo, /meteor])
;;
;;   


   IF n_params() EQ 0 THEN BEGIN 
      print, 'eps=solarabundance(input,error=err_variable,/photo,/meteor)'

      retall
   ENDIF 

   abundancefile = getenv('IDL_HOME')+ $
                   'idl/analysis/asplund2009_abundances.fits'

   a = mrdfits(abundancefile, 1, /silent)

   ;;Define output variables:
   bestabundance = fltarr(n_elements(input))
   besterr = bestabundance

   index = fltarr(n_elements(input))

   IF size(input, /type) EQ 7 THEN BEGIN 

      FOR j=0, n_elements(input)-1 DO $
         index[j] = $
         where(strtrim(a.element, 2) EQ strtrim(strlowcase(input[j]), 2)) 

   ENDIF ELSE BEGIN 

      FOR j=0, n_elements(input)-1 DO $
         index[j] = where(a.z EQ input[j])
      
   ENDELSE 

   FOR k=0, n_elements(index)-1 DO BEGIN 
      IF index[k] EQ -1 THEN BEGIN 

         bestabundance[k] = -30.
         besterr[k] = -30.

      ENDIF ELSE BEGIN 
         IF keyword_set(PHOTO) THEN BEGIN 

            bestabundance[k] = a[index[k]].photo
            besterr[k] = a[index[k]].photo_err
            
         ENDIF ELSE IF keyword_set(METEOR) THEN BEGIN 

            bestabundance[k] = a[index[k]].meteor
            besterr[k] = a[index[k]].meteor_err

         ENDIF ELSE BEGIN 

            bestabundance[k] = a[index[k]].best
            besterr[k] = a[index[k]].err
            
         ENDELSE 

      ENDELSE 
   ENDFOR


   bad = where(besterr EQ 30.)
   IF bad[0] NE -1 THEN BEGIN 
      bestabundance[bad] = -30.
      besterr[bad] = -30.
   ENDIF 
   
   return, bestabundance
END 
