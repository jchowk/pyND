PRO extinct_uv, wave, Rv
  ;IF n_params() EQ 0 THEN BEGIN
  ;  print, 'FUVextinct, wave, Rv'
  ;  retall
  ;ENDIF

   
   ;IF NOT keyword_set(Rv) THEN Rv = 3.1
   
   wave = wave*1.e-4
   x = 1./wave
   
   IF (x GT 3.3) AND (x LE 8) THEN BEGIN 
      
      IF (x LT 5.9) THEN BEGIN 
         fb = (fa = 0)
         
      ENDIF ELSE BEGIN 
         fa = -0.4473*(x-5.9)^2.-0.009779*(x-5.9)^3.
         
         fb = 0.2130*(x-5.9)^2.+0.1207*(x-5.9)^3
      ENDELSE 
      
      ax = 1.752-0.316*x-0.104/((x-4.67)^2.+0.341)+fa
      
      bx = -3.090+1.825*x+1.206/((x-4.62)^2.+0.263)+fb
      
   ENDIF ELSE IF (x GT 8) AND (x LT 11) THEN BEGIN 
      
      ax = -1.073-0.628*(x-8)+0.137*(x-8)^2.-0.070*(x-8)^3
      
      bx = 13.67+4.257*(x-8)-0.42*(x-8)^2.+0.374*(x-8)^3.
      
   ENDIF ELSE BEGIN 
      ax = 0.
      bx = 0.
      print, 'Wavelength not in the FUV'
   ENDELSE 
      
      awave2av = ax+bx/Rv
   
   print, 'A(wave)/A(V) = '+strtrim(string(awave2av), 2)
   print, '          (x = '+strtrim(string(x), 2)+' )'
   
END 
