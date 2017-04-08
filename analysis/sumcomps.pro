pro SumComps,Col,Err  
;############################################################
;
;	This procedure calculates the total column
;	density and errors for a series of components 
;	given in log space.
;       
;	Col = Array of input column density.
;	Err = Array of input errors.
;
;	11/10/97 -- Created by jch.
;############################################################

; IF n_params(0) EQ 0 THEN BEGIN
   
;;    print
;;    print
;;    print,  '  Calling Sequence:'
;;    print,  '      SumComps,Col,Err'
;;    print
;;    print, '          Col = Array of input column density.'
;;    print, '          Err = Array of input errors.'
;; retall   
;; END

   
   
sz1 = size(Col)

Sum    = 0
SumErr = 0

for i=0,sz1(1)-1 do begin

	Sum    = Sum + 10^(Col(i))
	SumErr = SumErr + (10^(Col(i))*(1-10D^(Err(i))))^2

end

SumErr = (SumErr)^(0.5)

LogSum = alog10(Sum)

LogErr = alog10(1+(SumErr/Sum))

;; LogSum = alog10(total(10.0D^Col))
;; LogErr = sqrt( total((10.0D^Col*2.3*Err)^2) )/total(10D^Col)/2.3


print,LogSum,LogErr

end
