;+
; Function:
;   scotrant2
;
; Purpose:
;   I calculate the transformation matrix T2 in Hapgood 1992,
;   doi: 10.1016/0032-0633(92)90012-D.
;   T2 is used for GSE = T2 * GEI.
;
; Parameters:
;   mjd, in, type = double, required.
;     The modified Julian date.
;   
; Keywords:
;   inverse = inverse, in, type = boolean, optional.
;     Set inverse to get the inverse matrix.
;   
; Return:
;   return, out, type = dblarr(3,3).
;     The transformation matrix T2.
;   
; Example:
;   t2 = scotrant2(mjd).
;   it2 = scotrant2(mjd, /inverse).
;   
; Dependence:
;   none.
;   
; Notes:
;   * Coefficients are from:
;     http://aa.usno.navy.mil/faq/docs/SunApprox.php.
;     
;   * The angles (e, l, etc) has an accuracy of ~1 arcminute
;     within two centuries of 2000.
;  
; Author:
;   Sheng Tian
; 
; History:
;   2012-07-07, Sheng Tian, create.
;-

function scotrant2, mjd, inverse = inverse

  compile_opt idl2
  
  on_error, 2
  
  sundir = sgetsundir(mjd)
  
  sine = sin(sundir[0])
  cose = cos(sundir[0])
  sinl = sin(sundir[1])
  cosl = cos(sundir[1])
  
;  ra = atan(cose *sinl /cosl)
;  dec = asin(sine *sinl)
  
  ; t2 = <lambda, z> * <epsilon, x>.
  t2 = [[ cosl, sinl *cose, sinl *sine], $
        [-sinl, cosl *cose, cosl *sine], $
        [   0D,      -sine,       cose]]
  
  if n_elements(inverse) eq 0 then return, t2 $
  else return, transpose(t2)
  
end
  
  ; <lambda, z> = [[ cosl, sinl, 0D], $ 
  ;                [-sinl, cosl, 0D], $
  ;                [   0D,   0D, 1D]]
  ;         
  ; <epsilon, x> = [[1D,    0D,   0D], $
  ;                 [0D,  cose, sine], $
  ;                 [0D, -sine, cose]]
