;+
; Function:
;   scotrant1
;
; Purpose:
;   I calculate the transformation matrix T1 in Hapgood 1992,
;   doi: 10.1016/0032-0633(92)90012-D.
;   T1 is used for GEO = T1 * GEI.
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
;     The transformation matrix T1.
;   
; Example:
;   t1 = scotrant1(mjd).
;   it1 = scotrant1(mjd, /inverse).
;   
; Dependence:
;   sgetgmst.
;   
; Notes:
;   none.
;  
; Author:
;   Sheng Tian
; 
; History:
;   2012-07-07, Sheng Tian, create.
;   2012-07-08, Sheng Tian, add inverse.
;-

function scotrant1, mjd, inverse = inverse

  compile_opt idl2
  
  on_error, 2
  
  ; get the GMST (Greenwich Mean Sidereal Time) in radian.
  theta = sgetgmst(mjd)
  
  cost = cos(theta)
  sint = sin(theta)

  ; t1 = <theta, z>
  t1 = [[ cost, sint, 0D], $
        [-sint, cost, 0D], $
        [   0D,   0D, 1D]]
        
  if n_elements(inverse) eq 0 then return, t1 $
  else return, transpose(t1)

end
