;+
; Function:
;   scotrant3
;
; Purpose:
;   I calculate the transformation matrix T3 in Hapgood 1992,
;   doi: 10.1016/0032-0633(92)90012-D.
;   T3 is used for GSM = T3 * GSE.
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
;     The transformation matrix T3.
;   
; Example:
;   t3 = scotrant3(mjd).
;   it3 = scotrant3(mjd, /inverse).
;   
; Dependence:
;   sdipoledir.
;   scotrant1.
;   scotrant2.
;   
; Notes:
;   none.
;  
; Author:
;   Sheng Tian
; 
; History:
;   2012-07-08, Sheng Tian, create.
;-

function scotrant3, mjd, inverse = inverse

  compile_opt idl2
  
  on_error, 2
  
  qg = sdipoledir(mjd)
  
  it1 = scotrant1(mjd, /inverse)
  t2 = scotrant2(mjd)
  
  qe = qg # it1 # t2
  
  ; psi = arctan(y_e/z_e)
  psi = atan(qe[1]/qe[2])
  
  sinp = sin(psi)
  cosp = cos(psi)
  
  ; t3 = <-psi, x>.
  t3 = [[1D,   0D,    0D], $
        [0D, cosp, -sinp], $
        [0D, sinp,  cosp]]
  
  if n_elements(inverse) eq 0 then return, t3 $
  else return, transpose(t3)

end
