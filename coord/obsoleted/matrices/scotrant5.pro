;+
; Function:
;   scotrant5
;
; Purpose:
;   I calculate the transformation matrix T5 in Hapgood 1992,
;   doi: 10.1016/0032-0633(92)90012-D.
;   T5 is used for MAG (GM) = T5 * GEO.
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
;     The transformation matrix T5.
;   
; Example:
;   t5 = scotrant5(mjd).
;   it5 = scotrant5(mjd, /inverse).
;   
; Dependence:
;   sdipoledir.
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

function scotrant5, mjd, inverse = inverse

  compile_opt idl2
  
  on_error, 2
  
  qg = sdipoledir(mjd)
  
  sinp = qg[2]
  cosp = -sqrt(1D - qg[2]^2)
  sinl = -qg[1] /cosp
  cosl = -qg[0] /cosp

  t5 = [[sinp *cosl, sinp *sinl, cosp], $
        [-sinl, cosl, 0], $
        [-cosp *cosl, -cosp *sinl, sinp]]

  if n_elements(inverse) eq 0 then return, t5 $
  else return, transpose(t5)
  
  ; t5 = <phi - 90, y> * <lambda, z>.
  ; or = <a, y> * <lambda, z>, a = phi - 90.
  ; cosa = sinp, sina = -cosp.
;  cosa = sinp
;  sina = -cosp
  
;  ry = [[cosa, 0, -sina], [0, 1, 0], [sina, 0, cosa]]
  
;  rz = [[cosl, sinl, 0], [-sinl, cosl, 0], [0, 0, 1]]
  
;  return, ry ## rz
  
end
