;+
; Function:
;   scotrant4
;
; Purpose:
;   I calculate the transformation matrix T4 in Hapgood 1992,
;   doi: 10.1016/0032-0633(92)90012-D.
;   T4 is used for SM = T4 * GSM.
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
;     The transformation matrix T4.
;   
; Example:
;   t4 = scotrant4(mjd).
;   it4 = scotrant4(mjd, /inverse).
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

function scotrant4, mjd, inverse = inverse

  compile_opt idl2
  
  on_error, 2
 
  qg = sdipoledir(mjd)
  
  it1 = scotrant1(mjd, /inverse)
  t2 = scotrant2(mjd)
  
  qe = qg # it1 # t2
  qe = t2 ## it1 ## qg
  
  ; mu is the dipole tilt angle.
  mu = asin(qe[0])
  
  sinm = sin(mu)
  cosm = cos(mu)

  ; t4 = <-mu, y>.
  t4 = [[ cosm, 0D, -sinm], $
        [   0D, 1D,    0D], $
        [ sinm, 0D,  cosm]]
      
  if n_elements(inverse) eq 0 then return, t4 $
  else return, transpose(t4)
  
end
