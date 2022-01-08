;+
; Type:
;   function.
; 
; Name:
;   ssm2gsm.
;
; Purpose:
;   Convert vector in SM coord to GSM.
;
; Parameters:
;   vec0, in, type = dblarr[3]/dblarr[n,3], required.
;       Vector(s) in SM.
;     
;   et, in, type = double/dblarr[n], required.
;       The UT epoch(s).
;   
; Keywords:
;   none.
;   
; Return:
;   return, out, type = dblarr[3]/dblarr[3,n].
;       Vector(s) in GSM. n is the num of records.
;   
; Example:
;   vec1 = ssm2gsm(vec0, epoch).
;   
; Dependence:
;   sdipoledir.
;   ssundir.
;   sgmst.
;   
; Notes:
;   Keep the unit of input vector in Cartesian coordinate.
;  
; Author:
;   Sheng Tian.
; 
; History:
;   2012-07-13, Sheng Tian, create.
;-

function ssm2gsm, vec0, et

    compile_opt idl2
    on_error, 2
    
    vec1 = double(vec0)
    n1 = n_elements(vec1)/3 & n2 = n1+n1 & n3 = n2+n1
    vx0 = vec1[0:n1-1]
    vy0 = vec1[n1:n2-1]
    vz0 = vec1[n2:n3-1]
    
    ; get transpose(t4).
    sdipoledir, et, qgx, qgy, qgz   ; qg in geo.
    sgmst, et, gmst, /radian        ; geo to gei.
    sing = sin(gmst)
    cosg = cos(gmst)
    qex = cosg*qgx - sing*qgy
    qey = sing*qgx + cosg*qgy
    ssundir, et, e, l               ; gei to gse.
    sine = sin(e)
    cose = cos(e)
    sinl = sin(l)
    cosl = cos(l)
    m = asin(cosl*qex+sinl*(cose*qey+sine*qgz))
    sinm = sin(m)
    cosm = cos(m)
    
    ; vectorized, so should be faster than matrix ##.
    vx1 =  cosm*vx0 + sinm*vz0
    vy1 =  vy0
    vz1 = -sinm*vx0 + cosm*vz0
    
    vec1[0:n1-1] = temporary(vx1)
    vec1[n1:n2-1] = temporary(vy1)
    vec1[n2:n3-1] = temporary(vz1)
    return, vec1

end

; t4 = <-m,y> 
;    = [[cosm, 0D, -sinm], $
;       [  0D, 1D,    0D], $
;       [sinm, 0D,  cosm]]
