;+
; Type:
;   function.
; 
; Name:
;   sgse2sm.
;
; Purpose:
;   Convert vector in GSE coord to SM.
;
; Parameters:
;   vec0, in, type = dblarr[3]/dblarr[n,3], required.
;       Vector(s) in GSE.
;     
;   et, in, type = double/dblarr[n], required.
;       The UT epoch(s).
;   
; Keywords:
;   none.
;   
; Return:
;   return, out, type = dblarr[3]/dblarr[3,n].
;       Vector(s) in SM. n is the num of records.
;   
; Example:
;   vec1 = sgse2sm(vec0, epoch).
;   
; Dependence:
;   sdipoledir.
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

function sgse2sm, vec0, et

    compile_opt idl2
    on_error, 2
    
    vec1 = double(vec0)
    n1 = n_elements(vec1)/3 & n2 = n1+n1 & n3 = n2+n1
    vx0 = vec1[0:n1-1]
    vy0 = vec1[n1:n2-1]
    vz0 = vec1[n2:n3-1]
    
    ; get t3.
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
    p = atan(-sinl*qex+cosl*(cose*qey+sine*qgz), -sine*qey+cose*qgz)
    sinp = sin(p)
    cosp = cos(p)
    m = asin(cosl*qex+sinl*(cose*qey+sine*qgz))
    sinm = sin(m)
    cosm = cos(m)
    
    ; vectorized, so should be faster than matrix ##.
    vx1 = vx0
    vy1 = cosp*vy0 - sinp*vz0
    vz1 = sinp*vy0 + cosp*vz0
    vx2 = cosm*vx1 - sinm*vz1
    vz2 = sinm*vx1 + cosm*vz1
    
    vec1[0:n1-1] = temporary(vx2)
    vec1[n1:n2-1] = temporary(vy1)
    vec1[n2:n3-1] = temporary(vz2)
    return, vec1

end

; t3 = <-p,x> 
;    = [[1D,   0D,    0D], $
;       [0D, cosp, -sinp], $
;       [0D, sinp,  cosp]]
