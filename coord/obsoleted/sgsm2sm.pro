;+
; Type: function.
; Purpose: Convert vector in GSM coord to SM.
; Parameters:
;   vec0, in, dblarr[3]/dblarr[n,3], req. Vector(s) in GSM.
;   et, in, double/dblarr[n], req. The UT epoch(s).
; Keywords: none.
; Return: dblarr[3]/dblarr[3,n]. Vector(s) in SM. n is the num of records.
; Dependence: slib/coord.
; Notes: Keep the unit of input vector in Cartesian coordinate.
; History:
;   2012-07-13, Sheng Tian, create.
;-

function sgsm2sm, vec0, et
    compile_opt idl2 & on_error, 2
    
    vec1 = double(vec0)
    n1 = n_elements(vec1)/3 & n2 = n1+n1 & n3 = n2+n1
    vx0 = vec1[0:n1-1]
    vy0 = vec1[n1:n2-1]
    vz0 = vec1[n2:n3-1]
    
    ; get t4.
    sdipoledir, et, qgx, qgy, qgz, /interp   ; qg in geo.
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
    vx1 = cosm*vx0 - sinm*vz0
    vy1 = vy0
    vz1 = sinm*vx0 + cosm*vz0
    
    vec1[0:n1-1] = temporary(vx1)
    vec1[n1:n2-1] = temporary(vy1)
    vec1[n2:n3-1] = temporary(vz1)
    return, vec1

end

; t4 = <-m,y> 
;    = [[cosm, 0D, -sinm], $
;       [  0D, 1D,    0D], $
;       [sinm, 0D,  cosm]]
