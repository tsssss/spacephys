;+
; Type:
;   function.
; 
; Name:
;   sgeo2mag.
;
; Purpose:
;   Convert vector in GEO coord to MAG (GM).
;
; Parameters:
;   vec0, in, type = dblarr[3]/dblarr[n,3]/dblarr[2]/dblarr[n,2], required.
;       Vector(s) in GEO. n is the num of records.
;     
;   et, in, type = double/dblarr[n], required.
;       The UT epoch(s).
;   
; Keywords:
;   degree = degree, in, type = boolean, optional.
;       Set to indicate vec0 and return value are in 
;       [lat, lon] in degree.
;
;   radian = radian, in, type = boolean, optional.
;       Set to indicate vec0 and return value are in 
;       [lat, lon] in radian.
;   
; Return:
;   return, out, type = dblarr[3]/dblarr[3,n]/dblarr[2]/dblarr[n,2].
;       Vector(s) in MAG (GM).
;   
; Example:
;   vec1 = sgeo2mag(vec0, epoch).
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

function sgeo2mag, vec0, et, degree = degree, radian = radian

    compile_opt idl2
    on_error, 2
    
    vec1 = double(vec0)
    
    if keyword_set(degree) or keyword_set(radian) then begin
        n1 = n_elements(vec1)/2 & n2 = n1+n1 & n3 = n2+n1
        dtor = !dpi/180D
        rtod = 180D/!dpi
        if keyword_set(degree) then begin
            lat = vec1[0:n1-1]*dtor
            lon = vec1[n1:n2-1]*dtor
        endif else begin
            lat = vec1[0:n1-1]
            lon = vec1[n1:n2-1]
        endelse
        vx0 = cos(lat)*cos(lon)
        vy0 = cos(lat)*sin(lon)
        vz0 = sin(lat)
    endif else begin
        n1 = n_elements(vec1)/3 & n2 = n1+n1 & n3 = n2+n1
        vx0 = vec1[0:n1-1]
        vy0 = vec1[n1:n2-1]
        vz0 = vec1[n2:n3-1]
    endelse
    
    ; get T5.
    sdipoledir, et, lat, lon, /radian
    sinp =  sin(lat)
    cosp = -cos(lat)
    sinl =  sin(lon)
    cosl =  cos(lon)
    
    ; vectorized, so should be faster than matrix ##.
    tmp =  cosl*vx0 + sinl*vy0
    vx1 =  sinp*tmp + cosp*vz0
    vy1 = -sinl*vx0 + cosl*vy0
    vz1 = -cosp*tmp + sinp*vz0
    
    if keyword_set(degree) or keyword_set(radian) then begin
        lat = asin(vz1)
        lon = atan(vy1, vx1)
        if keyword_set(degree) then begin
            lat *= rtod
            lon *= rtod
        endif
        vec1[0:n1-1] = lat
        vec1[n1:n2-1] = lon
        return, vec1
    endif
    
    vec1[0:n1-1] = temporary(vx1)
    vec1[n1:n2-1] = temporary(vy1)
    vec1[n2:n3-1] = temporary(vz1)
    return, vec1

end

;  t5 = [[ sinp*cosl,  sinp*sinl, cosp], $
;        [     -sinl,       cosl,    0], $
;        [-cosp*cosl, -cosp*sinl, sinp]]
