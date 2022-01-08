;+
; Function:
;   scat2sph
;
; Purpose:
;   I convert 3-d vector from Cartesian to spherical coord.
;   [x, y, z] to [r, latitude, longitude].
;
; Parameters:
;   vec0, in, type = dblarr(3) or dblarr(3,n), required.
;     Vector in Cartesian coordinate.
;   
; Keywords:
;   degree = degree, in, type = boolean, optional.
;     Set degree to get latitude and longitude in degree.
;   
; Return:
;   return, out, type = dblarr(3) or dblarr(3,n).
;     The vector in spherical coordinate.
;   
; Example:
;   sph = scat2sph(cat, /degree).
;   sph = scat2sph(cat).
;   
; Dependence:
;   none.
;   
; Notes:
;   * Latitude and longitude are in radian by default.
;  
; Author:
;   Sheng Tian
; 
; History:
;   2012-07-17, Sheng Tian, create.
;-

function scat2sph, vec0, degree = degree
  
  compile_opt idl2
  
  on_error, 2
  
  angcoeff = n_elements(degree) ? 57.2957795130823230D :1D
  
    
  ; magnitude.
  rad = sqrt(vec0[0,*]^2 + vec0[1,*]^2 + vec0[2,*]^2)
  
  ; if magnitude is 0, then components in vec0 are 0, so is vec1.
  idx = where(rad ne 0)
  if idx[0] ge 0 then begin
    lon = atan(vec0[1,idx], vec0[0,idx])  ; longitude in radian.
    lat = asin(vec0[2,idx] /rad)   ; latitude in radian.
  endif
  
  return, [rad, lat*angcoeff, lon*angcoeff]
  
  ; L shell, r = L [cos(lambda)]^2, r is magnitude, 
  ; lambda is geomagnetic latitude.
  ; lshell = rad /cos(lat)^2
  
;  ; invariant latitude Gamma, L *[cos(Gamma)]^2 = 1
;  if n_elements(ilat) ne 0 then begin
;    ilat = acos(cos(lat) /sqrt(rad))
;    ; give ilat sign, sign is same as lat.
;    idx = where(lat lt 0)
;    if idx[0] ge 0 then ilat[idx] = -ilat[idx]   
;    return, [rad, ilat*angcoeff, lon*angcoeff]
;  endif else return, [rad, lat*angcoeff, lon*angcoeff]
  
end

print, scat2sph([[1,0,1],[0,1,1]], /degree)

end

