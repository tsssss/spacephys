;+
; Type: function.
;
; Purpose: Convert vector in HOR coord to GEO.
;
; Parameters:
;   vec0, in, dblarr[3]/dblarr[n,3], required. Vector(s) in HOR.
;   et, in, double/dblarr[n], required. The UT epoch(s).
;
; Keywords: none.
;
; Return: same dimension as vec0. Vector(s) in GEO. n is the num of records.
;
; Notes: Keep the unit of input vector in Cartesian coordinate.
;
; Dependence: ssundir.
;  
; History:
;   2015-07-09, Sheng Tian, create.
;-

function shor2geo, vec0, glat0, glon0, degree = degree

    compile_opt idl2
    on_error, 2

    vec1 = double(vec0)
    n1 = n_elements(vec1)/3 & n2 = n1+n1 & n3 = n2+n1
    vx0 = vec1[0:n1-1]
    vy0 = vec1[n1:n2-1]
    vz0 = vec1[n2:n3-1]

    ; get rotation matrix.
    rad = !dpi/180d
    glat = keyword_set(degree)? glat0*rad: glat0
    glon = keyword_set(degree)? glon0*rad: glon0
    
    srot, vx0,vy0,vz0,-!dpi, 'z'
    srot, vx0,vy0,vz0,-(!dpi*0.5-glat), 'y'
    srot, vx0,vy0,vz0,-glon, 'z'
    
    vec1[0:n1-1] = temporary(vx0)
    vec1[n1:n2-1] = temporary(vy0)
    vec1[n2:n3-1] = temporary(vz0)
    return, vec1
end

; l = glon, b = glat.
; M = | cosl  sinl  0 | * | sinb  0 -cosb | * |-1  0  0 |
;     |-sinl  cosl  0 |   |    0  1     0 |   | 0 -1  0 |
;     |    0     0  1 |   | cosb  0  sinb |   | 0  0  1 |
glat = 54.720d
glon = 246.69d
vec0 = cv_coord(from_sphere = [glon,glat+10,1], /degree, /to_rect)
vec1 = shor2geo(vec0, glat, glon, /degree)
vec1 = cv_coord(from_rec = vec1, /to_sphere, /degree)
vec0 = cv_coord(from_sphere = [glon,glat,1], /degree, /to_rect)
vec2 = shor2geo(vec0, glat, glon, /degree)
vec2 = cv_coord(from_rec = vec2, /to_sphere, /degree)
print, vec1
print, vec2
end
