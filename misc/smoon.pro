;+
; Type: function.
;
; Purpose: Calc the Moon's altitude and azimuth for given time at glat/glon.
;
; Parameters:
;   et0, in, double, required. The epoch for given UTC.
;   glat, in, double, required. The geo lat in rad, in deg if set degree.
;   glon, in, double, required. The geo lon in rad, in deg if set degree.
;
; Keywords:
;   azimuth, out, double, optional. The azimuth in rad, in deg if set degree.
;   degree, in, boolean, optional. Set to return altitude and azimuth in deg.
;
; Return: doulbe. The moon's altitude in rad. In deg if set degree.
;
; Notes: none.
;
; Dependence: slib.
;
; History:
;   2015-07-02, Sheng Tian, create.
;-
function smoon, et0, glat0, glon0, azimuth = azm, degree = degree

    pi2 = 2*!dpi
    eps = 23.43929111d*!dtor    ; in rad.
    deg = 180d/!dpi
    rad = !dpi/180d
    arcs = 3600d*deg            ; arcs.
    
    mjd_t = sfmepoch(et0,'mjd')
    mjd_j2000 = 51544.5d    ; MJD J2000.
    
    t = mjd_t - mjd_j2000
    l0 = 0.606434d + 0.03660110129d*t   ; moon's mean longitude.
    l  = 0.374897d + 0.03629164709d*t   ; moon's mean anomaly.
    ls = 0.993133d + 0.0027377785d *t   ; sun's mean anomaly.
    d  = 0.827362d + 0.03386319198d*t   ; diff lon moon-sun.
    f  = 0.259091d + 0.03674819520d*t   ; dist from ascending node.
    
    l0 = l0 mod 1d
    l  = pi2*(l mod 1d)
    d  = pi2*(d mod 1d)
    f  = pi2*(f mod 1d)
    
    ; perturbations in lon and lat.
    dl = +22640d*sin(l) - 4586d*sin(l-2*d) + 2370d*sin(2*d) + 769d*sin(2*l) - $
        668d*sin(ls) - 412d*sin(2*f) - 212d*sin(2*l-2*d) - 206d*sin(l+ls-2*d) + $
        192d*sin(l+2*d) - 165d*sin(ls-2*d) - 125d*sin(d) - 110d*sin(l+ls) + $
        148d*sin(l-ls) - 55d*sin(2*f-2*d)
    s = f + (dl+412d*sin(2*f)+541d*sin(ls))/arcs
    h = f-2*d;
    n = -526d*sin(h) + 44d*sin(l+h) - 31d*sin(-l+h) - 23d*sin(ls+h) + $
        11d*sin(-ls+h) - 25d*sin(-2*l+f) + 21d*sin(-l+f);
        
    ; moon's ecliptic longitude and latitude in rad.
    lm = pi2*((l0+dl/1296d3) mod 1) ; lon.
    bm = (18520d*sin(s)+n)/arcs     ; iat.
    
    ; moon's ra and dec in rad.
    vec = transpose(cv_coord(from_sphere = [lm,bm,1d], /to_rect, /double))
    srotate, vec, eps, 0    ; rotate by -eps around x.
    ra = atan(vec[1],vec[0])
    dec = asin(vec[2])
    
    ; glat/glon.
    glat = glat0 & glon = glon0
    if keyword_set(degree) then begin
        glat = glat*rad & glon = glon*rad
    endif
    
    ; hour angle tau in rad, eqn (3.8).
    sgmst, et0, gmst, /radian       ; gmst in rad.
    tau = gmst + glon - ra
    
    ; from equatorial to horizontal coord.
    equ = transpose(cv_coord(from_sphere = [tau,dec,1d], /to_rect, /double))
    hor = equ & srotate, hor,-(!dpi*0.5-glat), 1
    azm = atan(hor[1],hor[0])       ; azimuth
    alt = asin(hor[2])              ; altitude in rad.

    if keyword_set(degree) then begin
        alt *= deg & azm *= deg
    end
    
    return, alt
    
end

ut0 = sfmdate('2015-07-01/12:50 CDT', '%Y-%m-%d/%H:%M:%S %Z')
et0 = stoepoch(ut0,'unix')
glat = 45*!dtor
glon =-93*!dtor
print, sfmepoch(et0)
print, smoon(et0, glat, glon, azimuth = azm, /degree), azm
end