;+
; Type:
;   procedure.
;   
; Name:
;   ssunlon.
;
; Purpose:
;   Calculate the Sun's longitude, declination, right ascension
;   equation of time.
;   Based on http://aa.usno.navy.mil/faq/docs/SunApprox.php.
;
; Parameters:
;   et, in, type = double/dblarr[n], required.
;       Time in epoch UTC.
;       
;   slon, out, type = double/dblarr[n], required.
;       Longtitude of sun, default in GEO in radian.
;       Can be in MAG (GM), and/or in degree.
;     
;   dec, out, type = double/dblarr[n], optional.
;       Declination, default in radian. Also be the latitude in GEO.
;   
;   ra, out, type = double/dblarr[n], optional.
;       Right ascension, default in radian.
;
;   eqt, out, type = double/dblarr[n], optional.
;       Equation of time.
;
; Keywords:
;   mag = mag, in, type = boolean, optional.
;       Set to return geomagnetic longitude.
;       
;   degree = degree, in, type = boolean, optional.
;       Set to return angles in degree.
;
; Return:
;   none.
;
; Example:
;   ssunlon, et, glon, glat, /degree.
;
; Dependence:
;   none.
;
; Notes:
;   All angles are default in radian, without range limitation.
;   Contain modified segments from sdipoledir and ssundir.
;
; Author:
;   Sheng Tian.
;
; History:
;   2013-07-15, Sheng Tian, create.
;-

pro ssunlon, et, slon, dec, ra, eqt, mag = mag, degree = degree

    ut = (et/86400000D mod 1)*24
    
    ; from ssundir.
    d = (1D/86400000D)*et-730485.5D
    e = 0.4090877233749509D - 6.2831853D-9*d
    g = 6.2400582213628066D + 0.0172019699945780D*d
    q = 4.8949329668507771D + 0.0172027916955899D*d
    l = q + 0.0334230551756914D*sin(g) + 0.0003490658503989D*sin(2*g)
    
    ; declination and right ascend.
    dec = asin(sin(e)*sin(l))          ; in radian.
    ra = atan(cos(e)*sin(l),cos(l))    ; in radian. IMPORTANT! atan(y,x).
    
    ; solar longitude.
    eqt = q - ra                       ; in radian.
    slon = !dpi *(1 - ut/12D) - eqt    ; in radian, in geo.

    ; short version of geo to mag (gm).
    if keyword_set(mag) then begin
        vx0 = cos(dec)*cos(slon)
        vy0 = cos(dec)*sin(slon)
        vz0 = sin(dec)
        ; from sdipoledir.
        t0 = (1D/86400000D)*et - 725007D
        lat = 1.3753194505715316D + 0.0000020466107099D*t0    ; in radian.
        lon = 5.0457468675156072D - 0.0000006751951087D*t0    ; in radian.
        sinl = sin(lon)
        cosl = cos(lon)
        slon = atan(-sinl*vx0+cosl*vy0, sin(lat)*(cosl*vx0+sinl*vy0)-cos(lat)*vz0)
    endif
    
    if keyword_set(degree) then begin
        rtod = 180D/!dpi
        slon = (slon*rtod) mod 360
        dec *= rtod
        ra *= rtod
        eqt *= rtod
    endif

end