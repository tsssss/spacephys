;+
; Type:
;   procedure.
;
; Name:
;   ssundir.
; 
; Purpose:
;   Calculate angles for GSE <-> GSM.
;   Based on http://aa.usno.navy.mil/faq/docs/SunApprox.php.
;
; Parameters:
;   et, in, type = double/dblarr[n], required.
;     Time in epoch UTC.
;     
;   e, out, type = double/dblarr[n], required.
;     Mean obliquity of the ecliptic.
;   
;   l, out, type = double/dblarr[n], required.
;     Geocentric apparent ecliptic longitude.
;   
;   g, out, type = double/dblarr[n], optional.
;     Mean anomaly of the Sun.  
;
;   q, out, type = double/dblarr[n], optional.
;     Mean longitude of the Sun.
;
; Keywords:
;   degree, in, type = boolean, optional.
;     Set to return angles in degree.
;
; Return:
;   none.
;
; Example:
;   ssundir, t0, e, l, g, q, /degree.
;
; Dependence:
;   none.
;
; Notes:
;   All angles are default in radian, without range limitation.
;
; Author:
;   Sheng Tian.
;
; History:
;   2013-07-15, Sheng Tian, document.
;-

pro ssundir, et, e, l, g, q, degree = degree

    ; d = mjd - mj2000
    d = (1D/86400000D)*et-730485.5D
    
    ; mean obliquity of the ecliptic, e.
    ; e = 23.439 - 0.00000036*d     ; in degree.
    e = 0.4090877233749509D - 6.2831853D-9*d
    
    ; mean anomaly of the Sun, g.
    ; g = 357.529 + 0.98560028*d    ; in degree.
    g = 6.2400582213628066D + 0.0172019699945780D*d
    
    ; mean longitude of the Sun, q.
    ; q = 280.459 + 0.98564736*d   ; in degree.
    q = 4.8949329668507771D + 0.0172027916955899D*d
    
    ; geocentric apparent ecliptic longitude, l.
    ; l = q + 1.915 sin g + 0.020 sin 2g    ; in degree.
    l = q + 0.0334230551756914D*sin(g) + 0.0003490658503989D*sin(2*g)
    
    if keyword_set(degree) then begin
        rtod = 180D/!dpi
        e *= rtod
        g *= rtod
        q *= rtod
        l *= rtod
    endif
    
end