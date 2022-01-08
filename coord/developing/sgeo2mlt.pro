; from glat, glon to mlat, mlt.

; glat, glon in degree.

pro sgeo2mlt, et, glat, glon, mlat, mlon, mlt

    mjd = (1D/86400000D)*et-678941D
    ut = (mjd - floor(mjd))*24D

    ; degree to radian.
    dtor = !dpi/180D
    rtod = 180D/!dpi
    
    glat *= dtor
    glon *= dtor

    ; glat, glon to [x,y,z].
    tmp = 0.5D*!dpi - glat
    gx = sin(tmp)*cos(glon)
    gy = sin(tmp)*sin(glon)
    gz = cos(tmp)
    
    ; dipole angles (approximated), see Hapgood 1992,
    ; p = phi-90 = 78.8-90 + 4.283e-2*(mjd-46066)/365.25
    ; l = lambda = 289.1 - 1.413e-2*(mjd-46066)/365.25
    t0 = mjd - 46066D
    p = 0.0000020466107099D*t0 - 0.1954768229601431D    ; in radian.
    l = 5.0457468675156072D - 0.0000006751951087D*t0    ; in radian.
        
    ; geo to mag.
    cosp = sin(p) & sinp = cos(p)
    cosl = cos(l) & sinl = sin(l)
    ; rotation <l,z> -> [mx,my,gz].
    mx = cosl*gx + sinl*gy
    my = cosl*gy - sinl*gx
    
    ; rotation <p,y> -> [gx,my,mz].
    gx = sinp*mx+cosp*gz
    mz = sinp*gz-cosp*mx
    
    ; get mlat, mlon
    mlat = asin(mz)             ; in radian.
    mlon = atan(my,gx)          ; in radian.
    
    ; get mlt.
    ; solar angles, see sdipoledir.
    t0 = mjd - 51544.5D   ; should be mjd-mj2000, adjust from t0.
    ; mean obliquity of the ecliptic, e, in radian.
    e = 0.4090877233749509D - 6.2831853D-9*t0
    ; mean anomaly of the Sun, g, in radian.
    g = 6.2400582213628066D + 0.0172019699945780D*t0
    ; mean longitude of the Sun, q, in radian.
    q = 4.8949329668507771D + 0.0172027916955899D*t0
    ; geocentric apparent ecliptic longitude, l, in radian.
    l = q + 0.0334230551756914D*sin(g) + 0.0003490658503989D*sin(2*g)
    ; declination and right ascend.
    dec = asin(sin(e)*sin(l))   ; in radian.
    ra = atan(cos(e)*tan(l))    ; in radian.
    ; solar longitude.
    eqtime = q - ra - !dpi
    slon = !dpi *(1 - ut/12D) - eqtime
    dec = 0.5D*!dpi - dec
    gx = sin(dec)*cos(slon)
    gy = sin(dec)*sin(slon)
    gz = cos(dec)
    print, 'sgeo2mlt: ', dec*rtod, 0d, eqtime*rtod mod 360, ut
    print, 'xtt, slon geo: ', slon*rtod mod 360

    slon = atan(sinp*(cosl*gx + sinl*gy) + cosp*gz, $
        cosl*gy - sinl*gx)      ; in radian.
    print, 'xtt, slon mag: ', slon*rtod

    mlt = ((mlon-slon)*12/!dpi + 12) mod 24
    print, 'geo2mlt: mag', slon/!dtor mod 180
    
    mlat *= rtod
    mlon *= rtod
    
end