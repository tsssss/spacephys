; test ways to get sun longitude, in geo.
; (a) image/fuv's get_local_time.
; (b) ssunlon.

pro test_sun_lon, et

    rtod = 180D/!dpi

    print, ''
    print, 'sun lon         ut (hr)         eqtime (deg)    dec (deg)       glon (deg)      mlon (deg)'

    ; (1) ssunlon.
    ut = (et/86400000D mod 1)*24
    ; declination and right ascend.
    ssundir, et, e, l, g, q
    dec = asin(sin(e)*sin(l))          ; in radian.
    ra = atan(cos(e)*sin(l),cos(l))    ; in radian.
    
    ; solar longitude.
    eqtime = q - ra
    slon = !dpi *(1 - ut/12D) - eqtime      ; in radian, in geo.
    ssunlon, et, slon, dec, ra, eqtime
    dir = sgeo2mag([dec, slon], et, /radian)
    
    print, 'ssunlon: ', ut, eqtime*rtod mod 360, dec*rtod, slon*rtod mod 360, dir[1]*rtod mod 360
    
    ; (2) from get_local_time.
    cdf_epoch, et, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    ut = hr + (mi + (sc + 1D-3*msc)/60D)/60D
    jd = julday(mo, dy, yr, 0)
    ephem, jd, ut, gha, dec, eqtime
    slon = (12D - ut - eqtime/15D)*15D
    apexfile = file_search('~/Dropbox/idl/idl82/Default/aurora/image/support/mlatlon.1997a.xdr')
    geo2apex, dec, slon, tmp, mlon
    print, 'fuv glt: ', ut, eqtime, dec, slon, mlon

end


et = stoepoch('1998-01-01 01:06:30')
et = stoepoch('2000-12-03 15:10:43.277','yyyy-mm-dd hh:mi:ss.msc')

print, ''
print, '**time: ', sfmepoch(et)

test_sun_lon, et

end