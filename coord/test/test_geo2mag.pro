; (a) image/fuv's geotoapex.
; (b) sgeo2mag.
; (c) geo_to_mag.
; (d) aacgm.

pro test_geo2mag, et

    rtod = 180D/!dpi
    ssunlon, et, slon, slat, /degree
    ssunlon, et, mlon, /degree, /mag
    print, slon, mlon

    print, ''
    print, 'sun dir         glat (deg)      glon (deg)      mlat (deg)      mlon (deg)'

    ; (1) sgeo2mag.
    dir = sgeo2mag([slat,slon], et, /degree)
    mlat = dir[0]
    mlon = dir[1]
    print, 'sgeo2mag: ', slat, slon, mlat, mlon
    
    ; (2) geotoapex.
    apexfile = file_search( $
        '~/Dropbox/idl/idl82/Default/aurora/image/support/mlatlon.1997a.xdr')
    geotoapex, slat, slon, apexfile, tmp, mlon
    print, 'geo2apex: ', slat, slon, mlat, mlon
    
    ; (3) geo_to_mag.
;    mlat = slat & mlon = slon
;    mag_to_geo, mlat, mlon, /degree, /mag
;   print, 'magtogeo: ', slat, slon, mlat[0], mlon[0]-360
    
    ; (4) aacgm.
    aacgmidl
    cnv_aacgm, slat, slon, 1, mlat, mlon
    print, 'aacgmdef: ', slat, slon, double(mlat), double(mlon)

end


et = stoepoch('1998-06-17 01:06:30')

print, ''
print, '**time: ', sfmepoch(et)

test_geo2mag, et

end