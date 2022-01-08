
pro test_arc_geometry

    site = 'atha'
    asc = sread_thm_asc(0, site, vars = ['glat','glon'])
    glat0 = asc.(0).glat
    glon0 = (asc.(0).glon+360) mod 360
    asc = sread_thm_asc(0, site, vars = ['glat','glon','elev','azim','alti'], type = 'asf')
    elevs = asc.(0).elev
    azims = asc.(0).azim
    altis = asc.(0).alti
    altid = 0
    glats = reform(asc.(0).glat[altid,*,*])
    glons = reform(asc.(0).glon[altid,*,*])
    alti0 = altis[altid]*1e-3
    re = !const.r_earth*1e-3
    
    npx = 256
    rad = !dpi/180
    
    glon1s = dblarr(npx,npx)
    glat1s = dblarr(npx,npx)
    
    angls = dblarr(npx,npx)
    
    et = stoepoch('2013-05-01/07:38')
    cdf_epoch, et, yr, mo, dy, hr, mi, sc, /breakdown_epoch
    geopack_recalc, yr, mo, dy, hr, mi, sc, /date
    for i = 0, npx-1 do begin
        for j = 0, npx-1 do begin
            telev = elevs[i,j]
            tazim = azims[i,j]*(-1)
            if finite(telev,/nan) then begin
                glon1s[i,j] = !values.d_nan
                glat1s[i,j] = !values.d_nan
                continue
            endif
            tmp = re/(re+alti0) & cosa = cos(telev*rad) & sina = sin(telev*rad)
            cosb = tmp*cosa^2 + sina*sqrt(1-tmp^2*cosa^2)
            tdist = sqrt(re^2+(re+alti0)^2-2*re*(re+alti0)*cosb)
            ; vector of aurora relative to camera in horizontal coord in rect.
            va_hor = cv_coord(from_sphere = [tazim,telev,tdist], /to_rect, /degree)
            ; transform it to geographic coord in rect.
            va_geo = shor2geo(va_hor, glat0, glon0, /degree)
            ; camera's pos in geographic coord in rect.
            vc_geo = cv_coord(from_sphere = [glon0,glat0,re], /to_rect, /degree)
            ; aurora pos in geographic coord in rect.
            va_rec = va_geo + vc_geo
            ; convert it to sphere.
            va_sph = cv_coord(from_rect = va_rec, /to_sphere, /degree)
            glon1s[i,j] = (va_sph[0]+360) mod 360
            glat1s[i,j] = va_sph[1]
            ; get magnetic field direction, in hor in rect.
            geopack_igrf_geo, va_sph[2]/re, 90-va_sph[1], va_sph[0]+360, /degree, $
                bu, bs, be
            vb_hor = [-bs,be,bu]
            angls[i,j] = sang(vb_hor, va_hor, /degree)
        endfor
    endfor
    tvscl, glon1s, 0, /nan
    tvscl, glat1s, 1, /nan
    tvscl, angls, 2, /nan
    tvscl, glons, 3, /nan
    tvscl, glats, 4, /nan
    print, minmax(glon1s), minmax(glons)
    print, minmax(glat1s), minmax(glats)
    stop
end