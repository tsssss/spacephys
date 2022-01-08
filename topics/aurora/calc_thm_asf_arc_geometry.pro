

pro calc_thm_asf_arc_geometry, t0, site0

    ; const.
    re = !const.r_earth
    rad = !dpi/180d
    red = 254
    coef = 3    ; magnify.
    blue = 20
    ct = 39
    alti0 = 110d    ; km.
    npx = 256
    
    ; read asf raw image.
    asi = sread_thg_asi(t0, site0[0], type = 'asf')
    img = float(asi.img)
    ut0 = asi.utsec
    et0 = stoepoch(ut0, 'unix')
    
    ; read calibration info.
    asc = sread_thg_asc(0, site0, vars = ['glat','glon'])
    glat0 = asc.(0).glat
    glon0 = (asc.(0).glon+360) mod 360
    asc = sread_thg_asc(0, site0, vars = ['elev','azim'], type = 'asf')
    elevs = asc.(0).elev
    azims = asc.(0).azim
    
    ; image processing.
    idx = where(img ne 0)
    minv = mean(img[0:10,0:10])
    maxv = 65535
    top = 254
    img[idx] = alog(img[idx]/minv)*(top/alog(maxv/minv))
    img = bytscl(img, min=40, max=150, top=top)
    edge = where(finite(elevs,/nan))
    img[edge] = !values.d_nan

     
    sgopen, shomedir()+'/arc_geomtry_eg.png', xsize = npx, ysize = npx
    sgindexcolor, 39
    sgtv, img
    xyouts, 1,1, /device, sfmepoch(et0)
    xyouts, 1,11,/device, strupcase(site0)
    sgclose
    
    glons = dblarr(npx,npx)     ; pixel glon in deg.
    glats = dblarr(npx,npx)     ; pixel glat in deg.
    angls = dblarr(npx,npx)     ; look direction and aurora b field in deg.
    
    cdf_epoch, et0, yr, mo, dy, hr, mi, sc, /breakdown_epoch
    geopack_recalc, yr, mo, dy, hr, mi, sc, /date
    for i = 0, npx-1 do begin
        for j = 0, npx-1 do begin
            telev = elevs[i,j]
            tazim = azims[i,j]*(-1)
            if finite(telev,/nan) then begin
                glons[i,j] = !values.d_nan
                glats[i,j] = !values.d_nan
                angls[i,j] = !values.d_nan
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
            glons[i,j] = (va_sph[0]+360) mod 360
            glats[i,j] = va_sph[1]
            ; get magnetic field direction, in hor in rect.
            geopack_igrf_geo, va_sph[2]/re, 90-va_sph[1], va_sph[0]+360, /degree, $
                bu, bs, be
            vb_hor = [-bs,be,bu]
            angls[i,j] = sang(vb_hor, va_hor, /degree)
        endfor
    endfor
    
    ; where look direct along b field.
    tmp = max(angls, idx, /nan)
    img[idx] = red
    ; draw image.
    isz = npx*coef
    timg = congrid(img, isz,isz)
    sgopen, 1, xsize = isz, ysize = isz
    sgindexcolor, ct
    sgtv, timg

    tmp = array_indices(angls, idx)
    
    stop    ; the following is wrong since the plot position should be involded.
    x0 = tmp[0]*coef & y0 = tmp[1]*coef
    tvcrs, x0, y0, /device
    button1 = 0     ; current button state.
    button0 = 0     ; previous button state.
    npts = 2
    xs = intarr(npts)  ; pixel x index.
    ys = intarr(npts)  ; pixel y index.
    i = 0
    while i lt npts do begin
        cursor, x1, y1, /change, /device
        button1 = !mouse.button
        if i eq 0 then begin
            sgtv, timg
            plots, [x0,x1], [y0,y1], /device, color = red
        endif
        case button1 of
            0: ; do nothing.
            1: begin
                xs[i] = x1 & ys[i] = y1
                i+= 1 & end
            2: ; do nothing?
            4: break
        endcase
    endwhile

    plots, xs, ys, /device, color = blue
    
    xs = xs/coef
    ys = ys/coef
    elevs = [elevs[xs[0],ys[0]],elevs[xs[1],ys[1]]]
    azims = [azims[xs[0],ys[0]],azims[xs[1],ys[1]]]
    vas = dblarr(2,3)
    for i = 0, 1 do begin
        telev = elevs[i]
        tazim = azims[i]
        tmp = re/(re+alti0) & cosa = cos(telev*rad) & sina = sin(telev*rad)
        cosb = tmp*cosa^2 + sina*sqrt(1-tmp^2*cosa^2)
        tdist = sqrt(re^2+(re+alti0)^2-2*re*(re+alti0)*cosb)
        ; vector of aurora relative to camera in horizontal coord in rect.
        va_hor = cv_coord(from_sphere = [tazim,telev,tdist], /to_rect, /degree)
        ; collect useful vars.
        vas[i,*] = va_hor
    endfor
    tmp = min(elevs, idx)
    val = reform(vas[idx,*])    ; arc's equatorward edge pos in hor coord.
    vah = reform(vas[1-idx,*])  ; poleward edge.
    
    ; get magnetic field direction at arc's equatorward edge, in hor in rect.
    ; transform arc pos to geographic coord in rect.
    val_geo = shor2geo(val, glat0, glon0, /degree)
    ; camera's pos in geographic coord in rect.
    vc_geo = cv_coord(from_sphere = [glon0,glat0,re], /to_rect, /degree)
    ; aurora pos in geographic coord in rect.
    val_geo = val_geo + vc_geo
    ; convert it to sphere.
    tmp = cv_coord(from_rect = val_geo, /to_sphere, /degree)
    ; get magnetic field direction, in hor in rect.
    geopack_igrf_geo, tmp[2]/re, 90-tmp[1], tmp[0]+360, /degree, bu, bs, be
    vb = [-bs,be,bu]
    
    ; do math.
    rl = snorm(val)
    rh = snorm(vah)
    a1 = sang(val, vah)
    w0 = sqrt(rl^2+rh^2-2*rl*rh*cos(a1))
    a2 = sang(val, vb)
    h0 = rl*sin(a1)/sin(!dpi-a1-a2)
    
    print, w0, h0
    
    img[xs[0],ys[0]] = red
    img[xs[1],ys[1]] = red
    sgopen, shomedir()+'/thm_arc.png', xsize = npx, ysize = npx
    sgindexcolor, ct
    sgtv, img
    xyouts, 1,1,/device, sfmepoch(et0)
    xyouts, 1,11,/device, strupcase(site0)
    xyouts, 1,21,/device, 'w0 = '+string(w0,format='(F7.2)')+' km', color = white
    xyouts, 1,31,/device, 'h0 = '+string(h0,format='(F7.2)')+' km', color = white
    sgclose
    
end

t0 = '2013-05-01/07:38:12'
site0 = 'atha'

t0 = '2013-06-07/04:56:06'
site0 = 'pina'
calc_thm_asf_arc_geometry, t0, site0
end