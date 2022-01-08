

pro calc_thg_asf_geometry, t0, site0

    ; const.
    re = !const.r_earth
    rad = !dpi/180d
    red = 254
    coef = 3    ; magnify.
    blue = 20
    ct = 39
    alti0 = 110000d    ; km.
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


    
    glons = dblarr(npx,npx)     ; pixel glon in deg.
    glats = dblarr(npx,npx)     ; pixel glat in deg.
    angls = dblarr(npx,npx)     ; look direction and aurora b field in deg.
    avecs = dblarr(npx,npx,3)   ; vector of each pixel in geo in xyz.
    bvecs = dblarr(npx,npx,3)   ; unit vec of B in geo in xyz.
    
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
            va_hor_xyz = cv_coord(from_sphere = [tazim,telev,tdist], /to_rect, /degree)
            ; transform it to geographic coord in rect.
            va_geo_xyz = shor2geo(va_hor_xyz, glat0, glon0, /degree)
            ; camera's pos in geographic coord in rect.
            vc_geo_xyz = cv_coord(from_sphere = [glon0,glat0,re], /to_rect, /degree)
            ; aurora pos in geographic coord in rect.
            va_geo_xyz = va_geo_xyz + vc_geo_xyz
            ; convert it to sphere.
            va_geo_sph = cv_coord(from_rect = va_geo_xyz, /to_sphere, /degree)
            glons[i,j] = (va_geo_sph[0]+360) mod 360
            glats[i,j] = va_geo_sph[1]
            ; get magnetic field direction, although the routine has geo in name,
            ; the returned b field is actually in hor in xyz. For example,
            ; the B at north pole is
            ; geopack_igrf_geo, 1, 0, 0, /degree, br, bt, bp
            ; br is the largest and positive, since B should be downward.
            ; the B at equator is
            ; geopack_igrf_geo, 1, 90, 0, /degree, br, bt, bp
            ; bt is the largest and negative, since B should be northward.
            ; br is radially out, bt (t for theta) is southward along
            ; constant meridian, bp (p for phi) is eastward.
            ; therefore in the last version, components are called bu, bs, be
            ; for up, southward and eastward.
            geopack_igrf_geo, va_geo_sph[2]/re, 90-va_geo_sph[1], va_geo_sph[0]+360, /degree, $
                bu, bs, be
            ; bu = br, bt = bs, bp = be
            vb_hor_xyz = [-bs,be,bu]
            
            ; angle b/w B and camara looking direction.
            angls[i,j] = sang(vb_hor_xyz, va_hor_xyz, /degree)
            bvecs[i,j,*] = sunitvec(vb_hor_xyz)
            avecs[i,j,*] = va_hor_xyz
        endfor
    endfor
    
    
    ; where look direct along b field.
    tmp = max(angls, idx, /nan) ; angle eq 180 means looking right oppose B.
    
    img[idx] = red
    ; draw image.
    isz = 600
    bgclr = 50
    
    ofn = shomedir()+'/arc_geom3.pdf'
;    ofn = 0
    sgopen, ofn, xsize = isz, ysize = isz
    tmp = findgen(21)/20*!dpi*2
    usersym, cos(tmp)*0.5, sin(tmp)*0.5, /fill
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    ; **** an imagined ray of aurora at given anchored at certain pixel.
    ;    stop
    ;    tis = [159,65]
    ;    tjs = [201,195]
    ;    cnts = [200,120]
    ;    nray = n_elements(tis)
    
    
    ; mask the image.
    idx = where(elevs ge -1)
    img[idx] = bgclr-10
    sgindexcolor, ct
    sgtv, img, position = [0,0,1,1]
    
    curtain0 = {glat:0d,glon:0d,dglat:0d,dglon:0d,bright:0d}
    curtains = []
    curtains = [curtains,{glat:52d,glon:260d,dglat:0.01,dglon:2.5,brightness:50}]
    curtains = [curtains,{glat:51d,glon:263d,dglat:0.01,dglon:2.5,brightness:80}]
    curtains = [curtains,{glat:50d,glon:266d,dglat:0.005,dglon:2.5,brightness:50}]
    ncurtain = n_elements(curtains)
    
    for k = 0, ncurtain-1 do begin
        tglat0 = curtains[k].glat
        tglon0 = curtains[k].glon
        dglat0 = curtains[k].dglat
        dglon0 = curtains[k].dglon
        bright = curtains[k].brightness
        rayidx = where(abs(glats-tglat0) le dglat0 and abs(glons-tglon0) le dglon0)
        nray = n_elements(rayidx)
        tmp = sort(glons[rayidx])
        rayidx = rayidx[tmp]
        
        tglats = glats & tglats[rayidx] = 255
        tglons = glons & tglons[rayidx] = 255
        
        
        for j = 0, nray-1 do begin
            tmp = array_indices(img, rayidx[j])
            ti = tmp[0]
            tj = tmp[1]
            tcnt = bright*sin(!dpi*j/nray)+150
            
            
            va_hor_xyz = reform(avecs[ti,tj,*])
            vb_hor_xyz = reform(bvecs[ti,tj,*])
            
            ; no need to figure out the geometry b/w B and each pixel's zenith
            ; this is b/c aurora ray is already along B field.
            dalti0 = 150*1e3
            daltis = smkarthm(0,dalti0,5e3,'dx')
            ndalti = n_elements(daltis)
            for i = 0, ndalti-1 do begin
                tdalti = daltis[i]
                tcolor = ((tcnt-50)*exp(-tdalti/dalti0*2.4)+bgclr)<254
                vc_hor_xyz = va_hor_xyz-vb_hor_xyz*tdalti
                vc_hor_sph = cv_coord(from_rec = vc_hor_xyz, /to_sphere, /degree)
                
                tazim = (-vc_hor_sph[0]+360) mod 360
                telev = vc_hor_sph[1]
                
                
                tmp = abs(elevs-telev)+abs(azims-tazim)
                idx = where(tmp eq min(tmp,/nan))
                tmp = array_indices(img, idx)
                pi = tmp[0]
                pj = tmp[1]
                plots, pi/256d, pj/256d, /normal, psym = 8, color = tcolor, symsize = (1-tdalti/dalti0*0.8)
            endfor
        endfor
        
        tmp = array_indices(img, rayidx[nray-1])
        ti = tmp[0]
        tj = tmp[1]        
        
        va_hor_xyz = reform(avecs[ti,tj,*])
        vc_hor_sph = cv_coord(from_rec = va_hor_xyz, /to_sphere, /degree)
        tazim = (-vc_hor_sph[0]+360) mod 360
        telev = vc_hor_sph[1]
        
        tmp = abs(elevs-telev)+abs(azims-tazim)
        idx = where(tmp eq min(tmp,/nan))
        tmp = array_indices(img, idx)
        pi = tmp[0]
        pj = tmp[1]
        
        xyouts, pi/256d, pj/256d, /normal, 'Curtain'+sgnum2str(k+1), color = 255
        xyouts, 0.01, 0.01+ychsz*1.2*(k+1), /normal, 'Curtain'+sgnum2str(k+1)+$
            ': height = '+string(alti0/1e3,format='(I0)')+'-'+string((alti0+dalti0)/1e3,format='(I0)')+' km'+$
            ', glat = '+string(tglat0,format='(F4.1)')+' deg'+$
            ', glon = '+string(tglon0-dglon0,format='(F5.1)')+'-'+string(tglon0+dglon0,format='(F5.1)')+' deg'+$
            ', max.cnt = '+string(bright,format='(I0)'), color = 255
    endfor
   
    
    ; the B zenith.
    tmp = max(angls, /nan, idx)
    tmp = array_indices(angls, idx)
    ti = tmp[0]
    tj = tmp[1]
    plots, ti/256d, tj/256d, /normal, psym = 6, color = red
    
    xyouts, 0.01, 0.01, /normal, time_string(ut0), color = 255
    sgclose
    
end

t0 = '2013-05-01/07:38:12'
site0 = 'atha'

t0 = '2013-06-07/04:56:06'
site0 = 'pina'
calc_thg_asf_geometry, t0, site0
end
