
pro find_arc_height_for_2013_0501_event

    et = '2013-05-01/07:38:00'
    site1 = 'atha'
    site2 = 'tpas'
    
    minlat = 50
    minlat = 50
    xr = [-90,90]
    yr = [minlat,90]
    ytickpos = -45
    white = 255
    xtickname = ['18','00','06']
    top = 255
    ct = 39
    
    pos1 = [0.0,0.5,0.5,1.0]
    pos2 = [0.5,0.5,1.0,1.0]
    pos3 = [0.0,0.0,0.5,0.5]
    pos4 = [0.5,0.0,0.75,0.5]
    pos5 = [0.75,0.0,1.0,0.5]
    
    ; conversion from co-lat to elevation.
    rad = !dpi/180d
    deg = 180d/!dpi
    r = !const.r_earth*1e-3 ; km.
    h = 110d                ; km.
    asc = sread_thm_asc(0, site1, vars = ['glat','glon'], type = 'asc')
    glat0 = asc.atha.glat[0]    ; camera's glat in deg.
    glon0 = asc.atha.glon[0]    ; camera's glon in deg.
    asc = sread_thm_asc(0, site1, vars = ['glat','glon','elev'], type = 'asf')
    glats = reform(asc.atha.glat[0,1,*,*])      ; pixel glats in deg.
    glons = reform(asc.atha.glon[0,1,*,*])      ; pixel glats in deg.
    elevs = reform(asc.atha.elev[0,*,*])        ; pixel elevation in deg.
    idx = where(~finite(elevs,/nan), cnt)
    elev0s = dblarr(cnt)
    televs = dblarr(cnt)
    for i = 0, cnt-1 do begin
        tidx = array_indices(elevs, idx[i])
        tglat = mean(glats[tidx[0]:tidx[0]+1,tidx[1]:tidx[1]+1])
        tglon = mean(glons[tidx[0]:tidx[0]+1,tidx[1]:tidx[1]+1])
        elev0s[i] = elevs[tidx[0],tidx[1]]
        vec0 = cv_coord(from_sphere=[glon0,glat0,1], /to_rect, /degree)
        vec1 = cv_coord(from_sphere=[tglon,tglat,1], /to_rect, /degree)
        b = acos(sdot(vec0,vec1))
        u = (r+h)*cos(b)-r
        d = sqrt(r^2+(r+h)^2-2*r*(r+h)*cos(b))
        televs[i] = asin(u/d)*deg
    endfor
    plot, elev0s
    oplot, televs, color = sgcolor('red')
    stop
    
    ; plot mlt, asf for atha and tpas.
    atha = sread_thm_mlt(et, site1, /half, /dark, minlat = minlat)
    tpas = sread_thm_mlt(et, site2, /half, /dark, minlat = minlat)
    asi1 = sread_thm_mlt(et, [site1,site2], /half, /dark, minlat = minlat)
    asi2 = sread_thm_asi(et, site1, type = 'asf')
    asi2.img = bytscl(asi2.img, top = top, /nan)
    asi3 = sread_thm_asi(et, site1, type = 'ast')
    asi3.img = bytscl(asi3.img, top = top, /nan)

    
    fn = shomedir()+'/thm_asi_arc_2013_05_01.png'
    sgopen, fn, xsize = 16, ysize = 8, /inch
    sgindexcolor, ct, /silent
    sgtv, atha.mltimg, position = pos1
    sgset_map, position = pos1, xrange = xr, yrange = yr, $
        yminor = 9, $
        ytickpos = ytickpos, xtickname = xtickname, color = white
    sgtv, atha.mltimg, position = pos2
    sgset_map, position = pos2, xrange = xr, yrange = yr, $
        ytickpos = ytickpos, xtickname = xtickname, color = white
    sgtv, asi1.mltimg,  position = pos3
    sgset_map, position = pos3, xrange = xr, yrange = yr, $
        ytickpos = ytickpos, xtickname = xtickname, color = white
    sgtv, asi2.img, position = pos4
    sgtv, asi3.img, position = pos5
    sgclose

    ; plot asf for atha in a series of times.    
    tr = ['2013-05-01/07:35:00','2013-05-01/07:40']
    rootdir = shomedir()+'/thm_'+site1+'_2013_0501'
    asi = sread_thm_asi(tr, site1, type = 'asf')
    ets = stoepoch(asi.utsec, 'unix')
    nrec = n_elements(ets)
    for i = 0, nrec-1 do begin
        timg = double(reform(asi.img[i,*,*]))
        
        idx = where(finite(timg,/nan), cnt)
        if cnt ne 0 then timg[idx] = 0
        idx = where(timg gt 0)
        timg[idx] = alog(timg[idx])
        timg[idx] = (timg[idx]-alog(2000))*top/4
;        print, minmax(timg[idx]), median(timg[idx])

        sz = size(timg, /dimensions)
        fn = rootdir+'/thm_asf_'+site1+sfmepoch(ets[i],'_YYYY_MMDD_hhmm_ss.png')
        sgopen, fn, xsize = sz[0], ysize = sz[1]
        sgindexcolor, ct, /silent
        tv, timg
        xyouts, 0, 0, /device, sfmepoch(ets[i]), color = white
        sgclose
    endfor

end