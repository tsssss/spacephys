
pro survey_polar_uvi_themis_asi, utr
    
    if n_elements(utr) eq 0 then utr = ['2008-02-12','2008-04-16']
;    if n_elements(utr) eq 0 then utr = ['2006-12-07','2008-04-16']
    ; '2006-01-01' to '2008-04-16'.
    tformat = (size(utr,/type) eq 7)? '': 'unix'
    etr = stoepoch(utr, tformat)
    
    ; settings.
    minlat = 50
    xr = [-90,90]
    yr = [minlat,90]
    ytickpos = -45
    white = 255
    xtickname = ['18','00','06']
    ct = 39 ; color table 39.
    xsz = 800
    ysz = 400
    uvipos = [0.1,0.1,0.9,0.9]
    uvitxtpos = [0.1,0.05]
    line1 = double(!d.y_ch_size)/ysz
    line2 = double(!d.y_ch_size)/228
    
    ; break down into days.
    etrs = sbreaktr(etr)
    utrs = sfmepoch(etrs,'unix')
    ntr = n_elements(utrs)/2
    for j = 0, ntr-1 do begin
        tr = reform(utrs[*,j])
        ; load ae and dst.
        omni = sread_omni(tr)
        ; load polar uvi data.
        uvi = sread_polar_uvi(tr, minlat = minlat, /half)
        if size(uvi,/type) ne 8 then continue   ; no polar data.
        rootdir = shomedir()+'/po_uvi_'+sfmepoch(uvi.epoch[0],'YYYY_MMDD')
        for i = 0, n_elements(uvi.epoch)-1 do begin
            tstr = sfmepoch(uvi.epoch[i],'YYYY_MMDD_hhmm_ss')+'.png'
            tae = interpol(omni.ae, omni.epoch, uvi.epoch[i])
            tdst = interpol(omni.symh, omni.epoch, uvi.epoch[i])
            ; mlt image.
            fn = rootdir+'/po_uvi_mlt_'+tstr
            sgopen, fn, xsize = xsz, ysize = ysz
            sgindexcolor, ct = 39, /silent
            tmp = bytscl(reform(uvi.mltimg[i,*,*]), top = 255)
            sgtv, tmp, position = uvipos
            sgset_map, xrange = xr, yrange = yr, pos = uvipos, $
                ytickpos = ytickpos, xtickname = xtickname, color = white
            xyouts, uvitxtpos[0], uvitxtpos[1], /normal, $
                'Polar UVI '+sfmepoch(uvi.epoch[i]), color = white
            xyouts, uvitxtpos[0], uvitxtpos[1]+1.5*line1, /normal, $
                'AE (nT) : '+string(tae,format='(I5)'), color = white
            xyouts, uvitxtpos[0], uvitxtpos[1]+3*line1, /normal, $
                'Dst (nT): '+string(tdst,format='(I5)'), color = white
            sgclose
            fn = rootdir+'/po_uvi_int_'+tstr
            sgopen, fn, xsize = 200, ysize = 228
            sgindexcolor, ct = 39, /silent
            tmp = bytscl(reform(uvi.img[i,*,*]), top = 255)
            sgtv, tmp
            xyouts, uvitxtpos[0], uvitxtpos[1], /normal, $
                'Polar UVI '+sfmepoch(uvi.epoch[i]), color = white
            xyouts, uvitxtpos[0], uvitxtpos[1]+1.5*line2, /normal, $
                'AE (nT) : '+string(tae,format='(I5)'), color = white
            xyouts, uvitxtpos[0], uvitxtpos[1]+3*line2, /normal, $
                'Dst (nT): '+string(tdst,format='(I5)'), color = white
            sgclose
        endfor
    endfor
end