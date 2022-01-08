
_2013_0607_load_data


; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1


; **** settings.
    top = 254
    maxcnt = 800
    maxcnt = 500
    mincnt = 50
    dr0 = 3d
    
    tmodel = 't01'
    probes = ['a','b']
    
    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 5


; **** load data.
    get_data, 'asf_info', tmp, info
    mlats = info.mlats
    mlts = info.mlts
    imgsz = info.imgsz
    minlat = info.minlat
    
    ; calc image coord.
    txs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & txs = txs-imgsz[0]/2
    tys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & tys = tys-imgsz[0]/2
    txs = txs/imgsz[0]*2
    tys = tys/imgsz[0]*2
    mlats = 90-sqrt(txs^2+tys^2)*(90-minlat)
    mlts = atan(tys, txs)*deg+90    ; in deg, 0 deg at midnight.

    rootdir = shomedir()+'/psbl_de_32hz'
    scasifn = rootdir+'/2013_0607_scasi_'+tmodel+'.tplot'
    tmp = ['','_mlat','_mlt']
    scasivars = tmodel+['_asf_equatorial'+tmp,'_asf_vertical'+tmp]
    vars = tmodel+['_grid_equatorial','_grid_vertical']
    
    
    get_data, vars[1], tmp, info
    scr0 = info.r0
    sctrng = reverse(minmax(info.ts)*deg/15)
    sczrng = minmax(info.zs)
    
    ; equatorial plane.
    get_data, vars[0], tmp, info
    tilt0 = info.z0
    eqrrng = minmax(info.rs)
    eqtrng = minmax(info.ts)/rad/15


    ; interpolate grid to every record and increase spatial resolution.
    load = 0
    if file_test(scasifn) eq 0 then load = 1 else tplot_restore, filename = scasifn
    foreach tvar, scasivars do if tnames(tvar) eq '' then load = 1
    ; load = 1
    if load then begin
        foreach tvar, vars do begin
            get_data, tvar, tuts, info
            utr = minmax(tuts)
            scfmlats = info.fmlat
            scfmlts = info.fmlt ; 0-24.
            scfmlts[where(scfmlts gt 12)]-= 24  ; mlt continuous around midnight.
            
            tmp = size(scfmlats,/dimensions)
            nrec0 = tmp[0]
            scfsz0 = tmp[1:*]
            
            
            ; **** interpolate to add grid density.
            nmag = 5
            scfsz1 = scfsz0*nmag-nmag+1
            dat0 = scfmlats
            dat1 = dblarr([nrec0,scfsz1])
            for i = 0, nrec0-1 do begin
                ; interpolate on columns at each row.
                tx0 = findgen(scfsz0[1])/(scfsz0[1]-1)
                tx1 = findgen(scfsz1[1])/(scfsz1[1]-1)
                tdat = dblarr(scfsz0[0],scfsz1[1])
                for j = 0, scfsz0[0]-1 do $
                    tdat[j,*] = interpol(reform(dat0[i,j,*]), tx0, tx1)
                ; interpolate on rows at each column.
                tx0 = findgen(scfsz0[0])/(scfsz0[0]-1)
                tx1 = findgen(scfsz1[0])/(scfsz1[0]-1)
                for j = 0, scfsz1[1]-1 do $
                    dat1[i,*,j] = interpol(reform(tdat[*,j]), tx0, tx1)
            endfor
            scfmlats = dat1
            
            ;dat2 = dblarr([nrec0,scfsz1])
            ;for i = 0, nrec0-1 do begin
            ;    dat2[i,*,*] = congrid(reform(dat0[i,*,*]), scfsz1[0], scfsz1[1], /center, /interp)
            ;endfor
            
            dat0 = scfmlts
            dat1 = dblarr([nrec0,scfsz1])
            for i = 0, nrec0-1 do begin
                ; interpolate on columns at each row.
                tx0 = findgen(scfsz0[1])/(scfsz0[1]-1)
                tx1 = findgen(scfsz1[1])/(scfsz1[1]-1)
                tdat = dblarr(scfsz0[0],scfsz1[1])
                for j = 0, scfsz0[0]-1 do $
                    tdat[j,*] = interpol(reform(dat0[i,j,*]), tx0, tx1)
                ; interpolate on rows at each column.
                tx0 = findgen(scfsz0[0])/(scfsz0[0]-1)
                tx1 = findgen(scfsz1[0])/(scfsz1[0]-1)
                for j = 0, scfsz1[1]-1 do $
                    dat1[i,*,j] = interpol(reform(tdat[*,j]), tx0, tx1)
            endfor
            scfmlts = dat1
            
            
            ; **** interpolate to add temporal resolution.
            uts = smkarthm(utr[0],utr[1],dr0,'dx')
            nrec1 = n_elements(uts)
            dat0 = scfmlats
            dat1 = dblarr([nrec1,scfsz1])
            for i = 0, scfsz1[0]-1 do begin
                for j = 0, scfsz1[1]-1 do begin
                    dat1[*,i,j] = interpol(dat0[*,i,j],tuts,uts,/quadratic)
                endfor
            endfor
            scfmlats = dat1
            
            dat0 = scfmlts
            dat1 = dblarr([nrec1,scfsz1])
            for i = 0, scfsz1[0]-1 do begin
                for j = 0, scfsz1[1]-1 do begin
                    dat1[*,i,j] = interpol(dat0[*,i,j],tuts,uts,/quadratic)
                endfor
            endfor
            scfmlts = dat1
    
            
            
            ; **** map aurora from ionosphere to in situ.
            ; convert scfmlats and scfmlts to scfxs and scfys.
            scfrs = (90-scfmlats)/(90-minlat)
            scfts = (24-scfmlts)*15*rad
            scfxs =-scfrs*sin(scfts)
            scfys =-scfrs*cos(scfts)
            txs = scfxs*imgsz[0]/2+imgsz[1]
            tys = scfys*imgsz[1]+imgsz[1]
            
            get_data, 'asf_mos', tmp, mos, pxidx
            idx = where(tmp ge utr[0] and tmp le utr[1])
            mos = mos[idx,*]
            
            ; the aurora images in tail.
            scasi = fltarr([nrec1,scfsz1])
            for i = 0, nrec1-1 do begin
                timg = fltarr(imgsz)
                timg[pxidx] = mos[i,*]
                for j = 0, scfsz1[0]-1 do begin
                    for k = 0, scfsz1[1]-1 do begin
                        tj = txs[i,j,k]
                        tk = tys[i,j,k]
                        scasi[i,j,k] = timg[tj,tk]
                    endfor
                endfor
            endfor
            
            
            tmp = strsplit(tvar,'_',/extract)
            tmp[1] = 'asf'
            var = strjoin(tmp,'_')
            store_data, var, uts, scasi
            store_data, var+'_mlat', uts, scfmlats
            store_data, var+'_mlt', uts, scfmlts
        endforeach
        
        tplot_save, scasivars, filename = scasifn
    endif


    ; get the upper and lower boundaries of auroral oval.
    vertlim0 = 140  ; threshold for upper limit.
    psbllim0 = 300  ; larger than 200 to be considered.
    
    tvar = tmodel+'_asf_vertical'
    get_data, tvar, uts, timg
    ; lines in situ images and lines in iono images.
    scfsz1 = size(reform(timg[0,*,*]),/dimensions)
    tinfo = {vert:dblarr(scfsz1[0]), psbl:dblarr(scfsz1[0]), equa:dblarr(scfsz1[0])}
    nrec = n_elements(uts)
    infos = replicate(tinfo,nrec)
    
    for i = 0, nrec-1 do begin
        tvar = tmodel+'_asf_vertical'
        get_data, tvar, uts, timg
        
    
        timg = reform(timg[i,*,*])
        xs = findgen(scfsz1[0])
        zs = dblarr(scfsz1[0])      ; the z-value of the mid line.
        midvs = dblarr(scfsz1[0])   ; the aurora brightness of the mid line.
        for j = 0, scfsz1[0]-1 do begin
            midvs[j] = max(timg[j,*],idx)
            zs[j] = idx
        endfor
        
        idx = where(midvs lt psbllim0, cnt)
        if cnt ne 0 then begin
            midvs[idx] = !values.d_nan
            zs[idx] = !values.d_nan
        endif
        infos[i].psbl = zs
    
        
        ; polaward boundary.
        for j = 0, scfsz1[0]-1 do begin
            idx = where(timg[j,*] ge vertlim0, cnt)
            if cnt eq 0 then zs[j] = !values.d_nan else zs[j] = idx[cnt-1]
        endfor
        infos[i].vert = zs
    
    
    
        tvar = tmodel+'_asf_equatorial'
        get_data, tvar, uts, timg
        timg = reform(timg[i,*,*])
        ys = findgen(scfsz1[0])
        xs = dblarr(scfsz1[0])      ; the x-value of the low line.
    
       
        ; earthward boundary.
        for j = 0, scfsz1[0]-1 do begin
            idx = where(timg[j,*] ge vertlim0, cnt)
            if cnt eq 0 then xs[j] = !values.d_nan else xs[j] = idx[0]
        endfor
        infos[i].equa = xs
    endfor



; **** save plots for movie.
    utr = minmax(uts)
    get_data, 'asf_mos', tmp, mos, pxidx
    idx = where(tmp ge utr[0] and tmp le utr[1])
    mos = mos[idx,*]
    
    for i = 0, n_elements(uts)-1 do begin
        
        
        ofn = rootdir+'/asi_iono2tail/asi_iono2tail_'+ $
            time_string(uts[i],tformat='YYYY_MMDD_hhmm_ss')+'.png'
    ;    ofn = 0
        sgopen, ofn, xsize = 12, ysize = 6, /inch
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        device, decomposed = 0
        loadct2, 43
        
        erase, 255
        
        
    ;    tmp = (pos0[3]-pos0[1])*!d.y_size/!d.x_size
    ;    pos0[0] = 0.5-tmp
    ;    pos0[2] = 0.5+tmp
        
        poss = sgcalcpos(2,region=[0.5,0,0.9,1], ypad = 5)
        poss[0,*] = 0.55
        poss[2,*] = 0.90
    
        
        ; saturated color.
        pos1 = [0.05,0.1,0.4,0.45]
        pos0 = [0.05,0.55,0.4,0.9]
    
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        timg = bytscl(timg, max=300, min=mincnt, top=top, /nan)
        sgtv, timg, position = pos1, ct = 1
        
        sgset_map, position = pos1, color = 255, xrange = [-90,90], $
            yrange = [minlat,90], ytickv = [55,65,75], $
            xtickname = ['18','00','06'], xticknudge = [[1.5,-.5],[0.2,1],[-1.5,-.5]]
        xyouts, (pos0[0]+pos0[2])*0.5, pos0[3]+ychsz*1, /normal, alignment = 0.5, $
            'Aurora images in 2 saturation schemes', charsize = 1.5
    
        xyouts, pos0[2]+xchsz*1, pos0[3]-ychsz*1.2, /normal, $
            'Dim', charsize = 1.2
        xyouts, pos1[2]+xchsz*1, pos1[3]-ychsz*1.2, /normal, $
            'Bright', charsize = 1.2
    
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        timg = bytscl(timg, max=maxcnt, min=mincnt, top=top, /nan)
        sgtv, timg, position = pos0, ct = 1
        
        sgset_map, position = pos0, color = 255, xrange = [-90,90], $
            yrange = [minlat,90], ytickv = [55,65,75], $
            xtickname = ['18','00','06'], xticknudge = [[1.5,-.5],[0.2,1],[-1.5,-.5]]
    
            
        xyouts, pos0[0]+xchsz*0.5, pos0[1]+ychsz*0.5, /normal, $
            time_string(uts[i], tformat='hh:mm:ss'), color = 255
            
    
        ; vertical plane.
        tvar = tmodel+'_asf_vertical'
        get_data, tvar+'_mlat', uts, scfmlats
        get_data, tvar+'_mlt', uts, scfmlts
        get_data, tvar, uts, timg
        timg = bytscl(reform(timg[i,*,*]), max=maxcnt, min=mincnt, top=top, /nan)
        
        scfrs = (90-scfmlats)/(90-minlat)
        scfts = (24-scfmlts)*15*rad
        scfxs =-scfrs*sin(scfts)
        scfys =-scfrs*cos(scfts)
        ; draw a box.
        plot, [-1,1], [-1,0], xstyle=5, ystyle=5, /iso, position = pos0, /nodata, /noerase
        plots, scfxs[i,0,*], scfys[i,0,*], psym = 3, color = 6
        plots, scfxs[i,-1,*], scfys[i,-1,*], psym = 3, color = 6
        plots, scfxs[i,*,0], scfys[i,*,0], psym = 3, color = 6
        plots, scfxs[i,*,-1], scfys[i,*,-1], psym = 3, color = 6
        
        
        
        
        xr = -sctrng*15*rad*scr0 ; from MLT to distance in Re.
        yr = sczrng
        ytitl = 'Z (Re)'
        xtitl = 'Azm.Dist (Re)'
        titl = ''
        tpos = poss[*,0]
        aspr = 1d/abs(xr[1]-xr[0])*abs(yr[1]-yr[0])*!d.x_size/!d.y_size
        tpos[3] = tpos[1]+(tpos[2]-tpos[0])*aspr
        tpos = tpos+[0,1,0,1]*(pos0[3]-tpos[3]) ; adjust to the same height as the upper left panel.
    
        xyouts, (tpos[0]+tpos[2])*0.5, tpos[3]+ychsz*1, /normal, alignment = 0.5, $
            'Project aurora to 2 planes in the tail', charsize = 1.5
        
        ; data coord.
        timg = congrid(timg, (tpos[2]-tpos[0])*!d.x_size, (tpos[3]-tpos[1])*!d.y_size)
        sgtv, timg, position = tpos, ct = 1
        
        plot, xr, yr, /nodata, /noerase, position = tpos, /iso, $
            xstyle = 1, xtitle = xtitl, xticklen = -0.02, $
            ystyle = 1, ytitle = ytitl, yticklen = -0.01
        foreach tprobe, probes do begin
            pre0 = 'rbsp'+tprobe+'_'
            scsym = (tprobe eq 'a')? 1: 7
            get_data, pre0+'pos_gsm', tuts, posgsm
            scz = interpol(posgsm[*,2],tuts, uts[i])
            get_data, pre0+'mlt', tuts, posmlt
            scmlt = interpol(posmlt,tuts, uts[i])
            scmlt = (scmlt-24)*15*rad*scr0
            plots, scmlt, scz, psym = scsym, color = 5, symsize = 0.8
        endforeach
        xyouts, 0.5, tpos[3]+ychsz*0.5, /normal, alignment = 0.5, titl
        
    
        ; device coord.
        xr = [0,scfsz1[0]-1]
        yr = [0,scfsz1[1]-1]
        plot, xr, yr, /noerase, /nodata, position = tpos, $
            xstyle = 5, ystyle = 5
        xs = findgen(scfsz1[0])
        ys = infos[i].vert
        plots, xs, ys, psym = 1, color = 6, symsize = 0.1
        ys = infos[i].psbl
        plots, xs, ys, psym = 1, color = 4, symsize = 0.1
    
        xyouts, tpos[2]+xchsz*1.5, tpos[3]-ychsz*1.2*1, /normal, 'Red Box', charsize = 1.2, color = 6
        xyouts, tpos[2]+xchsz*1.5, tpos[3]-ychsz*1.2*2, /normal, 'Vertical plane', charsize = 1
        xyouts, tpos[2]+xchsz*1.5, tpos[3]-ychsz*1.2*3, /normal, 'at L = 5', charsize = 1
        xyouts, tpos[2]+xchsz*3, (tpos[1]-ychsz*1.2*4+tpos[3])*0.5, /normal, 'Up'
        arrow, tpos[2]+xchsz*2, tpos[1], $
            tpos[2]+xchsz*2, tpos[3]-ychsz*1.2*4, /normal, /solid, hsize = 8
    
        
    
        
        ; equatorial plane.
        tvar = tmodel+'_asf_equatorial'
        get_data, tvar+'_mlat', uts, scfmlats
        get_data, tvar+'_mlt', uts, scfmlts
        get_data, tvar, uts, timg
        timg = bytscl(reform(timg[i,*,*]), max=maxcnt, min=mincnt, top=top, /nan)
    
        scfrs = (90-scfmlats)/(90-minlat)
        scfts = (24-scfmlts)*15*rad
        scfxs =-scfrs*sin(scfts)
        scfys =-scfrs*cos(scfts)
        ; draw a box.
        plot, [-1,1], [-1,0], xstyle=5, ystyle=5, /iso, position = pos0, /nodata, /noerase
        plots, scfxs[i,0,*], scfys[i,0,*], psym = 3, color = 4
        plots, scfxs[i,-1,*], scfys[i,-1,*], psym = 3, color = 4
        plots, scfxs[i,*,0], scfys[i,*,0], psym = 3, color = 4
        plots, scfxs[i,*,-1], scfys[i,*,-1], psym = 3, color = 4
        
        xr = eqtrng-24
        yr = eqrrng
        xtitl = 'MLT (hr)'
        ytitl = 'R (Re)'
        titl = ''
        tpos = poss[*,1]
        aspr = 1d/abs(xr[1]-xr[0])*abs(yr[1]-yr[0])*!d.x_size/!d.y_size
        tpos[3] = tpos[1]+(tpos[2]-tpos[0])*aspr
        tpos = tpos+[0,1,0,1]*(pos1[1]-tpos[1]) ; adjust to the same height as the upper left panel.
    
        plot, xr, yr, /nodata, /noerase, position = tpos, /iso, $
            xstyle = 1, xtitle = xtitl, xticklen = -0.02, $
            ystyle = 1, ytitle = ytitl, yticklen = -0.01
        xyouts, 0.5, tpos[3]+ychsz*0.5, /normal, alignment = 0.5, titl
        timg = congrid(timg, (tpos[2]-tpos[0])*!d.x_size, (tpos[3]-tpos[1])*!d.y_size)
        sgtv, timg, position = tpos, ct = 1
        
        
        xr = [0,scfsz1[0]-1]
        yr = [0,scfsz1[1]-1] & yr = (yr)
        plot, xr, yr, /noerase, /nodata, position = tpos, $
            xstyle = 5, ystyle = 5, yrange = yr
        xs = findgen(scfsz1[0])
        ys = infos[i].equa
        plots, xs, ys, psym = 1, color = 4, symsize = 0.1
        
    
    
        device, decomposed = 1
        xyouts, tpos[2]+xchsz*1.5, tpos[3]-ychsz*1.2*1, /normal, 'Green Box', charsize = 1.2, color = sgcolor('green')
        xyouts, tpos[2]+xchsz*1.5, tpos[3]-ychsz*1.2*2, /normal, 'Equatorial plane', charsize = 1
        xyouts, tpos[2]+xchsz*1.5, tpos[3]-ychsz*1.2*3, /normal, 'Tilt = 13 deg', charsize = 1
        xyouts, tpos[2]+xchsz*3, (tpos[3]-ychsz*1.2*4+tpos[1])*0.5, /normal, 'Earthward'
        arrow, tpos[2]+xchsz*2, tpos[3]-ychsz*1.2*4, $
            tpos[2]+xchsz*2, tpos[1], /normal, /solid, hsize = 8
        sgclose
    endfor


    ; pic to movie.
    vfn = rootdir+'/thg_2013_0607_insitu_v3.mp4'

    ext = 'png'
    spic2movie, rootdir+'/', vfn, ext

end
