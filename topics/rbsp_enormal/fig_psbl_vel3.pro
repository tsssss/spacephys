;----praparation.
    ; load data.
    _2013_0607_load_data
    
    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 4

    ; constants.
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1
    
    utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
    probes = ['a','b']
    tmodel = 't01'
    labchsz = 1.5       ; set label charsize larger.
    strsim = '!9'+string(126b)+'!X'



    ; settings.
    autr = time_double(['2013-06-07/04:46','2013-06-07/05:05'])
    butr = time_double(['2013-06-07/04:56','2013-06-07/05:00'])
    probes = ['a','b']
    red = sgcolor('red')
    





;----from fig_keo_ewo_insitu.pro
    ; settings.
    top = 254
    maxcnt = 800
    maxcnt = 500
    mincnt = 50
    dr0 = 3d
        
        
    ;--load data.
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
    
    
        scasifn = shomedir()+'/psbl_de_32hz/2013_0607_scasi_'+tmodel+'.tplot'
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
    
    
    
;----plot 2: aurora velocity in the perpendicular plane.



    ;--calc velocity of the west surge.
    figysz = 51     ; # of pixels in y.
    figxsz = 301    ; # of pixels in x.
    xr = -sctrng
    yr =  sczrng
    zr = [50,400]   ; color range for brightness.
    vert0 = 150     ; threshold for west boundary.
    figys = smkarthm(yr[0], yr[1], figysz, 'n')
    figxs = smkarthm(xr[0], xr[1], figxsz, 'n')


    ; load data.
    tvar = tmodel+'_asf_vertical'
    get_data, tvar, uts, timg
    dr0 = sdatarate(uts)

    
    
    ofn = shomedir()+'/fig_psbl_vel2.pdf'
    ;ofn = 0
    sgopen, ofn, xsize = 6, ysize = 5, /inch
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    hsize = (size(ofn,/type) eq 7)? 200: 10
    thick = (size(ofn,/type) eq 7)? 10: 4
    

    ; width are the same, heigh is determined by groups.
    ; 1&2 are aligned to their middle position.
    pos1 = [0.15,0.75,0.90,0]
    pos2 = [0.15,0,0.90,0.75]
    pos3 = [0.15,0.25,0.90,0]
    pos4 = [0.15,0,0.90,0.25]

    xtklen = -0.04
    xtickv = [-2,-1,0,1]
    xticks = n_elements(xtickv)-1
    xminor = 10
    xtickn = ['22','23','00','01']

    ytickv = [1.5,2.5,3.5,4.5]
    yticks = n_elements(ytickv)-1
    yminor = 4
    ytklen = -0.01



;---snapshot of the westward surge.
    tpos = pos1
    tpos[3] = tpos[1]+(tpos[2]-tpos[0])/figxsz*figysz*!d.x_size/!d.y_size
    tpos = tpos+[0,1,0,1]*ychsz*0.2

    ; settings.
    txr = [-3.5,-2]
    tutr = time_double(['2013-06-07/04:54','2013-06-07/04:56:30'])
    velutr = time_double(['2013-06-07/04:54:45','2013-06-07/04:55:24'])

    ; check between these y-locations.
    tys = smkarthm(3.6, 3.1, -0.1, 'dx')
    nty = n_elements(tys)
    idx = where(uts ge velutr[0] and uts le velutr[1])
    idx1 = median(idx)
    tut1 = uts[idx1]    ; middle time when u_west is calculated.

    ; draw snapshot.
    tmp = bytscl(reform(timg[idx1,*,*]), min = zr[0], max = zr[1], top = 254)>8
    sgtv, tmp, position = tpos, ct = 1, file = 'ct2'

    plot, xr, yr, /noerase, /nodata, position=tpos, $
        xstyle = 1, xticklen=xtklen, xtitle='', xtickformat='(A1)', $
        xminor=xminor, xtickv=xtickv, xticks=xticks, xtickname=xtickn, $
        ystyle = 1, yticklen=ytklen, ytitle = 'Z GSM (Re)', $
        yticks=yticks, yminor=yminor, ytickv=ytickv

    ; add positions of the slices.
    for i = 0, nty-1 do plots, xr, tys[i]+[0,0], color = sgcolor('red')
    plots, xr, tys[nty/2]+[0,0], color=sgcolor('red'), thick=thick*0.5
    xyouts, tpos[2]-xchsz*1, tpos[3]-ychsz*1.2, alignment = 1, /normal, $
        time_string(tut1), color = sgcolor('white')
    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, $
        alignment = 0, 'Aurora mapped to the vertical plane', charsize=labchsz


    ; draw arrow.
    u_west = 150
    tx = 22.2-24
    ty = 3.2
    tmp = convert_coord(tx, ty, /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    arrow, tx, ty, tx-xchsz*5, ty, /normal, /solid, hsize=hsize, $
        color=sgcolor('white'), thick=thick*0.5
    xyouts, tx-xchsz*2.5, ty+ychsz*0.6, /normal, color=sgcolor('white'), $
        '<u!Dwest!N> '+strsim+' '+sgnum2str(u_west,ndec=0)+' km/s', alignment=0.5, charsize=labchsz*0.8
    xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.5, /normal, 'a', color=sgcolor('white'), charsize=labchsz

;---snapshot of the poleward expansion.
    tpos = pos2
    tpos[1] = tpos[3]-(tpos[2]-tpos[0])/figxsz*figysz*!d.x_size/!d.y_size
    tpos = tpos-[0,1,0,1]*ychsz*0.2

    ; settings.
    tutr = time_double(['2013-06-07/04:54','2013-06-07/04:56:30'])
    velutr = time_double(['2013-06-07/04:54:54','2013-06-07/04:55:39'])

    ; check between these x-locations.
    txs = -smkarthm(2.4, 2.9, 0.1, 'dx')/(15*rad*scr0)
    ntx = n_elements(txs)
    idx = where(uts ge velutr[0] and uts le velutr[1])
    idx2 = median(idx)
    tut2 = uts[idx2]    ; middle time when u_up is calculated.

    ; draw snapshot.
    tmp = bytscl(reform(timg[idx2,*,*]), min = zr[0], max = zr[1], top = 254)>8
    sgtv, tmp, position = tpos, ct = 1, file = 'ct2'

    plot, xr, yr, /noerase, /nodata, position = tpos, $
        xstyle = 1, xticklen=xtklen, xtitle='MLT (hr)', $
        xminor=xminor, xtickv=xtickv, xticks=xticks, xtickname=xtickn, $
        ystyle = 1, yticklen=ytklen, ytitle = 'Z GSM (Re)', $
        yticks=yticks, yminor=yminor, ytickv=ytickv

    for i = 0, ntx-1 do plots, txs[i]+[0,0], yr, color = sgcolor('red')
    plots, txs[ntx/2]+[0,0], yr, color=sgcolor('red'), thick=thick*0.5
    xyouts, tpos[2]-xchsz*1, tpos[3]-ychsz*1.2, alignment = 1, /normal, $
        time_string(tut2), color = sgcolor('white')
    
    ; draw arrow.
    u_up = 130
    tx = 22.4-24
    ty = 3
    tmp = convert_coord(tx, ty, /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    arrow, tx, ty, tx, ty+ychsz*2, /normal, /solid, hsize=hsize, $
        color=sgcolor('white'), thick=thick*0.5
    polyfill, tx+xchsz*[1,15,15,1], ty+ychsz*[0,0,1.4,1.4], /normal, color=sgcolor('black')
    xyouts, tx+xchsz*1, ty+ychsz*0.5, /normal, color=sgcolor('white'), $
        '<u!Dup!N> '+strsim+' '+sgnum2str(u_up,ndec=0)+' km/s', charsize=labchsz*0.8
    xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.5, /normal, 'b', color=sgcolor('white'), charsize=labchsz


;---time series of the westward surge.
    device, decomposed=0
    loadct2, 1
    tpos = pos3
    tpos[3] = tpos[1]+(tpos[2]-tpos[0])/figxsz*figysz*!d.x_size/!d.y_size
    tpos = tpos+[0,1,0,1]*ychsz*0.2

    tty = tys[nty/2]
    tidx = (tty-yr[0])/(yr[1]-yr[0])*figysz
    tyr = xr[0]+[0,1.5]+24
    ewos = reform(timg[*,*,tidx])
    
;---new way use contour 2018-03-28.
    tty = 3.0
    dty = 0.07
    yys = smkarthm(yr[0],yr[1],51,'n')
    tyr = xr[0]+[0,1.5]+24
    tidx = where(yys ge tty-dty and yys le tty+dty)
    ewos = reform(total(timg[*,*,tidx],3)/n_elements(tidx))
    tvar = 'ewo_tmp'
    store_data, tvar, uts, ewos, figxs+24, $
        limits={ytitle:'MLT (hr)', spec:1, no_interp:1, $
        xstyle:1, xticklen:xtklen, $
        ystyle:1, yrange:tyr, yticklen:ytklen, $
        yticks:2, ytickv:[21.5,22.0,22.5], yminor:5, $
        ztitle:'Photon Count', zticks:3, ztickv:[1,2,3,4]*100, zrange:[50,400]}

    zrng = [50,250]
    levels = [0,smkarthm(zrng[0],zrng[1],20,'dx')]
    nlevel = n_elements(levels)
    colors = smkarthm(0,top, nlevel,'n')
    contour, ewos, uts, figxs+24, position=tpos, levels=levels, c_colors=colors, /fill, /noerase, $
        xstyle=1, xrange=tutr, xticklen=xtklen, $
        ystyle=1, yrange=tyr, yticklen=ytklen, yticks=2, ytickv=[21.5,22.0,22.5], yminor=5, xtickformat='(A1)'
    ttpos = tpos & ttpos[0] = tpos[2]+xchsz*1 & ttpos[2] = tpos[2]+xchsz*2
    sgcolorbar, colors, position=ttpos, zrange=zrng, zticks=2
    
    ; make a slice?
    tut = time_double(['2013-06-07/04:55'])
    idx = where(uts ge tut-3 and uts le tut+3)
    slice = total(ewos[idx,*,*],1)
    
    
    plot, tutr, tyr, position=tpos, /nodata, /noerase, xstyle=5, ystyle=5
    plots, tut1+[0,0], tyr, color=6

    velutr = time_double(['2013-06-07/04:54:45','2013-06-07/04:55:24'])
    verts = fltarr(nrec)
    for j = 0, nrec-1 do begin
        tdat = (reform(ewos[j,*]))   ; a slice at each time record.
        idx = where(tdat ge vert0, cnt) ; find the lowest point cross threshold.
        verts[j] = (cnt eq 0)? !values.f_nan: figxs[idx[0]]
    endfor
    pos0 = verts[where(uts eq velutr[0])]+24
    pos1 = verts[where(uts eq velutr[1])]+24
    tvel = -(pos1-pos0)*re/(velutr[1]-velutr[0])  ; in km/s.
    plots, velutr, [pos0,pos1], color = 6, thick=thick
    xyouts, tpos[0]+xchsz*2.5, tpos[1]+ychsz*1.5, /normal, $
        'Westward motion !CZ GSM = '+sgnum2str(tty, ndec=1)+' Re', color = 255
    xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*1.5, /normal, 'c', color=sgcolor('white'), charsize=labchsz

    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, $
        alignment = 0, 'Aurora motion at fixed locations in the vertical plane', charsize=labchsz
    


;---time series of the poleward expansion.
    tpos = pos4
    tpos[1] = tpos[3]-(tpos[2]-tpos[0])/figxsz*figysz*!d.x_size/!d.y_size
    tpos = tpos-[0,1,0,1]*ychsz*0.2

    ttx = txs[ntx/2]
    tyr = sczrng
    tidx = (ttx-xr[0])/(xr[1]-xr[0])*figxsz
    keos = reform(timg[*,tidx,*])
    tvar = 'keo_tmp'
    store_data, tvar, uts, keos, figys, $
        limits={ytitle:'Z GSM (Re)', spec:1, no_interp:1, $
        xstyle:1, xticklen:xtklen, $
        ystyle:1, yrange:tyr, yticklen:ytklen, $
        yticks:3, yminor:4, $
        ztitle:'Photon Count', zticks:3, ztickv:[1,2,3,4]*100, zrange:[50,400]}

    ;tplot, tvar, position=tpos, /noerase, /novtitle, trange=tutr
    
    zrng = [50,250]
    levels = [0,smkarthm(zrng[0],zrng[1],20,'dx')]
    nlevel = n_elements(levels)
    colors = smkarthm(0,top, nlevel,'n')
    contour, keos, uts, figys, position=tpos, levels=levels, c_colors=colors, /fill, /noerase, $
        xstyle=1, xrange=tutr, xticklen=xtklen, $
        ystyle=1, yrange=tyr, yticklen=ytklen, yticks=3, yminor=5, xtickformat='(A1)'
    ttpos = tpos & ttpos[0] = tpos[2]+xchsz*1 & ttpos[2] = tpos[2]+xchsz*2
    sgcolorbar, colors, position=ttpos, zrange=zrng, zticks=2
    
    plot, tutr, tyr, position=tpos, /nodata, /noerase, xstyle=5, ystyle=5
    plots, tut2+[0,0], tyr, color=6
    
    
    velutr = time_double(['2013-06-07/04:54:54','2013-06-07/04:55:39'])
    verts = fltarr(nrec)
    for j = 0, nrec-1 do begin
        tdat = (reform(keos[j,*]))   ; a slice at each time record.
        idx = where(tdat ge vert0, cnt) ; find the highest point cross threshold.
        verts[j] = (cnt eq 0)? !values.f_nan: figys[idx[cnt-1]]
    endfor
    pos0 = verts[where(uts eq velutr[0])]
    pos1 = verts[where(uts eq velutr[1])]
    tvel = -(pos1-pos0)*re/(velutr[1]-velutr[0])  ; in km/s.
    plots, velutr, [pos0,pos1], color = 6, thick=thick
    xyouts, tpos[0]+xchsz*2.5, tpos[3]-ychsz*1, /normal, $
        'Upward motion at!CMLT = '+sgnum2str(ttx+24, ndec=1)+' hr', color = 255
    xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*1, /normal, 'd', color=sgcolor('white'), charsize=labchsz


;---finish plot.
    sgclose

end
