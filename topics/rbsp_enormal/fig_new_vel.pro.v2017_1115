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


    ofn = 0
    ofn = shomedir()+'/fig_new_vel.pdf'
    sgopen, ofn, xsize=9, ysize=6, /inch

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    hsize = (size(ofn,/type) eq 7)? 200: 10
    thick = (size(ofn,/type) eq 7)? 10: 4

    ; plot for s/c position.
    pos1 = [0.1,0.15,0.4,0.9]
    pos2 = [0.45,0.15,0.92,0.38]
    pos3 = [0.45,0.52,0.92,0.9]


;----from fig_new_config.pro.
    device, decomposed=1
    tpos = pos1

    
    xr = [0,-6]
    xtickv = [0,-5]
    xticks = n_elements(xtickv)-1
    xminor = 5
    xtitle = 'X (Re)'
    xtklen = 0.01

    yr = [6,-6]
    ytickv = [-5,0,5]
    yticks = n_elements(ytickv)-1
    yminor = 5
    ytitle = 'Y (Re)'
    ytklen = 0.03

    plot, xr, yr, /nodata, /iso, /noerase, position = tpos, $
        xstyle = 1, xticks = 1, xrange = xr, xminor = xminor, xtitle = xtitle, xticklen = xtklen, xtickv = xtickv, $
        ystyle = 1, yticks = 1, yrange = yr, yminor = yminor, ytitle = ytitle, yticklen = ytklen, ytickv = ytickv
    
    ; add L=5.
    tmp = findgen(101)/100d
    txs = -cos(tmp*!dpi-!dpi*0.5)*5
    tys =  sin(tmp*!dpi-!dpi*0.5)*5
    plots, txs, tys, linestyle = 1
    xyouts, -1, 5, /data, 'R=5'
    
    ; add MLT=22.
    tmp = findgen(101)/100d*10+1
    txs = -tmp*cos(0*15*rad)
    tys =  tmp*sin(0*15*rad)
    oplot, txs, tys, linestyle = 1
    tmp = convert_coord(-2.5, 0.1, /data, /to_normal)
    xyouts, tmp[0], tmp[1]+ychsz*0.5, /normal, alignment=0, 'Mid-night'
    
    ; add earth.
    tmp = findgen(101)/100d
    txs = -cos(tmp*!dpi-!dpi*0.5)
    tys =  sin(tmp*!dpi-!dpi*0.5)
    ; polyfill, txs, tys, /line_fill, orientation=-45
    polyfill, txs, tys, /data, color=sgcolor('grey')
    plots, txs, tys
    
    ; add arrow.
    v_in = 23d
    v_azim = 25d
    coef = 0.06
    
    get_data, 'rbspa_pos_gsm', uts, rgsm
    ax0 = interpol(rgsm[*,0], uts, utr[0])
    ay0 = interpol(rgsm[*,1], uts, utr[0])
    get_data, 'rbspb_pos_gsm', uts, rgsm
    bx0 = interpol(rgsm[*,0], uts, utr[0])
    by0 = interpol(rgsm[*,1], uts, utr[0])
    
    tcolor = sgcolor('black')
    tx0 = ax0
    ty0 = ay0
    tx1 = tx0+(bx0-ax0)/snorm([bx0,by0]-[ax0,ay0])*v_azim*coef
    ty1 = ty0+(by0-ay0)/snorm([bx0,by0]-[ax0,ay0])*v_azim*coef
    arrow, tx0,ty0, tx1,ty1, /data, /solid, hsize=hsize, color=tcolor, thick=thick*0.5
    tmp = convert_coord(tx1,ty1, /data, /to_normal)
    xyouts, tmp[0], tmp[1]-ychsz, alignment=0.5, /normal, 'v!Dazim', color=tcolor, charsize=1.2
    
    tcolor = sgcolor('grey')
    tx0 = ax0
    ty0 = ay0
    tx1 = tx0+cos(atan(ay0,-ax0))*v_in*coef
    ty1 = ty0-sin(atan(ay0,-ax0))*v_in*coef
    arrow, tx0,ty0, tx1,ty1, /data, /solid, hsize=hsize, color=tcolor, thick=thick*0.5
    tmp = convert_coord(tx1,ty1, /data, /to_normal)
    xyouts, tmp[0], tmp[1]+ychsz*0.5, alignment=0, /normal, 'v!Din', color=tcolor, charsize=1.2


    ; plot s/c trajactory.
    tutr = utr[0]+[0,3]
    foreach tprobe, probes do begin
        tcolor = (tprobe eq 'a')? sgcolor('red'): sgcolor('blue')
        get_data, 'rbsp'+tprobe+'_pos_gsm', uts, rgsm
        idx = where(uts ge tutr[0] and uts le tutr[1])
        rgsm = rgsm[idx,*]
        plots, rgsm[*,0], rgsm[*,1], color = tcolor, thick = thick
        tmp = convert_coord(rgsm[0,0],rgsm[0,1], /data, /to_normal)
        xyouts, tmp[0]+xchsz*0.5, tmp[1]-ychsz*0.2, strupcase(tprobe), /normal, color=tcolor
        tmp = findgen(21)/20*2*!dpi
        txs = cos(tmp)
        tys = sin(tmp)
        usersym, txs, tys, /fill, color = tcolor
        plots, rgsm[0,0], rgsm[0,1], psym = 8, symsize = 0.35
    endforeach
    
    
    ; re-draw the box.
    plot, xr, yr, /nodata, /iso, /noerase, position = tpos, $
        xstyle = 1, xticks = xticks, xrange = xr, xminor = xminor, xtitle = xtitle, xticklen = xtklen, xtickv = xtickv, $
        ystyle = 1, yticks = yticks, yrange = yr, yminor = yminor, ytitle = ytitle, yticklen = ytklen, ytickv = ytickv

    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, 'a. Flow velocity', charsize=labchsz
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*labchsz, /normal, $
        'A schematic view in the plane Z!DGSM!N = 0'


;----from fig_2013_0607_vsc_corr.pro
    ; settings.
    autr = time_double(['2013-06-07/04:46','2013-06-07/05:05'])
    butr = time_double(['2013-06-07/04:56','2013-06-07/05:00'])
    probes = ['a','b']
    red = sgcolor('red')
    
    ; get Vsc.
    tvar = 'vsc'
    get_data, 'rbspa_'+tvar, uts, avsc
    get_data, 'rbspb_'+tvar, uts, bvsc
    idx = where(uts ge butr[0] and uts le butr[1], tnrec)
    tmp = bvsc[idx]
    tmp = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
    bvsc = tmp
    buts = uts[idx]

    ; settings for time shift.
    dr0 = sdatarate(uts)*10
    
    dtrg = [-200,200]           ; the time range for time shift.
    dtrg = [-200,50]           ; the time range for time shift.
    dnrg = ceil(dtrg/dr0)       ; corresponding # of records.
    dns = smkarthm(dnrg[0],dnrg[1],1,'dx')  ; each shift in # of record.
    dts = dns*dr0
    ndt = n_elements(dts)
    corrs = fltarr(ndt)         ; the cross correlation at each shift.


    
    for i = 0, ndt-1 do begin
        get_data, 'rbspa_'+tvar, uts, avsc
        tmp = avsc[where(uts ge butr[0]+dts[i] and uts le butr[1]+dts[i])]
        tmp = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
        avsc = tmp
        corrs[i] = c_correlate(avsc, bvsc, 0)
        corrs[i]*= stddev(avsc)
    endfor
    
    get_data, 'rbspa_'+tvar, uts, avsc
    corrs *= 1d/stddev(bvsc)

    maxcorr = max(corrs, idx)
    maxcorrdt = dts[idx]
    sigcorr = 2/sqrt(n_elements(bvsc))



    ;--plot correlation.
    xr = [200,2000]
    yr = [-50,5]
    yminor = 2
    tpos = pos2 & tpos[2] = tpos[0]+(pos2[2]-pos2[0])*0.48
    dvsc = -20  ; shift -B by -20 V.
    
    tutr = butr+[-1.5,0.5]*300
    get_data, 'rbspa_'+tvar, uts, avsc
    idx = where(uts ge tutr[0] and uts le tutr[1])
    txs = uts[idx]-uts[idx[0]]
    tys = avsc[idx]
    
    plot, txs, tys, /noerase, position = tpos, $
        xstyle=1, xticks=3, xminor=6, $
        xtitle='Second from '+time_string(uts[idx[0]]), $
        ystyle=9, yrange=yr, ytitle='RBSP-A Vsc (V)', yminor=yminor
    axis, max(txs), /yaxis, yrange = yr-dvsc, $
        ystyle=1, ytitle = 'RBSP-B Vsc (V)!Cshifted -20 V', yminor=yminor
    
    get_data, 'rbspb_'+tvar, uts, bvsc
    idx = where(uts ge tutr[0]-maxcorrdt and uts le tutr[1]-maxcorrdt)
    tys = bvsc[idx]+dvsc
    oplot, txs, tys, color = red
    
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.2-ychsz*1.0*1, /normal, 'RBSP-A'
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.2-ychsz*1.0*2, /normal, $
        'RBSP-B, dT = '+sgnum2str(maxcorrdt,nsgn=4)+' sec', color = red
    
    
    tpos = pos2 & tpos[0] = tpos[0]+(pos2[2]-pos2[0])*0.75
    tx = dts
    ty = corrs
    xtickv = [-150,-50,50]
    xticks = n_elements(xtickv)-1
    xminor = 2
    plot, tx, ty, position = tpos, /noerase, /ynozero, $
        ytitle = 'Cross Correlation', yrange = [0,1], yminor=2, $
        xstyle=1, xrange=dtrg, xtitle = 'Time shift on RBSP-B (sec)', $
        xtickv=xtickv, xticks=xticks, xminor=xminor
        
    yr = !y.crange
    plots, maxcorrdt+[0,0], yr, linestyle = 1
    plots, !x.crange, sigcorr+[0,0], linestyle = 2
    
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.2-ychsz*1.0*1, /normal, $
        'dT = '+sgnum2str(maxcorrdt,nsgn=4)+' sec'
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.2-ychsz*1.0*2, /normal, $
        'Max Corr = '+sgnum2str(maxcorr,nsgn=2)
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.2-ychsz*1.0*3, /normal, $
        'Sign.level = '+sgnum2str(sigcorr,nsgn=1)

    tpos = pos2
    xyouts, tpos[0], tpos[3]+ychsz*1, /normal, charsize=labchsz, $
        'd. v!Dazim!N = 25 km/s from timing b/w RBSP-A and -B'





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

    ;----calc velocity of the injections.
        figysz = 51     ; # of pixels in y.
        figxsz = 301    ; # of pixels in x.
        xr = eqtrng     ; MLT, in hr.
        yr = eqrrng     ; R, in Re.
        zr = [50,400]   ; color range for brightness.
        vert0 = 150     ; threshold for vertical boundary.
        figys = smkarthm(yr[0], yr[1], figysz, 'n')
        
        ; settings.
        tutr = time_double(['2013-06-07/04:54','2013-06-07/05:00'])
        tutr = utr
        velutr = time_double(['2013-06-07/04:55:27','2013-06-07/04:58:18'])
        ; check between these x-locations.
        txs = 24-smkarthm(1.6,2.1, 0.1, 'dx')
        ntx = n_elements(txs)

        
        ; load data.
        tvar = tmodel+'_asf_equatorial'
        get_data, tvar, uts, timg
        dr0 = sdatarate(uts)
    
    
        
        ;--plot overview and keogram locations.
        device, decomposed = 1
        xtklen = -0.04
        xtickv = [22,23,24,25]
        xticks = n_elements(xtickv)-1
        xminor = 10
        xtickn = ['22','23','00','01']
        ytklen = -0.01
        ytickv = [4.5,5.5,6.5]
        yticks = n_elements(ytickv)-1
        yminor = 4
        txidx = ntx/2       ; show this slice.

        ;--draw snapshot of the equatorial plane.
        tpos = pos3
        tpos[1] = tpos[3]-(tpos[2]-tpos[0])/figxsz*figysz*!d.x_size/!d.y_size
        posb = tpos

        idx = where(uts ge velutr[0] and uts le velutr[1])
        idx = median(idx)
        txut = uts[idx]

        tmp = bytscl(reform(timg[idx,*,*]), min = zr[0], max = zr[1], top = 254)>8
        sgtv, tmp, position = tpos, ct = 43, file = 'ct2'
        plot, xr, yr, /noerase, /nodata, position=tpos, $
            xstyle = 1, xticklen=xtklen, xtitle='MLT', $
            xminor=xminor, xtickv=xtickv, xticks=xticks, xtickname=xtickn, $
            ystyle = 1, yticklen=ytklen, ytitle='R (Re)', $
            yticks=yticks, yminor=yminor

        ; add positions of the slices.
        for i = 0, ntx-1 do plots, txs[i]+[0,0], yr, color = sgcolor('red')
        plots, txs[txidx]+[0,0], yr, color = sgcolor('red'), thick=thick
        xyouts, tpos[2]-xchsz*1, tpos[1]+ychsz*0.2, alignment = 1, /normal, $
            time_string(txut), color = sgcolor('white')
        
        xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, $
            alignment = 0, 'b. Aurora mapped to the equatorial plane', charsize=labchsz
   
            
    
        ;--plot a middle slice.
        device, decomposed=0
        loadct2, 43
        tpos = pos3 & tpos[3] = tpos[1]+(tpos[3]-tpos[1])*0.35
        
        vels = dblarr(ntx)
        
        for i = 0, ntx-1 do begin
            tx = txs[i]
            tidx = (tx-xr[0])/(xr[1]-xr[0])*301
            keos = reform(timg[*,tidx,*])
            tvar = 'inj_'+string(i,format='(I0)')
            store_data, tvar, uts, keos, figys, $
                limits = {ytitle:'R (Re)', spec:1, no_interp:1, $
                xstyle:1, xticklen:xtklen, $
                ystyle:1, yticks:yticks, ytickv:ytickv, yminor:yminor, yticklen:ytklen, $
                ztitle:'Photon Count', zticks:3, ztickv:[1,2,3,4]*100, zrange:[50,400]}
        
            if i eq txidx then begin
                tplot, tvar, position = tpos, /noerase, /novtitle, trange = tutr
                plot, tutr, yr, position = tpos, /nodata, /noerase, $
                    xstyle=5, ystyle=5
                plots, txut+[0,0], yr, color=6
            endif
        
            verts = fltarr(nrec)
            for j = 0, nrec-1 do begin
                tdat = (reform(keos[j,*]))   ; a slice at each time record.
                idx = where(tdat ge vert0, cnt) ; find the highest point cross threshold.
                verts[j] = (cnt eq 0)? !values.f_nan: figys[idx[0]]
            endfor
            tvar = 'inj_'+string(i,format='(I0)')+'_pos'
            store_data, tvar, uts, verts, limits = {yrange:yr}
        
            pos0 = verts[where(uts eq velutr[0])]
            pos1 = verts[where(uts eq velutr[1])]
            vels[i] = -(pos1-pos0)*re/(velutr[1]-velutr[0])  ; in km/s.
            ty = (pos0+pos1)*0.5
            if i eq txidx then begin
                plots, velutr, [pos0,pos1], color = 6, thick=10
                xyouts, tpos[0]+xchsz*1.5, tpos[1]+ychsz*0.2, /normal, $
                    'MLT = '+sgnum2str(tx, ndec=1), color = 255
            endif
        endfor
        
        xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, $
            alignment = 0, 'c. Earthward motion of the mapped aurora', charsize = labchsz


        ; go back and add avg velocity to panel b.
        v_in = mean(vels)
        tpos = posb
        arrow, tpos[0]+(tpos[2]-tpos[0])*0.3, tpos[3]-ychsz*1, $
            tpos[0]+(tpos[2]-tpos[0])*0.3, tpos[1]+ychsz*1, $
            /normal, color=sgcolor('white'), /solid, hsize=hsize, thick=thick*0.5
        xyouts, tpos[0]+(tpos[2]-tpos[0])*0.3+xchsz*1, (tpos[1]+tpos[3])*0.5, $
            /normal, color=sgcolor('white'), '<V!Din!N> = '+ $
            sgnum2str(v_in,ndec=0)+' km/s'


    ;--finish the plot.
    sgclose
    
    
    
    
    
;----plot 2: aurora velocity in the perpendicular plane.



    ;--calc velocity of the west surge.
    figysz = 51     ; # of pixels in y.
    figxsz = 301    ; # of pixels in x.
    xr = -sctrng
    yr =  sczrng
    zr = [50,400]   ; color range for brightness.
    vert0 = 150     ; threshold for west boundary.
    figys = smkarthm(yr[0], yr[1], figxsz, 'n')


    ; load data.
    tvar = tmodel+'_asf_vertical'
    get_data, tvar, uts, timg
    dr0 = sdatarate(uts)

    
    
    ofn = shomedir()+'/fig_new_vel2.pdf'
    ;ofn = 0
    sgopen, ofn, xsize = 6, ysize = 3, /inch
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    hsize = (size(ofn,/type) eq 7)? 200: 10
    thick = (size(ofn,/type) eq 7)? 10: 4
    

    pos1 = [0.15,0.45,0.95,0.75]
    pos2 = [0.15,0.15,0.95,0.45]

    xtklen = -0.04
    xtickv = [-2,-1,0,1]
    xticks = n_elements(xtickv)-1
    xminor = 10
    xtickn = ['22','23','00','01']

    ytickv = [1.5,2.5,3.5,4.5]
    yticks = n_elements(ytickv)-1
    yminor = 4
    ytklen = -0.01



    ;--westward surge.
    tpos = pos1
    tpos[1] = tpos[3]-(tpos[2]-tpos[0])/figxsz*figysz*!d.x_size/!d.y_size

    ; settings.
    txr = [-3.5,-2]
    tutr = time_double(['2013-06-07/04:54','2013-06-07/04:57'])
    velutr = time_double(['2013-06-07/04:54:45','2013-06-07/04:55:24'])

    ; check between these x-locations.
    tys = smkarthm(3.6, 3.1, -0.1, 'dx')
    nty = n_elements(tys)
    idx = where(uts ge velutr[0] and uts le velutr[1])
    idx1 = median(idx)
    tut1 = uts[idx1]    ; middle time when u_west is calculated.

    ; draw snapshot.
    tmp = bytscl(reform(timg[idx1,*,*]), min = zr[0], max = zr[1], top = 254)>8
    sgtv, tmp, position = tpos, ct = 43, file = 'ct2'

    plot, xr, yr, /noerase, /nodata, position=tpos, $
        xstyle = 1, xticklen=xtklen, xtitle='', xtickformat='(A1)', $
        xminor=xminor, xtickv=xtickv, xticks=xticks, xtickname=xtickn, $
        ystyle = 1, yticklen=ytklen, ytitle = 'Z GSM (Re)', $
        yticks=yticks, yminor=yminor, ytickv=ytickv

    ; add positions of the slices.
    for i = 0, nty-1 do plots, xr, tys[i]+[0,0], color = sgcolor('red')
    xyouts, tpos[2]-xchsz*1, tpos[3]-ychsz*1.2, alignment = 1, /normal, $
        time_string(tut1), color = sgcolor('white')
    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, $
        alignment = 0, 'Aurora mapped to the vertical plane', charsize=labchsz


    ; draw arrow.
    u_west = 150
    tx = 22.2-24
    ty = 3.5
    tmp = convert_coord(tx, ty, /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    arrow, tx, ty, tx-xchsz*5, ty, /normal, /solid, hsize=hsize, $
        color=sgcolor('white'), thick=thick*0.5
    xyouts, tx-xchsz*2.5, ty+ychsz*0.6, /normal, color=sgcolor('white'), $
        '<u!Dwest!N> '+strsim+' '+sgnum2str(u_west,ndec=0)+' km/s', alignment=0.5, charsize=labchsz*0.8
    xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.5, /normal, 'a', color=sgcolor('white'), charsize=labchsz

    ;--poleward expansion.
    tpos = pos2
    tpos[1] = tpos[3]-(tpos[2]-tpos[0])/figxsz*figysz*!d.x_size/!d.y_size

    ; settings.
    tutr = time_double(['2013-06-07/04:54','2013-06-07/04:57'])
    velutr = time_double(['2013-06-07/04:54:54','2013-06-07/04:55:39'])

    ; check between these x-locations.
    txs = -smkarthm(2.4, 2.9, 0.1, 'dx')/(15*rad*scr0)
    ntx = n_elements(txs)
    idx = where(uts ge velutr[0] and uts le velutr[1])
    idx2 = median(idx)
    tut2 = uts[idx2]    ; middle time when u_up is calculated.

    ; draw snapshot.
    tmp = bytscl(reform(timg[idx2,*,*]), min = zr[0], max = zr[1], top = 254)>8
    sgtv, tmp, position = tpos, ct = 43, file = 'ct2'

    plot, xr, yr, /noerase, /nodata, position = tpos, $
        xstyle = 1, xticklen=xtklen, xtitle='MLT', $
        xminor=xminor, xtickv=xtickv, xticks=xticks, xtickname=xtickn, $
        ystyle = 1, yticklen=ytklen, ytitle = 'Z GSM (Re)', $
        yticks=yticks, yminor=yminor, ytickv=ytickv

    for i = 0, ntx-1 do plots, txs[i]+[0,0], yr, color = sgcolor('red')
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


    sgclose

end







;        tpos = [0.15,0.35,0.95,0.75]
;        sgclose
;            
;        ; plot the details.
;        ofn = shomedir()+'/fig_ewo_vel.pdf'
;        ;ofn = 0
;        
;        
;        vels = dblarr(nty)
;        
;        for i = 0, nty-1 do begin
;            ty = tys[i]
;            tidx = (ty-yr[0])/(yr[1]-yr[0])*figysz
;            tvar = 'ewo_'+string(i,format='(I0)')
;            txs = smkarthm(xr[0],xr[1],figxsz,'n')
;            idx = where(txs ge txr[0] and txs le txr[1])
;            ewos = reform(timg[*,idx,tidx])
;            txs = txs[idx]
;            
;            store_data, tvar, uts, ewos, txs, $
;                limits = {ytitle:'Azm.D from Midn (Re)', spec:1, no_interp:1, $
;                zrange:zr, yrange:txr, ystyle:1, zticks:3, ztickv:[1,2,3,4]*100, $
;                ztitle:'Photon Count', yticklen:-0.01}
;                
;            nouttick = (i eq nty-1)? 0: 1
;            tplot, tvar, position = poss[*,i], /noerase, nouttick = nouttick, /novtitle, trange = tutr
;            plot, tutr, txr, position = poss[*,i], /nodata, /noerase, $
;                xstyle=5, ystyle=5
;        
;        
;            verts = fltarr(nrec)
;            for j = 0, nrec-1 do begin
;                tdat = (reform(ewos[j,*]))   ; a slice at each time record.
;                idx = where(tdat ge vert0, cnt) ; find the highest point cross threshold.
;                verts[j] = (cnt eq 0)? !values.f_nan: txs[idx[0]]
;            endfor
;            idx = where(finite(verts))
;            verts = interpol(verts[idx],uts[idx],uts)
;            tvar = 'ewo_'+string(i,format='(I0)')+'_psbl_pos'
;            store_data, tvar, uts, verts, limits = {yrange:txr}
;        
;            pos0 = verts[where(uts eq velutr[0])]
;            pos1 = verts[where(uts eq velutr[1])]
;            vels[i] = (pos1-pos0)*re/(velutr[1]-velutr[0])  ; in km/s.
;            oplot, velutr, [pos0,pos1], color = 6, thick=8
;            xyouts, poss[0,i]+xchsz*1.5, poss[1,i]+ychsz*0.5, /normal, $
;                'Z GSM = '+sgnum2str(ty, nsgn=2)+' Re, V = '+sgnum2str(vels[i],ndec=0)+' km/s', $
;                color = 255
;                
;        endfor
;        
;        xyouts, (poss[0,0]+poss[2,0])*0.5, poss[3,0]+ychsz*0.5, /normal, $
;            alignment = 0.5, '"Ewogram" in tail at several Z GSMs, <V> west = '+ $
;            sgnum2str(abs(mean(vels)),ndec=0)+' km/s', charsize = 1.2
;            
;        sgclose
;        
;        
;        
;    
;    
;    stop
;
;
;
;; **** calc velocity of the PSBL.
;    figysz = 51     ; # of pixels in y.
;    figxsz = 301    ; # of pixels in x.
;    xr = -sctrng*15*rad*scr0 ; from MLT to distance in Re.
;    yr = sczrng
;    zr = [50,400]   ; color range for brightness.
;    vert0 = 120     ; threshold for vertical boundary.
;    figys = smkarthm(yr[0], yr[1], figysz, 'n')
;    
;    ; settings.
;    tutr = time_double(['2013-06-07/04:54','2013-06-07/04:57'])
;    velutr = time_double(['2013-06-07/04:54:54','2013-06-07/04:55:39'])
;    ; check between these x-locations.
;    txs = -smkarthm(2.4, 2.9, 0.1, 'dx')
;    ntx = n_elements(txs)
;    
;    ; load data.
;    tvar = tmodel+'_asf_vertical'
;    get_data, tvar, uts, timg
;    dr0 = sdatarate(uts)
;    
;
;
;    ; plot the positions along with aurora image.
;    ofn = shomedir()+'/fig_keo_vel_img.pdf'
;    ; ofn = 0
;    sgopen, ofn, xsize = figxsz*3, ysize = figysz*6
;    
;    xchsz = double(!d.x_ch_size)/!d.x_size
;    ychsz = double(!d.y_ch_size)/!d.y_size
;    
;    device, decomposed = 1
;    idx = where(uts ge velutr[0] and uts le velutr[1])
;    idx = median(idx)
;    tpos = [0.15,0.35,0.95,0.75]
;    tmp = bytscl(reform(timg[idx,*,*]), min = zr[0], max = zr[1], top = 254)>8
;    sgtv, tmp, position = tpos, ct = 43, file = 'ct2'
;    plot, xr, yr, /noerase, /nodata, position = tpos, $
;        xstyle = 1, ystyle = 1, xtitle = 'Azim. Dist (Re)', ytitle = 'Z GSM (Re)'
;    for i = 0, ntx-1 do plots, txs[i]+[0,0], yr, color = sgcolor('red')
;    xyouts, (tpos[0]+tpos[2])*0.5, tpos[3]+ychsz*0.5, alignment = 0.5, /normal, $
;        'Positions where KEO velcities are calculated'
;    xyouts, tpos[2]-xchsz*1, tpos[3]-ychsz*1.2, alignment = 1, /normal, $
;        time_string(uts[idx]), color = sgcolor('white')
;    sgclose
;
;    
;    
;    ; plot the details.
;    ofn = shomedir()+'/fig_keo_vel.pdf'
;    ;ofn = 0
;    
;    sgopen, ofn, xsize = 8.5, ysize = 11, /inch
;    device, decomposed = 0
;    loadct2, 43
;
;    poss = sgcalcpos(ntx)
;    
;    xchsz = double(!d.x_ch_size)/!d.x_size
;    ychsz = double(!d.y_ch_size)/!d.y_size
;    
;    vels = dblarr(ntx)
;    
;    for i = 0, ntx-1 do begin
;        tx = txs[i]
;        tidx = (tx-xr[0])/(xr[1]-xr[0])*301
;        keos = reform(timg[*,tidx,*])
;        tvar = 'keo_'+string(i,format='(I0)')
;        store_data, tvar, uts, keos, smkarthm(yr[0],yr[1],figysz,'n'), $
;            limits = {ytitle:'Z GSM (Re)', spec:1, no_interp:1, $
;            zrange:zr, yrange:yr, ystyle:1, zticks:3, ztickv:[1,2,3,4]*100, $
;            yticks:2, ytickv:[2,3,4], yminor:5, ztitle:'Photon Count', yticklen:-0.01}
;        
;        nouttick = (i eq ntx-1)? 0: 1
;        tplot, tvar, position = poss[*,i], /noerase, nouttick = nouttick, /novtitle, trange = tutr
;        plot, tutr, yr, position = poss[*,i], /nodata, /noerase, $
;            xstyle=5, ystyle=5
;        
;        verts = fltarr(nrec)
;        for j = 0, nrec-1 do begin
;            tdat = (reform(keos[j,*]))   ; a slice at each time record.
;            idx = where(tdat ge vert0, cnt) ; find the highest point cross threshold.
;            verts[j] = (cnt eq 0)? !values.f_nan: figys[idx[cnt-1]]
;        endfor
;        tvar = 'keo_'+string(i,format='(I0)')+'_psbl_pos'
;        store_data, tvar, uts, verts, limits = {yrange:yr}
;        
;        pos0 = verts[where(uts eq velutr[0])]
;        pos1 = verts[where(uts eq velutr[1])]
;        vels[i] = (pos1-pos0)*re/(velutr[1]-velutr[0])  ; in km/s.
;        plots, velutr, [pos0,pos1], color = 6, thick=8
;        xyouts, poss[0,i]+xchsz*1.5, poss[3,i]-ychsz*1.2, /normal, $
;            'MLT = '+sgnum2str(tx/scr0*deg/15+24, nsgn=4)+', V = '+sgnum2str(vels[i],ndec=0)+' km/s', $
;            color = 255
;    
;    endfor
;    
;    xyouts, (poss[0,0]+poss[2,0])*0.5, poss[3,0]+ychsz*0.5, /normal, $
;        alignment = 0.5, '"Keogram" in tail at several MLTs, <V> up = '+ $
;        sgnum2str(mean(vels),ndec=0)+' km/s', charsize = 1.2
;    
;    sgclose
;
;
;
;
;stop
;
;
;
;
;
;end
