;+
; plot polar and fast ebratio spectrum, using wavelete result, color-code using altitude.
;-


;---Settings.
    rad = !dpi/180
    deg = 180/!dpi
    re = 6378d

    top = 254
    mindis = 1
    maxdis = 6
    psym = 1
    symsize = 0.7
    drec = 5    ; plot every 5 points.
    pres = ['po','fa']+'_'
    ct = 40
    
    ids = cusp_id('south_imf')
    ids = cusp_id('ion_outflow')
    
    
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then message, 'rootdir does not exist ...'
    
    datdir = rootdir+'/data/cusp'
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    
    
    ofn = 0
    ofn = shomedir()+'/cusp_fig_ebratio_spectrum.pdf'
    sgopen, ofn, xsize=8, ysize=4, /inch
    
    device, decomposed=0
    loadct, ct
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    

    tpos = [0.15,0.2,0.85,0.80]
    pos1 = tpos & pos1[1] = pos1[3]+ychsz*1 & pos1[3] = pos1[1]+ychsz*0.5
    poss = sgcalcpos(1,2, position=tpos, xpad=4)
    pos2 = poss[*,0]
    pos3 = poss[*,1]
    
    
    sgcolorbar, findgen(top), position=pos1, /horizontal, ztitle='R (Re)', $
        zrange=[mindis,maxdis], ct=ct, zcharsize=1
    

;---Spatial.
    xr = [1e-2,1e4]
    yr = [1e-1,1e5]
    tpos = pos3
    figlab = 'b.'
    
    tmp = alog10(xr)
    xtickv = smkarthm(tmp[0], tmp[1], 1, 'dx')
    xtickn = '10!U'+string(xtickv,format='(I0)')
    xticks = n_elements(xtickv)-1
    xtickv = 10d^xtickv
    xminor = 10
    xtickl = -0.02
    xtitle = 'Spatial Scale, (km)'
    
    tmp = alog10(yr)
    ytickv = smkarthm(tmp[0], tmp[1], 1, 'dx')
    ytickn = '10!U'+string(ytickv,format='(I0)')
    yticks = n_elements(ytickv)-1
    ytickv = 10d^ytickv
    yminor = 10
    ytickl = -0.02
    ytitle = 'E/B ratio (km/s)'

    plot, xr, yr, /nodata, position=tpos, $
        xstyle=1, xlog=1, xrange=xr, xtitle=xtitle, $
        xticks=xticks, xtickv=xtickv, xtickname=xtickn, xminor=xminor, xticklen=xtickl, $
        ystyle=1, ylog=1, yrange=yr, ytitle=ytitle, $
        yticks=yticks, ytickv=ytickv, ytickname=ytickn, yminor=yminor, yticklen=ytickl, $
        /noerase
    xyouts, tpos[0]-xchsz*7, tpos[3]-ychsz*0.8, /normal, figlab
    
    foreach tid, ids do begin        
        ; load all data.
        ifn1 = datdir+'/'+tid+'_1st_order_data.tplot'
        if file_test(ifn1) eq 0 then cusp_save_1st_order_data, tid, /save_data
        ifn2 = datdir+'/'+tid+'_2nd_order_data.tplot'
        if file_test(ifn2) eq 0 then cusp_save_2nd_order_data, tid, /save_data
        ifn3 = datdir+'/'+tid+'_3rd_order_data.tplot'
        if file_test(ifn3) eq 0 then cusp_save_3rd_order_data, tid, /save_data
        store_data, '*', /delete
        tplot_restore, filename=ifn1
        tplot_restore, filename=ifn2
        tplot_restore, filename=ifn3

        infofn = datdir+'/'+tid+'_scinfo.tplot'
        tplot_restore, filename=infofn
        

        get_data, 'scinfo', tutr, scinfo

        foreach pre0, pres do begin
            tinfo = (pre0 eq 'po_')? scinfo.polar: scinfo.fast
            dis0 = tinfo.dis0
            dilat = abs(tinfo.cusp_ilat[1]-tinfo.cusp_ilat[0])
            dtime = abs(tinfo.cusp_time[1]-tinfo.cusp_time[0])
            cc = dilat*rad*dis0*re/dtime    ; convert sec to km.
            
            get_data, pre0+'de_mor_fft_info', 0, einfo
            get_data, pre0+'db_mor_fft_info', 0, binfo
            ps = einfo.ps
            ebratio = sqrt(einfo.gws/binfo.gws)*1e3

            txs = ps*cc
            tys = ebratio
            txs = txs[0:*:drec]
            tys = tys[0:*:drec]
            
            tcolor = (dis0-mindis)/(maxdis-mindis)*top
            oplot, txs, tys, color=tcolor, psym=psym, symsize=symsize
        endforeach
    endforeach


;---Temporal.
    xr = [1e-2,1e4]
    yr = [1e-1,1e5]
    tpos = pos2
    figlab = 'a.'

    tmp = alog10(xr)
    xtickv = smkarthm(tmp[0], tmp[1], 1, 'dx')
    xtickn = '10!U'+string(xtickv,format='(I0)')
    xticks = n_elements(xtickv)-1
    xtickv = 10d^xtickv
    xminor = 10
    xtickl = -0.02
    xtitle = 'Temporal Scale, (sec)'
    
    tmp = alog10(yr)
    ytickv = smkarthm(tmp[0], tmp[1], 1, 'dx')
    ytickn = '10!U'+string(ytickv,format='(I0)')
    yticks = n_elements(ytickv)-1
    ytickv = 10d^ytickv
    yminor = 10
    ytickl = -0.02
    ytitle = 'E/B ratio (km/s)'

    plot, xr, yr, /nodata, position=tpos, $
        xstyle=1, xlog=1, xrange=xr, xtitle=xtitle, $
        xticks=xticks, xtickv=xtickv, xtickname=xtickn, xminor=xminor, xticklen=xtickl, $
        ystyle=1, ylog=1, yrange=yr, ytitle=ytitle, $
        yticks=yticks, ytickv=ytickv, ytickname=ytickn, yminor=yminor, yticklen=ytickl, $
        /noerase
    xyouts, tpos[0]-xchsz*7, tpos[3]-ychsz*0.8, /normal, figlab

        
    foreach tid, ids do begin
        ; load all data.
        ifn1 = datdir+'/'+tid+'_1st_order_data.tplot'
        if file_test(ifn1) eq 0 then cusp_save_1st_order_data, tid, /save_data
        ifn2 = datdir+'/'+tid+'_2nd_order_data.tplot'
        if file_test(ifn2) eq 0 then cusp_save_2nd_order_data, tid, /save_data
        ifn3 = datdir+'/'+tid+'_3rd_order_data.tplot'
        if file_test(ifn3) eq 0 then cusp_save_3rd_order_data, tid, /save_data
        store_data, '*', /delete
        tplot_restore, filename=ifn1
        tplot_restore, filename=ifn2
        tplot_restore, filename=ifn3
        
        infofn = datdir+'/'+tid+'_scinfo.tplot'
        tplot_restore, filename=infofn
        
        
        get_data, 'scinfo', tutr, scinfo
        
        foreach pre0, pres do begin
            tinfo = (pre0 eq 'po_')? scinfo.polar: scinfo.fast
            dis0 = tinfo.dis0
            dilat = abs(tinfo.cusp_ilat[1]-tinfo.cusp_ilat[0])
            dtime = abs(tinfo.cusp_time[1]-tinfo.cusp_time[0])
            cc = dilat*rad*dis0*re/dtime    ; convert sec to km.
            
            get_data, pre0+'de_mor_fft_info', 0, einfo
            get_data, pre0+'db_mor_fft_info', 0, binfo
            ps = einfo.ps
            ebratio = sqrt(einfo.gws/binfo.gws)*1e3
            
            txs = ps;*cc
            tys = ebratio
            txs = txs[0:*:drec]
            tys = tys[0:*:drec]
            
            tcolor = (dis0-mindis)/(maxdis-mindis)*top
            oplot, txs, tys, color=tcolor, psym=psym, symsize=symsize
        endforeach
    endforeach


    sgclose

end
