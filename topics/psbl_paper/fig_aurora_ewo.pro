;---Constants.
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1

    test = 0



;---Settings.
    id = '2013_0607_event_info'
    if tnames(id) eq '' then load = 1 else load = 0
    if load eq 1 then _2013_0607_load_data3
    edge_elev = 7d  ; deg.

    get_data, id, tmp, einfo
    utr0 = time_double(['2013-06-07/04:54:30','2013-06-07/04:56:00'])
    utr1 = time_double(['2013-06-07/04:55:00','2013-06-07/04:55:30'])
    time_ticks = smkarthm(utr0[0], utr0[1], 30, 'dx')
    ofn = join_path([srootdir(),'fig_aurora_ewo.pdf'])
    if test eq 1 then ofn = 0


    site = einfo.sites
    pre0 = 'thg_'+site+'_asf_'
    mlonimg_var = pre0+'mlonimg'
    get_data, mlonimg_var, times, mlonimgs
    mlon_bins = get_var_data(pre0+'new_mlon_bins')
    mlat_bins = get_var_data(pre0+'new_mlat_bins')
    elev_2d = get_var_data(pre0+'new_elevs')
    edge_index = where(elev_2d le edge_elev)

    sgopen, ofn, xsize=5, ysize=4, /inch

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    thick = (size(ofn,/type) eq 7)? 8: 2
    ticklen = -0.02
    tplot_options, 'xticklen', ticklen
    tplot_options, 'yticklen', ticklen*0.5
    label_size = 0.8

;---Put some snapshots on.
    snapshots = time_double('2013-06-07/04:55:'+['00','09','18','27'])
    count_range = [0,1500]  ; the count range to scale the images.
    mlon_range = [-40,-25]
    mlat_range = [61,67]
    mlon_ticks = [-40,-35,-30,-25]
    mlat_ticks = [61,66]
    mlat_range_for_ewogram = [63,64.5]

    nsnapshot = n_elements(snapshots)
    poss = sgcalcpos(1,nsnapshot, region=[0,0,1,0.5], lmargin=2, rmargin=6, tmargin=2, bmargin=4, xpad=0.5)
    labels = ['a','b','c','d','e']
    for ii=0, nsnapshot-1 do begin
        index = where(times eq snapshots[ii])
        current_time = times[index]
        timg = reform(mlonimgs[index,*,*])
        timg = bytscl(timg, min=count_range[0], max=count_range[1], top=254)
        index = where(mlon_bins ge mlon_range[0] and mlon_bins le mlon_range[1])
        timg = timg[index,*]
        index = where(mlat_bins ge mlat_range[0] and mlat_bins le mlat_range[1])
        timg = timg[*,index]
        timg = transpose(timg)  ; make mlat the x-axis.
        timg = reverse(timg,1)

        tpos = poss[*,ii]
        ytitle = 'MLon (deg)'
        xtitle = 'MLat!C(deg)'
        xrange = reverse(mlat_range)
        yrange = mlon_range
        xtickv = float(mlat_ticks)
        ytickv = float(mlon_ticks)
        xticks = n_elements(xtickv)-1
        yticks = n_elements(ytickv)-1
        sgtv, timg, position=tpos, /resize, ct=1
        plot, mlat_bins, mlon_bins, /nodata, position=tpos, $
            /noerase, ticklen=ticklen, $
            xstyle=1, xrange=xrange, xminor=5, xtickformat='(A1)', xtickv=xtickv, xticks=xticks, xtitle='', $
            ystyle=1, yrange=yrange, yminor=5, ytickformat='(A1)', ytickv=ytickv, yticks=yticks, ytitle=''
        xyouts, tpos[0]+xchsz*1.5, tpos[1]+ychsz*0.2, /normal, orientation=90, $
            labels[ii]+'. '+time_string(current_time,tformat='hh:mm:ss UT'), color=sgcolor('white')
        txs = (xtickv-xrange[0])/(xrange[1]-xrange[0])*(tpos[2]-tpos[0])+tpos[0]
        tys = tpos[1]+fltarr(xticks+1)
        for jj=0, xticks do xyouts, txs[jj]+xchsz*0.5, tys[jj]-xchsz*3, /normal, sgnum2str(xtickv[jj],ndec=0), orientation=90, alignment=0
        xyouts, (tpos[0]+tpos[2])*0.5, tpos[1]-xchsz*5, /normal, alignment=0.5, xtitle, orientation=90
        foreach tval, mlat_range_for_ewogram do plots, tval+[0,0], yrange, color=sgcolor('white'), linestyle=1

        if ii eq 0 then begin
            foreach tval, mlat_range_for_ewogram do begin
                tmp = [tval, yrange[1]]
                tmp = convert_coord(tmp, /data, /to_normal)
                tx = tmp[0]+xchsz*0.4
                ty = tmp[1]+ychsz*0.2
                xyouts, tx, ty, /normal, alignment=0, orientation=90, sgnum2str(tval), charsize=label_size
            endforeach
        endif

        ; Add xticks.
        if ii eq nsnapshot-1 then begin
            txs = tpos[2]+fltarr(yticks+1)
            tys = (ytickv-yrange[0])/(yrange[1]-yrange[0])*(tpos[3]-tpos[1])+tpos[1]
            for jj=0, yticks do xyouts, txs[jj]+xchsz*2.4, tys[jj]+ychsz*(-0.2), /normal, alignment=0.5, sgnum2str(ytickv[jj],ndec=0), orientation=90
            xyouts, tpos[2]+xchsz*4.5, (tpos[1]+tpos[3])*0.5, /normal, alignment=0.5, ytitle, orientation=90
        endif
    endfor


;---Get the EWOgram.
    mlat_range = mlat_range_for_ewogram
    mlat_index = where(mlat_bins ge mlat_range[0] and mlat_bins le mlat_range[1], nmlat_index)
    ntime = n_elements(times)
    nmlon_bin = n_elements(mlon_bins)
    ewogram = fltarr(ntime,nmlon_bin)
    for ii=0, ntime-1 do begin
        for jj=0, nmlon_bin-1 do begin
            ewogram[ii,jj] = total(mlonimgs[ii,jj,mlat_index],3)/nmlat_index
            ewogram[ii,jj] = max(mlonimgs[ii,jj,mlat_index])
        endfor
    endfor

    ; Save it to tplot.
    yrange = mlon_range
    zrange = count_range
    store_data, pre0+'ewogram', times, ewogram, mlon_bins, limits={$
        spec:1, no_interp:1, $
        zrange:zrange, ztitle:'Photon Count', $
        yrange:yrange, ystyle:1, ytitle:'MLon (deg)', yticks:3, yminor:5}


;---Plot the EWOgram.
    tpos = sgcalcpos(1, region=[0,0.5,1,1], tmargin=1, bmargin=4, lmargin=10, rmargin=6)
    colorbar_pos = tpos & colorbar_pos[2] = tpos[0]-xchsz*1 & colorbar_pos[0] = colorbar_pos[2]-xchsz*0.8
    loadct, 1
    device, decomposed=0
    tplot_options, 'version', 2
    tvar = pre0+'ewogram'
    options, tvar, 'xtickformat', '(A1)'
    options, tvar, 'ytickformat', '(A1)'
    options, tvar, 'ztickformat', '(A1)'
    options, tvar, 'ytitle', ''
    options, tvar, 'ztitle', ''
    options, tvar, 'zposition', colorbar_pos
    tplot, pre0+'ewogram', trange=utr0, position=tpos, /novtitle, /noerase

    loadct2, 43
    plot, utr0, mlon_range, /nodata, /noerase, position=tpos, $
        xstyle=5, ystyle=5, xrange=utr0, yrange=mlon_range
    for ii=0, nsnapshot-1 do begin
        tmp = convert_coord(snapshots[ii],mlon_range[0], /data, /to_normal)
        plots, tmp[0]+[0,0], tmp[1]+ychsz*[-0.25,0.25], /normal, color=6, thick=thick
        xyouts, tmp[0]+xchsz*0.3, tmp[1]+ychsz*1, /normal, alignment=0.5, labels[ii], color=sgcolor('white'), orientation=90
    endfor

    ; Add x and y ticks.
    xrange = utr0
    xtickv = time_ticks
    xticks = n_elements(xtickv)-1
    xtick0 = xtickv[0]
    xtitle = 'Seconds from '+time_string(xtick0,tformat='YYYY-MM-DD/hh:mm:ss UT')
    txs = (xtickv-xrange[0])/(xrange[1]-xrange[0])*(tpos[2]-tpos[0])+tpos[0]
    tys = tpos[1]+fltarr(xticks+1)
    for jj=0, xticks do xyouts, txs[jj]+xchsz*0.5, tys[jj]-xchsz*2, /normal, sgnum2str(xtickv[jj]-xtick0,ndec=0), orientation=90, alignment=0.5
    xyouts, (tpos[0]+tpos[2])*0.5, tpos[1]-xchsz*6, /normal, alignment=0.5, xtitle, orientation=0

    txs = tpos[2]+fltarr(yticks+1)
    tys = (ytickv-yrange[0])/(yrange[1]-yrange[0])*(tpos[3]-tpos[1])+tpos[1]
    for jj=0, yticks do xyouts, txs[jj]+xchsz*2.4, tys[jj]+ychsz*(-0.2), /normal, alignment=0.5, sgnum2str(ytickv[jj],ndec=0), orientation=90
    xyouts, tpos[2]+xchsz*4.5, (tpos[1]+tpos[3])*0.5, /normal, alignment=0.5, ytitle, orientation=90

    ; Add the points and linear fit.
    xs = []
    ys = []
    for ii=0, ntime-1 do begin
        tdat = reform(mlonimgs[ii,*,mlat_index])
        index = where(tdat ge 1200, count)
        if count eq 0 then continue
        if times[ii] lt utr1[0] or times[ii] gt utr1[1] then continue
        tmp = array_indices(tdat, index)
        txs = times[ii]
        tys = mean(mlon_bins[tmp[0,*]])
        plots, txs+fltarr(count), mlon_bins[tmp[0,*]], psym=3, color=4
        plots, txs, tys, psym=6, color=6, symsize=0.5
        xs = [xs, txs]
        ys = [ys, tys]
    endfor

    res = linfit(xs, ys)
    txs = minmax(xs)+[-1,1]*6
    tys = res[0]+minmax(txs)*res[1]
    plots, txs, tys, color=6


    ; Add labels for the color bar.
    xyouts, colorbar_pos[0]-xchsz*1, colorbar_pos[1], /normal, orientation=90, alignment=0, sgnum2str(count_range[0],ndec=0), charsize=label_size
    xyouts, colorbar_pos[0]-xchsz*1, colorbar_pos[3], /normal, orientation=90, alignment=1, sgnum2str(count_range[1],ndec=0), charsize=label_size
    xyouts, colorbar_pos[0]-xchsz*2, (colorbar_pos[1]+colorbar_pos[3])*0.5, /normal, orientation=90, alignment=0.5, 'Photon Count', charsize=label_size

    vel = res[1]*60 ; deg/min.
    xyouts, xchsz*3, tpos[1], /normal, orientation=90, $
        labels[nsnapshot]+'. EWOgram '+sgnum2str(mlat_range[0])+'-'+sgnum2str(mlat_range[1])+' deg'
    xyouts, tpos[0]+xchsz*1.5, tpos[1]+ychsz*0.2, /normal, color=255, charsize=0.8, orientation=90, '|Lin.fit slope|='+sgnum2str(abs(vel),ndec=1)+' deg/min'

    sgclose
end
