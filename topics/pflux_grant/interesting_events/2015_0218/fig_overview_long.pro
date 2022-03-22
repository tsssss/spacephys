;+
; Plot over 9 h to show general plasma environment.
;-


;---Load data and settings.
    event_info = _2015_0218_02_load_data()
test = 0


    probe = event_info['probe']
    prefix = event_info['prefix']

    ps_times = time_double(['2015-02-18/00:02','2015-02-18/06:24'])
    event_time = time_double(['2015-02-18/02:10'])
    streamer_time = time_double(['2015-02-18/02:26'])
    onset_time = time_double(['2015-02-18/02:42'])
    cloak_time = time_double(['2015-02-18/00:40','2015-02-18/05:40'])

    long_time_range = time_double(['2015-02-17/22:35','2015-02-18/07:25'])


;---Settings for the plot.
    root_dir = srootdir()
    plot_file1 = join_path([root_dir,'fig_overview_long1.pdf'])
    plot_file2 = join_path([root_dir,'fig_overview_long2.pdf'])
    if keyword_set(test) then plot_file1 = 0
    if keyword_set(test) then plot_file2 = 1

    ps_color = sgcolor('black')
    cloak_color = ps_color

    rgb = constant('rgb')
    xyz = constant('xyz')
    label_size = 0.8
    hsize = keyword_set(test)? 6: 120
    xticklen_chsz = -0.2
    yticklen_chsz = -0.5


    ; Load en spec data.
    rbsp_read_en_spec, long_time_range, probe=probe
    rbsp_read_pa_spec, long_time_range, probe=probe

    ; Load density data.
    boom_pair = '12'
    efw_var = prefix+'efw_density'
    rbsp_efw_phasef_read_density, long_time_range, probe=probe, boom_pair=boom_pair, dmin=0
    rename_var, prefix+'density_12', to=efw_var
    options, efw_var, 'yrange', [1,1e4]
    options, efw_var, 'labflag', -1
    options, efw_var, 'labels', 'EFW Density'

    emfisis_var = prefix+'emfisis_density'
    rbsp_read_emfisis_density, long_time_range, probe=probe
    interp_time, emfisis_var, to=efw_var


    var = prefix+'density'
    rename_var, emfisis_var, to=var
    store_data, var, limits={$
        labels: ['Emfisis'], $
        ytitle: '(cm!U-3!N)', $
        labflag: -1, $
        yrange: [1,5e3], $
        ystyle: 1, $
        ylog: 1 }
    store_data, var, dlimit=0
;    var = prefix+'density'
;    stplot_merge, [efw_var,emfisis_var], newname=var
;    store_data, var, limits={$
;        labels: ['EFW','Emfisis'], $
;        colors: ['red','blue'], $
;        ytitle: '(cm!U-3!N)', $
;        labflag: -1, $
;        yrange: [1,5e3], $
;        ystyle: 1, $
;        ylog: 1 }


;    ; Load Vsc.
;    rbsp_efw_phasef_read_vsvy, long_time_range, probe=probe
;    vsc_var = prefix+'efw_vsvy'
;    vsvy = get_var_data(vsc_var, times=times)
;    vsc_var = prefix+'vsc'
;    vsc = total(vsvy[*,0:1],2)*0.5
;    store_data, vsc_var, times, vsc
;    add_setting, vsc_var, /smart, dictionary($
;        'display_type', 'scalar', $
;        'short_name', 'Vsc', $
;        'unit', 'V', $
;        'yrange', [-10,1], $
;        'ytickv', [0,-5,-10], $
;        'yticks', 2, $
;        'yminor', 5, $
;        'boom_pair', '12')


    ; Load AE.
    omni_read_index, long_time_range
    store_data, 'ae', limits={$
        ytitle: '(nT)', $
        labels: 'AE', $
        yrange: [0,1000], $
        ytickv: [0,500,1000], $
        yticks: 2, $
        yminor: 5 }


;---Prepare the plot.
    margins = [12,4,10,2]
    poss = panel_pos(nxpan=1, nypan=4, pansize=[4,0.5], ypans=[1,2,2,1.5], fig_size=fig_size1, margins=margins)

    vars = ['ae',prefix+['e_en_spec','p_en_spec','density']]
    nvar = n_elements(vars)

    options, prefix+'e_en_spec', 'ztickv', 10d^[5,7,9]
    options, prefix+'e_en_spec', 'zticks', 2
    options, prefix+'e_en_spec', 'ytitle', 'Energy (eV)'

    options, prefix+'p_en_spec', 'ztickv', 10d^[4,6,8]
    options, prefix+'p_en_spec', 'zticks', 2
    options, prefix+'p_en_spec', 'ytitle', 'Energy (eV)'

    sgopen, plot_file1, xsize=fig_size1[0], ysize=fig_size1[1], xchsz=xchsz, ychsz=ychsz

    tplot_options, 'num_lab_min', 4
    tplot_options, 'versions', 3
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    options, vars, 'xticklen', xticklen
    options, vars, 'yticklen', yticklen

    tplot_options, 'charsize', 1d
    tplot_options, 'xcharsize', 1d
    tplot_options, 'ycharsize', 1d
    tplot_options, 'zcharsize', 1d
    tplot_options, 'thick', 1d
    tplot_options, 'zticklen', 0


    tplot, vars, trange=long_time_range, position=poss

    fig_labels = letters(nvar)+') '+['AE','e-','H+','N']
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]-ychsz*0.7
        xyouts, tx,ty,/normal, fig_labels[ii]
    endfor

    ; Add PS and PP.
    tpos = poss[*,3]
    tpos[3] = poss[3,1]
    xrange = long_time_range
    yrange = [0,1]
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos
    txs = [ps_times]
    foreach tx, txs do begin
        plots, tx+[0,0], yrange, color=sgcolor('black')
    endforeach

    ; Add ps.
    tpos = poss[*,1]
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos
    x0 = mean(txs)
    tmp = convert_coord(x0,0, /data, /to_normal)
    x0 = tmp[0]
    y0 = tpos[3]-ychsz*3
    xyouts, x0,y0, /normal, alignment=0.5, 'Plasma sheet', color=ps_color
    foreach tx, txs, ii do begin
        tmp = convert_coord(tx,0, /data, /to_normal)
        x1 = tmp[0]
        y1 = y0+ychsz*0.3
        if ii eq 0 then begin
            tx = x0-xchsz*5
        endif else begin
            tx = x0+xchsz*5
        endelse
        arrow, tx,y1, x1,y1, /normal, solid=1, hsize=hsize, color=ps_color
    endforeach

    ; Add pp.
    tpos = poss[*,3]
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos
    x0 = mean(txs)
    y0 = 0.8
    tmp = convert_coord(x0,y0, /data, /to_normal)
    x0 = tmp[0]
    y0 = tmp[1]
    xyouts, x0,y0, /normal, alignment=0.5, 'Plasmapause'
    foreach tx, txs, ii do begin
        ty = 0.6
        tmp = convert_coord(tx,ty, /data, /to_normal)
        x1 = tmp[0]
        y1 = tmp[1]
        if ii eq 0 then begin
            tx = x0-xchsz*5
        endif else begin
            tx = x0+xchsz*5
        endelse
        arrow, tx,y0, x1,y1, /normal, solid=1, hsize=hsize
    endforeach


    ; Add cloak.
    cloak_times = [$
        mean([ps_times[0],cloak_time[0]]), $
        mean([ps_times[1],cloak_time[1]]), $
        event_time ]

    tpos = poss[*,2]
    xrange = long_time_range
    yrange = [0,1]
    plot, xrange, [0,1], $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos
    txs = [cloak_time]
    foreach tx, txs do begin
        plots, tx+[0,0], yrange, color=sgcolor('black')
    endforeach

    x0 = mean(txs)
    y0 = 0.8
    tmp = convert_coord(x0,y0, /data, /to_normal)
    x0 = tmp[0]-xchsz*1
    y0 = tmp[1]
    xyouts, x0,y0, /normal, alignment=0.5, 'Cloak', color=cloak_color
    foreach tx, cloak_times, ii do begin
        ty = 0.5
        tmp = convert_coord(tx,ty, /data, /to_normal)
        x1 = tmp[0]
        y1 = tmp[1]
        if ii eq 0 then begin
            tx = x0-xchsz*2
            ty = y0
        endif else if ii eq 2 then begin
            tx = x0
            ty = y0-ychsz*0.2
        endif else begin
            tx = x0+xchsz*2
            ty = y0
        endelse
        arrow, tx,ty, x1,y1, /normal, solid=1, hsize=hsize, color=cloak_color
    endforeach


    ; Add event times.
;    event_times = [event_time,streamer_time,onset_time]
;    event_labels = ['This event','Streamer','Approx. onset']

    event_times = [event_time]
    event_labels = ['This event: auroral beads & pseudo breakup']

    tpos = poss[*,0]
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos


    foreach tx, event_times do begin
        plots, tx+[0,0], yrange, color=sgcolor('black')
    endforeach


    dy = ychsz*0.2
    foreach tx, event_times, ii do begin
        tmp = convert_coord(tx,0, /data, /to_normal)
        x1 = tmp[0]
        y1 = tpos[3]
        plots, x1+[0,0], y1+[0,dy], /normal
;        if ii eq 0 then begin
;            x0 = x1-xchsz*5
;            dx = xchsz*3
;        endif else if ii eq 1 then begin
;            x0 = x1+xchsz*2
;            dx = -xchsz*1
;        endif else begin
;            x0 = x1+xchsz*9
;            dx = -xchsz*4
;        endelse
        x0 = x1
        dx = xchsz*2
        y0 = y1+ychsz*0.5
        plots, [x1,x0+dx], [y1+dy,y0], /normal
        msg = event_labels[ii]
        xyouts, x0+dx,y0+ychsz*0.2,/normal, msg, alignment=0.5, charsize=label_size
    endforeach

    if keyword_set(test) then stop
    sgclose







;---Plot sc location.
    fig_labels2 = (letters(nvar+2))[nvar:nvar+1]
    rbsp_read_orbit, long_time_range, probe=probe
    get_data, prefix+'r_gse', times, r_gse
    dis_var = prefix+'dis'
    store_data, dis_var, times, snorm(r_gse)
    add_setting, dis_var, /smart, dictionary($
        'short_name', 'RBSP-'+strupcase(probe), $
        'yrange', [1,6], $
        'ytickv', [2,4,6], $
        'yticks', 2, $
        'yminor', 2, $
        'unit', 'Re', $
        'display_type', 'scalar' )

    r_sm = cotran(r_gse, times, 'gse2sm')
    store_data, prefix+'r_sm', times, r_sm
    add_setting, prefix+'r_sm', /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'R', $
        'unit', 'Re', $
        'coord', 'SM', $
        'coord_labels', ['x','y','z'] )
;    r_mag = cotran(r_gse, times, 'gse2mag')
    mlts = atan(r_sm[*,1],r_sm[*,0])*constant('deg')/15+12

    ztickv = make_bins(long_time_range, 3600*2, /inner)
    zticks = n_elements(ztickv)-1
    ztickn = time_string(ztickv,tformat='hh:mm')
    ztickn[0] = time_string(ztickv[0],tformat='hh:mm!CMTH DD')

    xrange = [1.5,-6]
    yrange = [3.5,-2.2]
    xpan = abs(total(xrange*[-1,1]))
    ypan = abs(total(yrange*[-1,1]))

    poss = panel_pos(plot_file2, 'fit', nxpan=2, nypan=1, xpans=[1.5,1], pansize=[xpan,ypan], $
        xpad=8, margins=[6,4,1.5,0.5], xsize=fig_size1[0], fig_size=fig_size2)
    sgopen, plot_file2, xsize=fig_size2[0], ysize=fig_size2[1], xchsz=xchsz, ychsz=ychsz


    ; Set up coord.
    tpos = poss[*,0]
    plot, [0,1], [0,1], nodata=1, noerase=0, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, iso=1



    tmp = smkarthm(0,2*!dpi,20, 'n')
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs, tys, fill=1

    tmp = smkarthm(0,2*!dpi,50, 'n')
    txs = cos(tmp)
    tys = sin(tmp)


    ; Add earth and circles.
    oplot, xrange, [0,0], color=sgcolor('silver'), linestyle=1
    oplot, [0,0], yrange, color=sgcolor('silver'), linestyle=1
    polyfill, txs<0, tys, color=sgcolor('silver')
    polyfill, txs>0, tys, color=sgcolor('white')
    plots, txs, tys

    foreach rr, [2,3,4,5,6] do begin
        oplot, txs*rr, tys*rr, color=sgcolor('silver'), linestyle=1
    endforeach

    plots, r_sm[*,0], r_sm[*,1]

    ; Add time ticks.
    foreach ztick, ztickv, ii do begin
        tr = interp(r_sm, times, ztick)
        tmp = convert_coord(tr[0],tr[1], /data, /to_normal)
        tx = tmp[0]
        ty = tmp[1]
        plots, tx,ty,/normal, psym=8, symsize=0.5
        if ii eq 0 then begin
            xyouts, tx,ty+ychsz*1.8,ztickn[ii], /normal, alignment=0.5
        endif else begin
            xyouts, tx,ty-ychsz*1,ztickn[ii], /normal, alignment=0.5
        endelse
    endforeach

    ; Add event time.
    foreach ztick, event_time, ii do begin
        tr = interp(r_sm, times, ztick)
        tmp = convert_coord(tr[0],tr[1], /data, /to_normal)
        tx = tmp[0]
        ty = tmp[1]
        plots, tx,ty,/normal, psym=8, color=sgcolor('red'), symsize=0.8
;        msg = '('+strjoin(strtrim(string(tr,format='(F5.1)'),2),',')+') Re'
        mlt = interp(mlts, times, ztick)
        msg = 'MLT = '+string(mlt,format='(F4.1)')+' h'
        msg = msg+'!C|R| = '+strtrim(string(snorm(tr),format='(F5.1)'),2)+' Re'
        text_x1 = tpos[2]-xchsz*6
        text_y1 = tmp[1]+ychsz*1.5
        xyouts, text_x1,text_y1,/normal, alignment=0.5, msg, color=sgcolor('red')
    endforeach

    ; Add cloak times.
    text_x2 = text_x1-xchsz*5
    text_y2 = text_y1+ychsz*2.5
    x0 = text_x2
    y0 = text_y2
    xyouts, x0,y0, 'Outer boundary!Cof Cloak', /normal, alignment=0.5
    foreach ztick, cloak_time, ii do begin
        tr = interp(r_sm, times, ztick)
        tmp = convert_coord(tr[0],tr[1], /data, /to_normal)
        x1 = tmp[0]
        y1 = tmp[1]
        if ii eq 0 then begin
            tx = x0-xchsz*1
            ty = y0-ychsz*1.5
        endif else begin
            tx = x0;+xchsz*5
            ty = y0+ychsz*1
        endelse
        arrow, tx,ty,x1,y1, /normal, solid=1, hsize=hsize
    endforeach

    ; Add PS times.
    tr = sinterpol(r_sm, times, ps_times)
    x0 = mean([0,tr[*,0]])
    y0 = mean([0,tr[*,1]])
    tmp = convert_coord(x0,y0, /data, /to_normal)
    x0 = tmp[0]+xchsz*2
    y0 = tmp[1]-ychsz*0
    text_x3 = text_x2-xchsz*6
    text_y3 = text_y2+ychsz*2.5
    x0 = text_x3
    y0 = text_y3
    xyouts, x0,y0, 'Inner boundary of PS!C & Plasmapause', /normal, alignment=0.5
    foreach ztick, ps_times, ii do begin
        tr = interp(r_sm, times, ztick)
        tmp = convert_coord(tr[0],tr[1], /data, /to_normal)
        x1 = tmp[0]
        y1 = tmp[1]
        if ii eq 0 then begin
            tx = x0-xchsz*1
            ty = y0-ychsz*1.5
        endif else begin
            tx = x0;+xchsz*5
            ty = y0+ychsz*1
        endelse
        arrow, tx,ty,x1,y1, /normal, solid=1, hsize=hsize
    endforeach


    ; Add axes.
    xstep = 2
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = 4
    xtitle = 'SM X (Re)'

    ystep = 2
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 4
    ytitle = 'SM Y (Re)'

    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    plot, xrange, yrange, nodata=1, noerase=1, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, $
        ystyle=1, yrange=yrange, ytitle=ytitle, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, $
        position=tpos, iso=1

    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = fig_labels2[0]+') Orbit'
    xyouts, tx,ty,/normal, msg


;---ASI.
    sites = event_info['site']
    asi_time = time_double('2015-02-18/02:10:42')
    mlt_range = [-2.5,-0.5]
    mlat_range = [59,69]
    asi_ct = 1
    top_color = 254

    time_range = asi_time+[-1,1]*30
    themis_read_mltimg, time_range, sites=sites
    mltimg_var = 'thg_mltimg'
    get_data, mltimg_var, times, mltimgs
    time_index = where(times eq asi_time)
    mltimg = reform(mltimgs[time_index,*,*])
    mlt_bins = get_setting(mltimg_var, 'mlt_bins')
    mlt_index = lazy_where(mlt_bins, '[]', mlt_range)
    mlt_bins = mlt_bins[mlt_index]
    mltimg = mltimg[mlt_index,*]
    mlat_bins = get_setting(mltimg_var, 'mlat_bins')
    mlat_index = lazy_where(mlat_bins, '[]', mlat_range)
    mlat_bins = mlat_bins[mlat_index]
    mltimg = mltimg[*,mlat_index]

    tpos = poss[*,1]
    tpos[3] = tpos[3]-ychsz*2.5
    cbpos = tpos
    cbpos[1] = tpos[3]+ychsz*0.2
    cbpos[3] = cbpos[1]+ychsz*0.5

    zrange = [0,250]
    ztitle = 'Count (#)'
    zzs = bytscl(mltimg, min=zrange[0],max=zrange[1], top=top_color)
    sgtv, zzs, ct=asi_ct, position=tpos, resize=1
    sgcolorbar, findgen(top_color), horizontal=1, ztitle=ztitle, zrange=zrange, ct=asi_ct, position=cbpos

    xtitle = 'MLT (h)'
    xrange = mlt_range+24
    xstep = 0.5
    xtickv = make_bins(xrange,xstep,/inner)
    xticks = n_elements(xtickv)-1
    xminor = 5

    ytitle = 'MLat (deg)'
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange,ystep,/inner)
    yticks = n_elements(ytickv)-1
    yminor = 5

    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ; Add grid.
    plot, xrange, yrange, $
        xstyle=1, xlog=0, xrange=xrange, xtitle='', xtickformat='(A1)', xtickv=xtickv, xticks=xticks, xminor=1, xticklen=1, xgridstyle=1, $
        ystyle=1, ylog=0, yrange=yrange, ytitle='', ytickformat='(A1)', ytickv=ytickv, yticks=yticks, yminor=1, yticklen=1, ygridstyle=1, $
        position=tpos, nodata=1, noerase=1, ynozero=1, color=sgcolor('gray')

    ; Draw axes.
    plot, xrange, yrange, $
        xstyle=1, xlog=0, xrange=xrange, xtitle=xtitle, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, $
        ystyle=1, ylog=0, yrange=yrange, ytitle=ytitle, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
        position=tpos, nodata=1, noerase=1, ynozero=1

    ; Add labels.
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = fig_labels2[1]+') Auroral beads'
    xyouts, tx,ty,/normal, msg, color=sgcolor('white')

    ty = tpos[1]+ychsz*0.2
    msg = strupcase(sites[0])+' '+time_string(asi_time,tformat='hh:mm:ss')+' UT'
    xyouts, tx,ty,/normal, msg, color=sgcolor('white')

    model_setting = event_info['model_setting']
    model = model_setting['model']
    fmlt = get_var_data(prefix+'fmlt_'+model, at=asi_time)+24
    fmlat = get_var_data(prefix+'fmlat_'+model, at=asi_time)
    rb_color = sgcolor('red')
    tmp = convert_coord(fmlt, fmlat, /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    plots, tx,ty,/normal, psym=6, symsize=0.5, color=rb_color
    msg = 'RBSP-'+strupcase(probe)
    xyouts, tx,ty+ychsz*0.5,/normal, msg, color=rb_color, alignment=0.5

    ; Add wavelength.
    bead_color = sgcolor('white')
    bead_mlts = smkarthm(22.10,22.85,5,'n')
    bead_mlat = 61.8
    plots, minmax(bead_mlts), bead_mlat+[0,0], data=1, color=bead_color
    foreach bead_mlt, bead_mlts do begin
        tmp = convert_coord(bead_mlt,bead_mlat, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        plots, tx+[0,0], ty+[-1,1]*ychsz*0.25, normal=1, color=bead_color
    endforeach
    dtheta_b = abs(total(bead_mlts[0:1]*[-1,1]))
    msg = tex2str('theta')+'!Db!N = '+string(dtheta_b*15,format='(F3.1)')+' deg'
    tx = mean(minmax(bead_mlts))
    ty = bead_mlat
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]-ychsz*1.5
    xyouts, tx,ty,normal=1,alignment=0.5, msg, color=bead_color

    if keyword_set(test) then stop
    sgclose


end
