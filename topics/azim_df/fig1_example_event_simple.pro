;+
; Long time range:
;    panel a is Dst/AE.
;    panel b is MLT-UT spectrogram.
; Short time range:
;   panel c is MLT-UT spectrogram.
;   panel d is X-UT spectrogram.
;   panel e are data after time lags removed.
; Panel f is the xy-plane with arrows of good triads.
;-



pro fig1_example_event_simple, event, project=project, $
    short_time_range=short_time_range


test = 0

;---Check inputs.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(short_time_range) ne 2 then short_time_range = event.time_range+[-1,1]*1800.


    azim_df_load_basic_data, project=project


;---Sort based on obs_time.
    df_list = event.df_list
    ndf = df_list.length
    obs_times = dblarr(ndf)
    foreach df, df_list, ii do obs_times[ii] = df.obs_time
    index = sort(obs_times)
    obs_times = obs_times[index]
    df_list = df_list[index]
    probes = strarr(ndf)
    foreach df, df_list, ii do probes[ii] = df.probe
    sorted_probes = probes
    obs_mlts = fltarr(ndf)
    foreach df, df_list, ii do obs_mlts[ii] = df.obs_mlt
    obs_rxys = fltarr(ndf)
    foreach df, df_list, ii do obs_rxys[ii] = df.obs_rxy
    df_widths = fltarr(ndf)
    foreach df, df_list, ii do df_widths[ii] = df.width
    df_mean_width = mean(df_widths)



;---Prepare colors.
    nprobe = n_elements(sorted_probes)
    probe_colors = azim_df_probe_colors()


;---Plot settings.
    ; Cut off data outside the range.
    region_name = (strsplit(event.region,'%',/extract))[0]
    mlt_range = (region_name eq 'post_midn')? [-3,9]: [-9,3]
    probe_xys = [-1,1]
    foreach df, event.df_list do probe_xys = [probe_xys,df.obs_rxy]
    xy_range = minmax(make_bins(probe_xys,5))
    rxy_range = [4.,30]
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 2
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=xchsz, ychsz=ychsz
    sgclose, /wdelete
    abs_xchsz = xchsz
    abs_ychsz = ychsz
    ; Get the sizes of the left panels.
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1,0.6]
    nypanel = n_elements(ypans)
    panel_ypads = [0.4,3]
    panel_x_ratio = 3.0     ; inch per hour.
    panel_xsize = abs(total(short_time_range*[-1,1]))/3600*panel_x_ratio
    ; Get the size of the right panel.
    right_pan_ysize = total(panel_ypads)*abs_ychsz+total(ypans)*panel_ysize
    right_pan_xsize = right_pan_ysize


    margins = [12.,5,4,2]
    panel_xpad = 14

    fig_ysize = right_pan_ysize+total(margins[[1,3]])*abs_ychsz
    fig_xsize = panel_xsize+panel_xpad*abs_xchsz+right_pan_xsize+total(margins[[0,2]])*abs_xchsz
    event_id = time_string(mean(event.time_range),tformat='YYYY_MMDD_hh')
    file = join_path([project.plot_dir,'paper_events','fig_'+event_id+'_timing_result_simple.pdf'])
    if keyword_set(test) then file = test
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch

    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    lower_pos = pos
    poss = sgcalcpos(1,2, position=lower_pos, xpans=[panel_xsize,right_pan_xsize], xpad=panel_xpad)
    left_pos = poss[*,0]
    left_poss = sgcalcpos(nypanel, position=left_pos, ypans=ypans, ypad=panel_ypads, xchsz=xchsz, ychsz=ychsz)
    right_pos = poss[*,1]

;    for ii=0, upper_nypanel-1 do plot, [0,1], [0,1], position=upper_poss[*,ii], /noerase, /nodata
;    for ii=0, nypanel-1 do plot, [0,1], [0,1], position=left_poss[*,ii], /noerase, /nodata
;    plot, [0,1], [0,1], position=right_pos, /noerase, /nodata
;    stop

;---Common x-axis.
    xrange = short_time_range
    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    xminor = 10 ; min.
    xtickv = make_bins(xrange, 60*xminor, /inner) ; make time line up.
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    secofday = 86400d
    for ii=0, xticks do begin
        the_time = xtickv[ii]
        xtickn[ii] = time_string(the_time,tformat='hh:mm')
        date = time_string(the_time,tformat='YYYY-MM-DD')
        if ii eq 0 then begin
            xtickn[ii] += '!C'+date
            continue
        endif
        if the_time mod secofday ne 0 then continue
        xtickn[ii] += '!C'+date
    endfor



;---Panel a. Dst/AE.
    tpos = left_poss[*,0]
    xtickformat = '(A1)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    the_var = 'dst'
    ystep = 50.
    constants = [-50,0]
    yys = get_var_data(the_var, in=short_time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    ytitle = '(nT)'

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Add data.
    oplot, xxs, yys

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=9, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    the_var = 'ae'
    ystep = 500.
    constants = [0,500]
    yys = get_var_data(the_var, in=short_time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    ytitle = '(nT)'

    ; Set up coord.
    ae_color = sgcolor('red')
    axis, yaxis=1, /save, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        color=ae_color

    ; Add data.
    oplot, xxs, yys, color=ae_color
    foreach constant, constants do oplot, xrange, constant+[0,0], color=ae_color, linestyle=1

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'a. Dst/'
    xyouts, tx,ty,/normal, '           AE', color=ae_color



;---Panel b. UT-MLT short.
    tpos = left_poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    ytitle = 'MLT (hr)'
    yrange = mlt_range
    yminor = 3
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0

    data_var_suffix = 'theta'
    ztitle = 'Detrended tilt '+tex2str('theta')+' (deg)'
    tilt_range = [-1,1]*64
    zrange = tilt_range
    ztickv = [-64,-16,-4,0,4,16,64]
    ztickn = string(ztickv,format='(I0)')
    ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
    zrange = alog10(abs(zrange))*[-1,1]
    zticks = n_elements(ztickv)-1

    ; Symbol.
    psym = 8
    symsize = 0.5
    usersym, [1,1,-1,-1,1]*0.1, [-1,1,1,-1,-1]*1
    ct = 70
    reverse_ct = 1

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase

    pos_var_suffix = 'pseudo_mlt'
    flags = list()
    foreach probe, sorted_probes do begin
        prefix = probe+'_'
        pos_var = prefix+pos_var_suffix
        data_var = prefix+data_var_suffix
        zzs = get_var_data(data_var, in=short_time_range, times=xxs)
        yys = get_var_data(pos_var, at=xxs)
        mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
        zzs = azim_df_normalize_theta(zzs, zrange=tilt_range, ct=ct, /reverse_ct)

        ; Remove data out of range.
        r_sm = get_var_data(prefix+'r_sm', at=xxs)
        rxy = snorm(r_sm[*,0:1])
        ; Exclude data outside the MLT range.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = lazy_where(rxy, '][', rxy_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(r_sm)
        index = where(magn_flags eq 0, count)
        if count ne 0 then xxs[index] = !values.f_nan

        index = where(finite(xxs), count)
        the_flag = (count gt 0)
        flags.add, the_flag
        if count eq 0 then continue

        nxx = count
        xxs = xxs[index]
        yys = yys[index]
        zzs = zzs[index]
        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach
    oplot, xrange, constant+[0,0], linestyle=constant_linestyle

    ; Linear-fit result.
    xxs = obs_times
    yys = obs_mlts
    fit_result = linfit(xxs,yys, yfit=yfit)
    ss_res = total((yys-yfit)^2)
    ss_tot = total((yys-mean(yys))^2)
    r_square = 1-ss_res/ss_tot
    omega_fit = -fit_result[1]*15*60
    fit_times = short_time_range
    fit_mlts = fit_times*fit_result[1]+fit_result[0]
    oplot, fit_times, fit_mlts, linestyle=constant_linestyle


    arrival_times = obs_times
    arrival_mlts = obs_mlts
    foreach df, df_list, ii do plots, df.obs_time,df.obs_mlt,/data, psym=6, symsize=0.5
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1.2
    msg = tex2str('omega')+'!Dfit!N = '+strmid(string(omega_fit,format='(F6.1)'),2)+' deg/min, r!U2!N = '+string(r_square,format='(F4.2)')
    xyouts, tx,ty,/normal, msg


    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'b. UT-'+strmid(ytitle,0,strpos(ytitle,' '))


;---Time lag removed.
    tpos = left_poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    xtitle = ' '
    data_thick = 1
    arrow_hsize = keyword_set(test)? 4:80
    arrow_solid = 1

    data_var_suffix = 'theta'
    var_unit = 'deg'
    strlength = 3
    ty_label = tpos[3]+ychsz
    arrival_time0 = min(arrival_times)

    ; Set up coord.
    plot, xrange, [0,1], /noerase, /nodata, position=tpos, $
        xrange=xrange, xstyle=5, $
        yrange=yrange, ystyle=5

    ; Calculate position of sc short_name in normal coord.
    label_txs = fltarr(ndf)
    foreach df, df_list, ii do label_txs[ii] = (convert_coord(df.obs_time,0, /data, /to_normal))[0]
    min_dx = (strlength+0.5)*label_size*xchsz
    for ii=1,ndf-1 do begin
        if label_txs[ii]-label_txs[ii-1] ge min_dx then continue
        if ii eq 1 then begin
            label_txs[ii-1] = label_txs[ii]-min_dx
        endif else begin
            label_txs[ii] = label_txs[ii-1]+min_dx
        endelse
    endfor

    ; Calculate common_scale for height.
    yranges = fltarr(ndf,2)
    foreach df, df_list, ii do begin
        prefix = df.probe+'_'
        yys = get_var_data(prefix+data_var_suffix, in=df.time_range+[-1,1]*df.width*2.5, times=xxs)
        yrange = minmax(yys)
;        yrange = (yrange[1]+yrange[0])*0.5+[-1,1]*(yrange[1]-yrange[0])*0.75
        yrange = (yrange[1]+yrange[0])*0.5+[-1,1]*(yrange[1]-yrange[0])*1.1
        yranges[ii,*] = yrange
    endforeach
    ysteps = yranges[*,1]-yranges[*,0]
    min_ystep = min(ysteps)
    common_scale = ceil(min_ystep/2)*2

    foreach df, df_list, ii do begin
        probe = df.probe
        prefix = probe+'_'
        arrival_time = arrival_times[ii]
        time_lag = arrival_time-arrival_time0
        color = probe_colors[ii]
        yrange = reform(yranges[ii,*])

        yys = get_var_data(prefix+data_var_suffix, in=short_time_range, times=xxs)
        xxs -= time_lag

        ; Plot the data and print y-axis info.
        plot, xrange, yrange, /noerase, /nodata, position=tpos, $
            xrange=xrange, xstyle=5, $
            yrange=yrange, ystyle=5
        oplot, xxs, yys, color=color, thick=data_thick

        ; Add the common arrival_time.
        if ii eq 0 then oplot, arrival_time+[0,0], yrange, linestyle=1

        ; Add markers to indicate the time lag.
        ty = tpos[3]
        tmp = convert_coord(arrival_time,0, /data, /to_normal)
        tx = tmp[0]
        arrow, tx,ty+0.5*ychsz, tx,ty, /normal, color=color, hsize=arrow_hsize*2, /solid, thick=data_thick*2
        probe_info = resolve_probe(probe)
        short_name = probe_info.short_name
        plots, [tx,label_txs[ii]], $
            [ty+0.5*ychsz,ty_label-ychsz*0.2], /normal, color=color
        xyouts, label_txs[ii],ty_label,/normal, alignment=0.5, strupcase(short_name), charsize=label_size;, color=color

        ; Add y-labels for the yrange.
        tx = tpos[0]-xchsz*(ii+1)*0.75
        tmp = convert_coord(0,0, /data, /to_normal)
        y0 = tmp[1]
        tmp = convert_coord(0,common_scale, /data, /to_normal)
        y1 = tmp[1]
        dy = y1-y0
        tys = mean(tpos[[1,3]])+[-1,1]*0.5*dy
        plots, tx+[0,0], tys, /normal, color=color
        foreach ty,tys do plots, tx+[-1,1]*ychsz*0.1/fig_xsize*fig_ysize, ty, /normal, color=color

        if ii eq nprobe-1 then begin
            tx = tpos[0]-xchsz*(ii+2)*0.75
            ty = mean(tpos[[1,3]])-ychsz*0.3*label_size
            xyouts, tx,ty,/normal,alignment=1, sgnum2str(common_scale)+' '+var_unit, charsize=label_size
        endif
    endforeach

    axis, xaxis=0, xrange=xrange, xtitle=xtitle, xstyle=1, xminor=xminor, xticks=xticks, xtickv=xtickv, xtickname=xtickn, xticklen=xticklen

    tx = label_xshift*xchsz
    xyouts, tx,ty_label,/normal, 'c. Detrended '+tex2str('theta')+', time lag removed'


;---Right panel.
    tpos = right_pos
    panel_xrange = -xy_range
    panel_yrange = (region_name eq 'post_midn')? -xy_range: reverse(xy_range)

    dis_scale = 4.
    vel_scale = 40
    xtitle = 'SM X (Re)'
    ytitle = 'SM Y (Re)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    plot, panel_xrange, panel_yrange, xstyle=5, ystyle=5, position=tpos, $
        xrange=panel_xrange, yrange=panel_yrange, /nodata, /noerase, /isotropic, $
        xtitle=xtitle, ytitle=ytitle, xticklen=xticklen, yticklen=yticklen



    ; Add figure label.
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1.2
    label = 'd. SM x-y plane'
    xyouts, tx,ty,/normal,alignment=0, label


;---Draw shape.
    rad = constant('rad')
    deg = constant('deg')
    nangle = 50

    ; Gather info.
    triad_list = event.triad_list
    vmags = fltarr(triad_list.length)
    foreach triad, triad_list, ii do vmags[ii] = triad.vmag_obs_time
    mean_vmag = mean(vmags)

    omegas = fltarr(triad_list.length)
    foreach triad, triad_list, ii do omegas[ii] = triad.omega_obs_time
    omega_2d = mean(omegas)
    omega_2d_error = stddev(omegas)

    rxy_extent = minmax(obs_rxys)
    mlt_extent = minmax(obs_mlts)
    angle_extent = (mlt_extent*15+180)*rad
    if obs_mlts[0] gt obs_mlts[ndf-1] then angle_extent = reverse(angle_extent)
    mean_rxy = mean(rxy_extent)
    mean_angle = mean(angle_extent)

    ww = df_mean_width*mean_vmag/constant('re')
    ww_angle = abs(omega_2d)*min(df_widths)/60*rad


    ; Shape.
    mlt_extent_color = sgcolor('blue')
    shape_color = sgcolor('light_cyan')
    rr = 0
    the_angle = angle_extent[0]*(1-rr)+angle_extent[1]*rr
    angles = smkarthm(the_angle-ww_angle*0.5, the_angle+ww_angle*0.5,nangle, 'n')
    trs = smkarthm(rxy_extent[1],rxy_extent[0],nangle,'n')
    txs = [rxy_extent[1]*cos(angles),trs*cos(angles[-1]),rxy_extent[0]*reverse(cos(angles)),reverse(trs)*cos(angles[0])]
    tys = [rxy_extent[1]*sin(angles),trs*sin(angles[-1]),rxy_extent[0]*reverse(sin(angles)),reverse(trs)*sin(angles[0])]
    txs = min(panel_xrange)>txs<max(panel_xrange)
    tys = min(panel_yrange)>tys<max(panel_yrange)
    polyfill, txs, tys, /data, color=shape_color
    plots, txs,tys, /data, color=mlt_extent_color


    ; MLT extent.
    rr = 0.8
    the_rxy = rxy_extent[0]*(1-rr)+rxy_extent[1]*rr
    angles = smkarthm(angle_extent[0], angle_extent[1], nangle, 'n')
    txs = the_rxy*cos(angles)
    tys = the_rxy*sin(angles)
    plots, txs, tys, /data, color=mlt_extent_color
    arrow, txs[-2],tys[-2], txs[-1],tys[-1], /data, $
        hsize=arrow_hsize*1.5, solid=arrow_solid, color=mlt_extent_color
    tmp = convert_coord(0,ychsz, /normal,/to_data)-convert_coord(0,0, /normal,/to_data)
    drxy = tmp[1]*0.3
    tmp = the_rxy+[-1,1]*drxy
    plots, tmp*cos(angles[0]), tmp*sin(angles[0]), /data, color=mlt_extent_color

    tx = the_rxy*cos(angles[0])
    ty = the_rxy*sin(angles[0])
    tmp = convert_coord(tx,ty, /data,/to_normal)
    tx = tmp[0]+xchsz*1
    ty = tmp[1]
    ww_center = ww_angle*mean(rxy_extent)
    ww_10re = ww_angle*10
    xyouts, tx,ty,/normal, $
        'W = '+string(ww_angle*deg,format='(I0)')+' deg!C  ~ '+$
        sgnum2str(ww_10re,ndec=1)+' Re', charsize=label_size



;---Add earth.
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    plots, circle_x, circle_y
    foreach r, make_bins(minmax(abs(xy_range)),5,/inner) do oplot, circle_x*r, circle_y*r, linestyle=1



;---Draw v_2d and triads.
    ntriad_vertex = 3
    ntriad = triad_list.length
    triad_color = sgcolor('silver')
;    if ntriad le 3 then begin
;        foreach triad, triad_list do begin
;            r_sms = triad.r_sms
;            txs = r_sms[[0,1,2,0],0]
;            tys = r_sms[[0,1,2,0],1]
;            plots, txs, tys, color=triad_color
;        endforeach
;    endif

    foreach triad, triad_list do begin
        center_r_sm = triad.center_r_sm
        v_mag = triad.vmag_obs_time
        v_hat = triad.vhat_obs_time

        x0 = center_r_sm[0]
        y0 = center_r_sm[1]
        scale = dis_scale/vel_scale*v_mag
        x1 = x0+v_hat[0]*scale
        y1 = y0+v_hat[1]*scale
        arrow, x0,y0,x1,y1,/data, solid=arrow_solid, hsize=arrow_hsize*1.5
    endforeach


    square_symbol = 8
    ncirc = 21
    angles = findgen(ncirc)/(ncirc-1)*2*!dpi
    circ_x = cos(angles)
    circ_y = sin(angles)
    usersym, circ_x, circ_y, /fill
    square_size = label_size
    foreach df, df_list, ii do begin
        probe = df.probe
        rsm = df.obs_r_sm
        plots, rsm[0],rsm[1],/data, psym=square_symbol, symsize=square_size, color=probe_colors[ii]
        probe_info = resolve_probe(probe)
        short_name = probe_info.short_name
        tmp = convert_coord(rsm[0:1], /data, /to_normal)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]
        xyouts, tx,ty,/normal, strupcase(short_name);, color=probe_colors[ii]
    endforeach

    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*2.4
    msg = tex2str('omega')+'!D2D!N = '+strmid(string(omega_2d,format='(F6.1)'),2)+tex2str('pm')+strmid(string(omega_2d_error,format='(F6.1)'),2)+' deg/min'
    xyouts, tx,ty,/normal, alignment=0, msg


    ; Add scale.
    tx = tpos[2]-xchsz*2
    ty = tpos[3]-ychsz*2.5
    tmp = convert_coord(tx,ty, /normal, /to_data)
    tmp = convert_coord(tmp[0]-dis_scale,tmp[1], /data, /to_normal)
    dx = abs(tmp[0]-tx)
    xs = tx+[-1,0]*dx
    plots, xs, ty+[0,0], /normal
    foreach tx, xs do plots, tx+[0,0], ty+[-1,1]*ychsz*0.1, /normal
    tx = mean(xs)
    ty = tpos[3]-ychsz*2
    xyouts, tx,ty,/normal, alignment=0.5, sgnum2str(vel_scale)+' km/s'


    ; Add box.
    plot, panel_xrange, panel_yrange, xstyle=1, ystyle=1, position=tpos, $
        xrange=panel_xrange, yrange=panel_yrange, /nodata, /noerase, /isotropic, $
        xtitle=xtitle, ytitle=ytitle, xticklen=xticklen, yticklen=yticklen

    if keyword_set(test) then stop
    sgclose

end

;---Test events.
    test_event_times = time_double([$
        '2007-11-20/17:18:10', $
        '2008-01-07/08:37:20', $
        '2008-01-09/11:27:45', $
        '2008-02-29/08:26:50', $
        '2008-03-20/12:06:00', $
        '2009-04-07/07:05:15', $
        '2014-12-26/01:05:25', $
        '2016-10-13/12:22:35', $
        '2017-03-28/03:00:40'])


;---Load all df_group.
    project = azim_df_load_project()
    candidates = azim_df_find_dfgroup(project=project)


;---Select test_list.
    test_candidates = list()
    foreach test_event_time, test_event_times do begin
        lprmsg, 'Find event '+time_string(test_event_time)+' ...'
        foreach candidate, candidates do begin
            if candidate.time_range[0] eq test_event_time then begin
                lprmsg, 'Found the wanted event ...'
                test_candidates.add, candidate
                break
            endif
        endforeach
    endforeach


    foreach candidate, test_candidates do begin
        fig1_example_event_simple, candidate
    endforeach
end
