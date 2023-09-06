;+
; panel a is MLT-UT spectrogram.
; panel b are data after time lags removed.
; panel c is the xy-plane with arrows of good triads.
;-
pro azim_df_plot_timing_result, event_time_range, event_id=event_id, project=project, xy_range=xy_range

    test = 1


;---Check inputs.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(event_time_range) ne 2 then begin
        errmsg = handle_error('No input event_time_range ...')
        return
    endif

;---Find the event.
    events = azim_df_search_all_events(project=project)
    foreach event, events do if product(event.time_range-event_time_range) lt 0 then break
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
    nprobe = ndf
    foreach df, df_list, ii do probes[ii] = df.probe
    sorted_probes = probes
    obs_mlts = fltarr(ndf)
    foreach df, df_list, ii do obs_mlts[ii] = df.obs_mlt


;---Prepare colors.
    probe_colors = sgcolor(['firebrick','orange','sea_green','deep_sky_blue','medium_blue','dark_violet'])


;---Plot settings.
    ; Cut off data outside the range.
    region_name = event.region
    event_id = time_string(event.time_range[0],tformat='YYYY_MMDD_hh')
    if n_elements(mlt_range) ne 2 then mlt_range = (region_name eq 'post_midn')? [-3,9]: [-9,3]
    if n_elements(xsm_range) ne 2 then xsm_range = [-20,4]
    if n_elements(dis_range) ne 2 then dis_range = minmax(abs(xsm_range))
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
    panel_ysize = 1.5    ; inch.
    ypans = [1,0.6]
    nypanel = n_elements(ypans)
    panel_ypad = 3
    panel_x_ratio = 3.     ; inch per hour.
    panel_xsize = abs(total(event_time_range*[-1,1]))/3600*panel_x_ratio
    ; Get the size of the right panel.
    right_pan_ysize = panel_ypad*(nypanel-1)*abs_ychsz+total(ypans)*panel_ysize
    right_pan_xsize = right_pan_ysize

    margins = [12.,5,2,2]
    panel_xpad = 15
    fig_ysize = right_pan_ysize+total(margins[[1,3]])*abs_ychsz
    fig_xsize = panel_xsize+panel_xpad*abs_xchsz+right_pan_xsize+total(margins[[0,2]])*abs_xchsz
    file = join_path([project.plot_dir,'time_lag','fig_'+time_string(event.time_range[0],tformat='YYYY_MMDD_hh')+'_timing_result.pdf'])
    if keyword_set(test) then file = test
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch

    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    poss = sgcalcpos(1,2, position=pos, xpans=[panel_xsize,right_pan_xsize], xpad=panel_xpad)
    left_pos = poss[*,0]
    left_poss = sgcalcpos(nypanel, position=left_pos, ypans=ypans, ypad=panel_ypad, xchsz=xchsz, ychsz=ychsz)
    right_pos = poss[*,1]

    ;for ii=0, nypanel-1 do plot, [0,1], [0,1], position=left_poss[*,ii], /noerase, /nodata
    ;plot, [0,1], [0,1], position=right_pos, /noerase, /nodata

;---Panel a. EWOgram.
    xrange = event_time_range
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
    tpos = left_poss[*,0]
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
        get_data, data_var, xxs, zzs
        zzs = get_var_data(data_var, in=event_time_range, times=xxs)
        yys = get_var_data(pos_var, in=event_time_range, times=xxs)
        mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
        zzs = azim_df_normalize_theta(zzs, mlt, ct=ct, reverse_ct=reverse_ct, zrange=tilt_range)

        ; Remove data out of range.
        r_sm = get_var_data(prefix+'r_sm', in=event_time_range, times=times)
        xsm = r_sm[*,0]
        dis = snorm(r_sm)
        ; Exclude data outside the MLT range.
        index = where_pro(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = where_pro(dis, '][', dis_range, count=count)
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
    omega_fit = fit_result[1]*15*60

    txs = minmax(xxs)
    tys = txs*fit_result[1]+fit_result[0]
    oplot, txs, tys, linestyle=constant_linestyle
    txs = dblarr(nprobe)
    tys = fltarr(nprobe)
    foreach df, df_list, ii do begin
        txs[ii] = df.obs_time
        tys[ii] = df.obs_mlt
    endforeach
    squre_symbol = 6
    squre_size = label_size
    foreach probe, sorted_probes, ii do plots, txs[ii],tys[ii],/data, color=probe_colors[ii], psym=squre_symbol, symsize=squre_size
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1.2
    msg = tex2str('omega')+'!Dfit!N = '+strmid(string(omega_fit,format='(F6.1)'),2)+' deg/min, r!U2!N = '+string(r_square,format='(F4.2)')
    xyouts, tx,ty,/normal, msg


    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'a. UT-'+strmid(ytitle,0,strpos(ytitle,' '))

    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;---Time lag removed.
    tpos = left_poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    xtitle = ' '
    data_thick = 1
    arrow_hsize = keyword_set(test)? 4:80
    common_scale = 8.  ; deg.
    if event_id eq '2008_0204_12' then common_scale = 5

    data_var_suffix = 'theta'
    var_unit = 'deg'
    strlength = nprobe*3
    ty_label = tpos[3]+ychsz
    foreach df, df_list, ii do begin
        probe = df.probe
        prefix = probe+'_'
        time_lag = obs_times[ii]-obs_times[0]
        arrival_time = df.obs_time
        color = probe_colors[ii]

        get_data, prefix+data_var_suffix, xxs, yys
        xxs -= time_lag

        yrange = minmax(yys[index])
        yrange = (yrange[1]+yrange[0])*0.5+[-1,1]*(yrange[1]-yrange[0])*0.75

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
        xyouts, tx,ty_label,/normal, alignment=0.5, strupcase(short_name), color=color, charsize=label_size

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
    xyouts, tx,ty_label,/normal, 'b. '+tex2str('theta')+' time lag removed'


;---Right panel.
    tpos = right_pos
    if n_elements(xy_range) ne 2 then xy_range = [5,-15]  ; sun to the left.
    panel_xrange = xy_range
    panel_yrange = (region_name eq 'post_midn')? xy_range: -reverse(xy_range)
    if event_id eq '2008_0204_12' then begin
        panel_xrange = [5,-20]
        panel_yrange = [10,-15]
    endif
    dis_scale = 4.
    vel_scale = 40
    hsize = keyword_set(test)? 8: 160
    xtitle = 'SM X (Re)'
    ytitle = 'SM Y (Re)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    plot, panel_xrange, panel_yrange, xstyle=1, ystyle=1, position=tpos, $
        xrange=panel_xrange, yrange=panel_yrange, /nodata, /noerase, /isotropic, $
        xtitle=xtitle, ytitle=ytitle, xticklen=xticklen, yticklen=yticklen

    ; Add scale.
    tx = tpos[2]-xchsz*1
    ty = tpos[3]-ychsz*3.5
    tmp = convert_coord(tx,ty, /normal, /to_data)
    tmp = convert_coord(tmp[0]-dis_scale,tmp[1], /data, /to_normal)
    dx = abs(tmp[0]-tx)
    xs = tx+[-1,0]*dx
    plots, xs, ty+[0,0], /normal
    foreach tx, xs do plots, tx+[0,0], ty+[-1,1]*ychsz*0.1, /normal
    tx = mean(xs)
    ty = tpos[3]-ychsz*3
    xyouts, tx,ty,/normal, alignment=0.5, sgnum2str(vel_scale)+' km/s', charsize=label_size

    ; Add figure label.
    tx = tpos[2]-xchsz*1
    ty = tpos[3]-ychsz*1
    label = 'c. SM x-y plane'
    xyouts, tx,ty,/normal,alignment=1, label

    ; Add earth.
    tmp = 50
    tmp = findgen(tmp)/(tmp-1)*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    polyfill, xs<0, ys, /line_fill, orientation=45, spacing=ychsz*2
    plots, xs, ys
    foreach r, make_bins(minmax(abs(xy_range)),5,/inner) do oplot, xs*r, ys*r, linestyle=1


    triad_list = event.triad_list
    omegas = list()
    foreach triad, triad_list do begin
        rsm_center = triad.center_r_sm
        v_mag = triad.vmag_obs_time
        v_hat = triad.vhat_obs_time
        omegas.add, triad.omega_obs_time

        x0 = rsm_center[0]
        y0 = rsm_center[1]
        scale = dis_scale/vel_scale*v_mag
        x1 = x0+v_hat[0]*scale
        y1 = y0+v_hat[1]*scale
        arrow, x0,y0,x1,y1,/data, /solid, hsize=hsize
    endforeach
    omegas = omegas.toarray()
    omega_2d = mean(omegas)
    omega_2d_error = stddev(omegas)

    foreach df, df_list, ii do begin
        probe = df.probe
        rsm = df.obs_r_sm
        plots, rsm[0],rsm[1],/data, psym=squre_symbol, symsize=squre_size, color=probe_colors[ii]
        probe_info = resolve_probe(probe)
        short_name = probe_info.short_name
        tmp = convert_coord(rsm[0:1], /data, /to_normal)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]
        xyouts, tx,ty,/normal, strupcase(short_name), color=probe_colors[ii]
    endforeach

    tx = tpos[2]-xchsz*1
    ty = tpos[3]-ychsz*2
    msg = tex2str('omega')+'!D2D!N = '+strmid(string(omega_2d,format='(F6.1)'),2)+tex2str('pm')+strmid(string(omega_2d_error,format='(F6.1)'),2)+' deg/min'
    xyouts, tx,ty,/normal, alignment=1, msg

    if keyword_set(test) then stop
    sgclose
end

time_range = time_double(['2014-08-28/09:30','2014-08-28/11:30'])
xy_range = !null
;time_range = time_double(['2008-02-29/08:00','2008-02-29/09:30'])
;xy_range = [5,-20]
;time_range = time_double(['2008-02-04/14:30','2008-02-04/16:00'])
;xy_range = [5,-20]
;time_range = time_double(['2017-03-28/02:00','2017-03-28/05:00'])
azim_df_plot_timing_result, time_range, project=project, xy_range=xy_range
end
