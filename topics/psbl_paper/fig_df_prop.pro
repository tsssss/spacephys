;+
; Show the DF westward motion.
;-


test = 0

;---Find the event using results in the azim_df project.
    time_range = time_double(['2013-06-07/04:30','2013-06-07/05:30'])
    plot_time_range = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
    the_probes = ['g13','g15','rbspa','rbspb']
    local_root = srootdir()
    data_file = join_path([local_root,'data','fig_df_prop.tplot'])
    common_data_rate = 10.

    if file_test(data_file) eq 0 then begin
        vars = ['r_sm','theta']
        if n_elements(project) eq 0 then project = azim_df_load_project()

        save_vars = list()
        foreach probe, the_probes do begin
            prefix = probe+'_'
            data_time_range = mean(time_range)+[-1,1]*3600*5

            foreach var, vars do begin
                azim_df_read_data, var, time_range=data_time_range, probe=probe, project=project
                the_var = prefix+var
                uniform_time, the_var, common_data_rate
                save_vars.add, the_var
            endforeach

            mlt_var = prefix+'mlt'
            r_sms = get_var_data(prefix+'r_sm', times=times)
            mlt = azim_df_calc_pseudo_mlt(r_sms)
            store_data, mlt_var, times, mlt
            add_setting, mlt_var, /smart, {$
                unit: 'hr', $
                display_type: 'scalar', $
                short_name: 'MLT'}
            save_vars.add, mlt_var
        endforeach

        ; AE and Dst.
        omni_read_index, time_range
        save_vars.add, ['ae','dst'], /extract

        save_vars = save_vars.toarray()
        tplot_save, save_vars, filename=data_file
        lprmsg, 'Data saved to '+data_file+' ...'
    endif

    if file_test(data_file) then tplot_restore, filename=data_file


;---Determine the DF times.
    short_time_range = time_range
    sorted_probes = ['g13','rbspa','rbspb','g15']
    boxcar_window = 120.    ; sec.
    boxcar_ratio = 0.8
    dtheta_nsigma = 1
    df_min_duration = 60.   ; sec.
    df_list = list()
    foreach probe, sorted_probes do begin
        prefix = probe+'_'
        theta_var = prefix+'theta'
        azim_df_smooth_theta, theta_var, time_range, smooth_window=boxcar_window, stat_ratio=boxcar_ratio
        theta = get_var_data(theta_var, in=short_time_range, times=times)
        theta_smooth_combo_var = theta_var+'_smooth_combo'
        tmp = get_var_data(theta_smooth_combo_var)
        theta_smooth = tmp[*,0]
        theta_stddev = tmp[*,1]
        dtheta = tmp[*,2]
        dtheta_stddev = tmp[*,3]

        ; Pick out the times when dtheta has significant peaks.
        index = where(dtheta gt dtheta_stddev*dtheta_nsigma, count)
        if count eq 0 then begin
            lprmsg, tab+tab+tab+'No ramp with significant positive slope ...'
            continue
        endif
        time_ranges = time_to_range(times[index], time_step=common_data_rate)
        durations = time_ranges[*,1]-time_ranges[*,0]
        index = where(durations gt df_min_duration, ntime_range)
        if ntime_range eq 0 then begin
            lprmsg, tab+tab+tab+'No ramp lasts long enough ...'
            continue
        endif
        time_ranges = time_ranges[index,*]

        ; Pick out the ramps that have nodes.
        probe_df_list = list()
        for ii=0, ntime_range-1 do begin
            the_time_range = reform(time_ranges[ii,*])
            ; Exclude the time ranges that are on the edges.
            if the_time_range[0] eq time_range[0] then continue
            if the_time_range[1] eq time_range[1] then continue
            ; Find if there is a node.
            index = lazy_where(times, '[]', the_time_range, count=npoint)
            min_value_index = index[0]
            max_value_index = index[npoint-1]
            min_value = theta_smooth[min_value_index]
            max_value = theta_smooth[max_value_index]
            min_value_stddev = theta_stddev[min_value_index]
            max_value_stddev = theta_stddev[max_value_index]
            if min_value ge  min_value_stddev then continue
            if max_value le -max_value_stddev then continue
            if (max_value-min_value) le (min_value_stddev+max_value_stddev) then continue   ; value change is too small.
            the_value_range = [min_value, max_value]
                            arrival_time = mean(the_time_range)
            width = total(the_time_range*[-1,1])
            height = total(the_value_range*[-1,1])

            arrival_r_sm = get_var_data(prefix+'r_sm', at=arrival_time)
            df_info = dictionary($
                'probe', probe, $
                'time_range', the_time_range, $
                'value_range', the_value_range, $
                'arrival_time', arrival_time, $
                'arrival_r_sm', arrival_r_sm, $
                'width', width, $
                'height', height)
            probe_df_list.add, df_info
        endfor
        df_list.add, probe_df_list[0]   ; only 1 DF per probe.
    endforeach

;---Xcor.
    xcor_section_ratio = [-1.5,1.5]
    foreach the_df, df_list, jj do begin
        if jj eq 0 then begin
            xcor_max = 0.
            xcor_err = 0.
            xcor_time_lag = 0.
            xcor_time = the_df.arrival_time
        endif else begin
            pre_df = df_list[jj-1]
            pre_df_width = total(pre_df.time_range*[-1,1])
            the_df_width = total(the_df.time_range*[-1,1])
            the_wdith = max([pre_df_width,the_df_width])
            adjust_time_range = xcor_section_ratio*the_df_width
            pre_data = get_var_data(pre_df.probe+'_theta', in=pre_df.arrival_time+adjust_time_range, times=pre_time)
            the_data = get_var_data(the_df.probe+'_theta', in=the_df.arrival_time+adjust_time_range, times=the_time)
            adjust_time_lag = the_df.arrival_time-pre_df.arrival_time
            nlag = floor(the_df_width/common_data_rate)
            lags = findgen(nlag)-round(nlag/2)
            xcors = c_correlate(pre_data, the_data, lags)
            xcor_max = round(max(xcors, index)*10)/10.
            xcor_time_lag = lags[index]*common_data_rate
            if xcor_max gt 0 and finite(xcor_max) then begin
                xx = pre_data
                del_x = xx-mean(xx)
                dx_dt = deriv(xx)/common_data_rate
                nrec = n_elements(xx)
                xcor_err = sqrt(1./(nrec-1)*(1-xcor_max)/xcor_max*2*mean(del_x^2)/mean(dx_dt^2))
            endif else begin
                xcor_max = 2.
                xcor_err = 0.
                xcor_time_lag = 0.
            endelse
            xcor_time = pre_df.arrival_time+xcor_time_lag+adjust_time_lag
        endelse
        the_df['xcor_max'] = xcor_max
        the_df['xcor_err'] = xcor_err
        the_df['xcor_time_lag'] = xcor_time_lag
        the_df['xcor_time'] = xcor_time
    endforeach



;---Linear fit.
    hr_per_sec_2_deg_per_min = 15.*60
    df_times = dblarr(df_list.length)
    df_mlts = fltarr(df_list.length)
    foreach df, df_list, ii do begin
        df_times[ii] = 0.5*(df.arrival_time+df.xcor_time)
        df_mlts[ii] = azim_df_calc_pseudo_mlt(df.arrival_r_sm)
    endforeach
    xxs = df_times-df_times[0]
    yys = df_mlts
    fit_result = linfit(xxs,yys, yfit=yfit)
    linear_fit_slope = fit_result[1]*hr_per_sec_2_deg_per_min   ; in deg/min.
    ss_res = total((yys-yfit)^2)
    ss_tot = total((yys-mean(yys))^2)
    linear_fit_r_square = 1-ss_res/ss_tot



;---Plot settings.
    nprobe = n_elements(sorted_probes)
    probe_colors = smkarthm(100,250,nprobe, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)


    ; Cut off data outside the range.
    mlt_range = [-6,2]
    xy_range = [-5,10]
    if n_elements(project) eq 0 then project = azim_df_load_project()
    rxy_range = [4,30]
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
    ypans = [0.6,1]
    nypanel = n_elements(ypans)
    panel_ypads = [0.4,3]
    panel_x_ratio = 1.2     ; inch per hour.
    panel_xsize = abs(total(short_time_range*[-1,1]))/3600*panel_x_ratio
    ; Get the size of the right panel.
    right_pan_ysize = total(panel_ypads)*abs_ychsz+total(ypans)*panel_ysize
    right_pan_xsize = right_pan_ysize

    lower_left_pan_xsize = panel_xsize
    lower_left_pan_ysize = right_pan_ysize
    lower_right_pan_xsize = right_pan_xsize
    lower_right_pan_ysize = right_pan_ysize
    lower_panel_ypad = 5


    margins = [12.,5,8,2]
    panel_xpad = 12

    fig_ysize = right_pan_ysize+total(margins[[1,3]])*abs_ychsz+lower_right_pan_ysize+lower_panel_ypad*abs_ychsz
    fig_xsize = panel_xsize+panel_xpad*abs_xchsz+right_pan_xsize+total(margins[[0,2]])*abs_xchsz
    event_id = time_string(mean(time_range),tformat='YYYY_MMDD')
    file = join_path([srootdir(),'fig_df_prop.pdf'])
    if keyword_set(test) then file = test
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch

    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    poss = sgcalcpos(2, position=pos, ypad=lower_panel_ypad)
    upper_pos = poss[*,0]
    upper_poss = sgcalcpos(nypanel, position=upper_pos, ypans=ypans, ypad=panel_ypads, xchsz=xchsz, ychsz=ychsz)
    lower_pos = poss[*,1]
    poss = sgcalcpos(1,2, position=upper_pos, xpans=[panel_xsize,right_pan_xsize], xpad=panel_xpad)
    left_pos = poss[*,0]
    left_poss = sgcalcpos(nypanel, position=left_pos, ypans=ypans, ypad=panel_ypads, xchsz=xchsz, ychsz=ychsz)
    right_pos = poss[*,1]
    poss = sgcalcpos(1,2, position=lower_pos, xpans=[1.5,1], xpad=panel_xpad)
    lower_left_pos = poss[*,0]
    lower_left_poss = sgcalcpos(nypanel, position=lower_left_pos, ypad=panel_ypads)
    lower_right_pos =poss[*,1]


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
    tpos = upper_poss[*,0]
    xtickformat = '(A1)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    the_var = 'dst'
    ystep = 50.
    constants = [-50]
    yys = get_var_data(the_var, in=short_time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    ytitle = '(nT)'
    yminor = 5

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Add data.
    oplot, xxs, yys
    foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=9, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    the_var = 'ae'
    ystep = 1000.
    constants = [1000]
    yys = get_var_data(the_var, in=short_time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    ytitle = '(nT)'
    yminor = 5

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
    tpos = upper_poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = ''
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
    usersym, [1,1,-1,-1,1]*0.25, [-1,1,1,-1,-1]*1, /fill
    ct = 70
    reverse_ct = 1

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase

    pos_var_suffix = 'mlt'
    flags = list()
    foreach probe, sorted_probes, probe_id do begin
        prefix = probe+'_'
        pos_var = prefix+pos_var_suffix
        data_var = prefix+data_var_suffix
        get_data, data_var, xxs, zzs
        get_data, pos_var, xxs, yys
        if n_elements(zzs) ne n_elements(xxs) then begin
            errmsg = handle_error(data_var+': inconsistent data and time ...')
            stop
        endif
        mlt = get_var_data(prefix+'mlt')
        zzs = azim_df_normalize_theta(zzs, zrange=tilt_range, ct=ct, /reverse_ct)

        ; Remove data out of range.
        get_data, prefix+'r_sm', times, r_sm
        rxy = snorm(r_sm[*,0:1])
        ; Exclude data outside the MLT range.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = lazy_where(rxy, '][', rxy_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside time range.
        index = lazy_where(xxs, '][', short_time_range, count=count)
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
        
        tx = xxs[0]
        ty = yys[0]
        tmp = convert_coord(tx,ty,/data, /to_normal)
        tx = tpos[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.3
        if probe eq 'rbspb' then ty = tmp[1]-ychsz*0.9
        probe_info = resolve_probe(probe)
        short_name = probe_info.short_name
        xyouts, tx,ty,/normal, strupcase(short_name);, color=probe_colors[probe_id];, charsize=label_size;
    endforeach
    oplot, xrange, constant+[0,0], linestyle=constant_linestyle

    ; Linear-fit result.
    fit_times = short_time_range
    fit_mlts = (fit_times-df_times[0])*fit_result[1]+fit_result[0]
    oplot, fit_times, fit_mlts, linestyle=constant_linestyle

    squre_symbol = 6
    squre_size = label_size
    arrival_times = df_times
    arrival_mlts = df_mlts
    foreach probe, sorted_probes, ii do plots, arrival_times[ii],arrival_mlts[ii],/data, color=probe_colors[ii], psym=squre_symbol, symsize=squre_size
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1.2
    omega_fit = linear_fit_slope
    r_square = linear_fit_r_square
    msg = tex2str('omega')+'!Dfit!N = '+strmid(string(omega_fit,format='(F6.1)'),2)+' deg/min, r!U2!N = '+string(r_square,format='(F4.2)')
    xyouts, tx,ty,/normal, msg


    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


    ; Add bar.
    bar_thick = keyword_set(test)? 4: 8
    ;plots, plot_time_range, yrange[0]+[0,0], /data, thick=bar_thick
    tmp = convert_coord(plot_time_range[0], yrange[0], /data, /to_normal)
    ttpos = lower_left_poss[*,0]
    plots, [tmp[0],ttpos[0]], [tmp[1],ttpos[3]], /normal, linestyle=0
    tmp = convert_coord(plot_time_range[1], yrange[0], /data, /to_normal)
    plots, [tmp[0],ttpos[2]], [tmp[1],ttpos[3]], /normal, linestyle=0


    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'b. UT-'+strmid(ytitle,0,strpos(ytitle,' '))



;---Right panel.
    tpos = lower_right_pos
    tpos[2] = 1-2*xchsz
    xsize = (tpos[2]-tpos[0])*fig_xsize
    ysize = (tpos[3]-tpos[1])*fig_ysize
    the_size = min([xsize,ysize])
    xcenter = 0.5*([tpos[2]+tpos[0]])
    ycenter = 0.5*([tpos[3]+tpos[1]])
    tpos = [$
        xcenter-the_size/fig_xsize*0.5, $
        ycenter-the_size/fig_ysize*0.5, $
        xcenter+the_size/fig_xsize*0.5, $
        ycenter+the_size/fig_ysize*0.5]
    if n_elements(xy_range) ne 2 then xy_range = [5,-15]  ; sun to the left.
    panel_xrange = -xy_range
    panel_yrange = reverse(xy_range)

    dis_scale = 4.
    vel_scale = 40
    hsize = keyword_set(test)? 8: 160
    xtitle = 'SM X (Re)'
    ytitle = 'SM Y (Re)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xminor = 5
    xstep = 5
    xtickv = make_bins(panel_xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    yminor = 5
    ystep = 5
    ytickv = make_bins(panel_yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    plot, panel_xrange, panel_yrange, position=tpos, $
        xstyle=1, xrange=panel_xrange, xtitle=xtitle, xticklen=xticklen, xminor=xminor, xticks=xticks, xtickv=xtickv, $
        ystyle=1, yrange=panel_yrange, ytitle=ytitle, yticklen=yticklen, yminor=yminor, yticks=yticks, ytickv=ytickv, $
        /nodata, /noerase, /isotropic


    constant_mlts = [22,18,24,6,12]
    dis_range = [1,50]
    foreach constant_mlt, constant_mlts do begin
        angle = (constant_mlt*15+180)*constant('rad')
        oplot, dis_range*cos(angle), dis_range*sin(angle), linestyle=1
    endforeach
    the_mlt = 22
    the_dis = 8.5
    angle = (the_mlt*15+180)*constant('rad')
    tx = the_dis*cos(angle)
    ty = the_dis*sin(angle)
    tmp = convert_coord(tx,ty, /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    xyouts, tx,ty,/normal, string(the_mlt,format='(I0)')+' MLT', alignment=0.5


    ; Add arrow.
    dis_scale = 4.
    vel_scale = 20
    hsize = keyword_set(test)? 8: 160
    the_mlt = 22.
    the_dis = 5.
    angle = (the_mlt*15+180)*constant('rad')
    x0 = the_dis*cos(angle)
    y0 = the_dis*sin(angle)
    v_hat = [-sin(angle),cos(angle)]
    v_mag = omega_fit*constant('rad')/60*the_dis*constant('re')
    scale = dis_scale/vel_scale*v_mag
    x1 = x0+v_hat[0]*scale
    y1 = y0+v_hat[1]*scale
    arrow, x0,y0,x1,y1,/data, /solid, hsize=hsize
    tx = mean([x0,x1])
    ty = mean([y0,y1])
    tmp = convert_coord(tx,ty, /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    xyouts, tx,ty,/normal, sgnum2str(abs(v_mag),ndec=1)+' km/s', alignment=1;, charsize=constant('label_size')


    ; Add figure label.
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1.2
    label = 'e. SM x-y plane'
    xyouts, tx,ty,/normal,alignment=0, label

    ; Add earth.
    tmp = 50
    tmp = findgen(tmp)/(tmp-1)*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    polyfill, xs<0, ys, /line_fill, orientation=45, spacing=ychsz*6
    polyfill, xs<0, ys, color=sgcolor('silver')
    plots, xs, ys
    foreach r, make_bins(minmax(abs(xy_range)),5,/inner) do oplot, xs*r, ys*r, linestyle=1

    ; Add probes.
    foreach probe, sorted_probes, ii do begin
        rsm = df_list[ii].arrival_r_sm
        plots, rsm[0],rsm[1],/data, psym=squre_symbol, symsize=squre_size, color=probe_colors[ii]
        probe_info = resolve_probe(probe)
        short_name = probe_info.short_name
        tmp = convert_coord(rsm[0:1], /data, /to_normal)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]
        if short_name eq 'rbb' then ty = ty-ychsz*0.4
        xyouts, tx,ty,/normal, strupcase(short_name), color=probe_colors[ii]
    endforeach


;---Panel for RBSP.
    id = '2013_0607_event_info'
    if tnames(id) eq '' then _2013_0607_load_data3

    ; Common x-axis.
    xminor = 2
    xstep = 120. ; sec.
    xrange = plot_time_range
    xtickv = make_bins(xrange, xstep, /inner) ; make time line up.
    xlog = 0
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

    thick = keyword_set(test)? 2: 4
    theta_color = sgcolor('red')
    density_colors = [sgcolor('black'),sgcolor('gray')]
    density_thicks = [thick,1]

    theta_yrange = [-2,10]
    theta_yticks = 2
    theta_yminor = 3
    theta_ytitle = 'Detr. tilt '+tex2str('theta')+' (deg)'



    ; RBSP-A.
    tpos = lower_left_poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    get_data, 'rba_n_combine', times, density, limits=lim

    xtickformat = '(A1)'
    yrange = [0.01,1]
    ylog = 1
    yminor = 10
    ytickv = smkgmtrc(yrange[0],yrange[1],10,'dx')
    yticks = n_elements(ytickv)-1
    ytitle = 'Density (cm!U-3!N)'


    ; Draw the box.
    plot, xrange, yrange, $
        xstyle=5, xlog=xlog, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, ylog=ylog, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
        /nodata, /noerase, position=tpos

    for ii=0,1 do oplot, times, density[*,ii], color=density_colors[ii], thick=density_thicks[ii]

    ; Draw the box.
    plot, xrange, yrange, $
        xstyle=1, xlog=xlog, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=9, ylog=ylog, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase, position=tpos


    theta = get_var_data('rbspa_theta', in=xrange, times=times)
    axis, /yaxis, /save, yrange=theta_yrange, ylog=0, yticklen=yticklen, $
        ytitle=theta_ytitle, yticks=theta_yticks, yminor=theta_yminor, ystyle=1, color=theta_color
    oplot, times, theta, color=theta_color, thick=thick

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'c. RB-A'

    ; RBSP-B.
    tpos = lower_left_poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    get_data, 'rbb_n_combine', times, density, limits=lim

    xtickformat = ''
    yrange = [0.01,1]
    ylog = 1
    yminor = 10
    ytickv = smkgmtrc(yrange[0],yrange[1],10,'dx')
    yticks = n_elements(ytickv)-1
    ytitle = 'Density (cm!U-3!N)'


    ; Draw the box.
    plot, xrange, yrange, $
        xstyle=5, xlog=xlog, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, ylog=ylog, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
        /nodata, /noerase, position=tpos

    for ii=0,1 do oplot, times, density[*,ii], color=density_colors[ii], thick=density_thicks[ii]

    ; Draw the box.
    plot, xrange, yrange, $
        xstyle=1, xlog=xlog, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=9, ylog=ylog, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase, position=tpos


    theta = get_var_data('rbspb_theta', in=xrange, times=times)
    axis, /yaxis, /save, yrange=theta_yrange, ylog=0, yticklen=yticklen, $
        ytitle=theta_ytitle, yticks=theta_yticks, yminor=theta_yminor, ystyle=1, color=theta_color
    oplot, times, theta, color=theta_color, thick=thick

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'd. RB-B'



    if keyword_set(test) then stop
    sgclose

end
