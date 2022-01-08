;+
; Plot SC position, AE/Dst, and theta.
;
; events.
; project=.
; dirname=.
;-

pro azim_df_event_survey_plot2, project=project, event, dirname=dirname

test = 1
    secofhour = 3600.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(event) eq 0 then stop
    if n_elements(dirname) eq 0 then stop

    if isa(event,'list') then foreach tmp, event do azim_df_event_survey_plot1, project=project, dirname=dirname, tmp


    short_time_range = event.time_range+[-1,1]*secofhour*0.5
;    duration = total(event.time_range*[-1,1])
;    short_time_range = event.time_range+[-1,1]*duration
    df_list = event.df_list
    ndf = df_list.length
    probes = strarr(ndf)
    obs_times = dblarr(ndf)
    obs_mlts = fltarr(ndf)
    obs_rxys = fltarr(ndf)
    foreach df, event.df_list, ii do begin
        probes[ii] = df.probe
        obs_times[ii] = df.obs_time
        obs_mlts[ii] = df.obs_mlt
        obs_rxys[ii] = df.obs_rxy
    endforeach
    index = sort(obs_times)
    sorted_probes = probes[index]
    nprobe = n_elements(sorted_probes)
    probe_colors = smkarthm(100,250,nprobe, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)

;---Plot settings.
    ; Cut off data outside the range.
    region_name = event.region
    mlt_range = minmax(make_bins(make_bins(minmax(obs_mlts),1),3))
    xy_range = minmax(make_bins(minmax([-1,1,obs_rxys]),5))
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
    ypans = [0.6,1,1,0.6]
    nypanel = n_elements(ypans)
    panel_ypads = [0.4,0.4,3]
    panel_x_ratio = 3.0     ; inch per hour.
    panel_xsize = abs(total(short_time_range*[-1,1]))/3600*panel_x_ratio
    ; Get the size of the right panel.
    right_pan_ysize = total(panel_ypads)*abs_ychsz+total(ypans)*panel_ysize
    right_pan_xsize = right_pan_ysize


    margins = [12.,5,4,2]
    panel_xpad = 14

    fig_ysize = right_pan_ysize+total(margins[[1,3]])*abs_ychsz
    fig_xsize = panel_xsize+panel_xpad*abs_xchsz+right_pan_xsize+total(margins[[0,2]])*abs_xchsz
    event_id = time_string(mean(event.time_range),tformat='YYYY_MMDD')
    file = join_path([project.plot_dir,'time_lag','fig_'+event_id+'_timing_result_simple.pdf'])
    if keyword_set(test) then file = test
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch

    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    lower_pos = pos
    poss = sgcalcpos(1,2, position=lower_pos, xpans=[panel_xsize,right_pan_xsize], xpad=panel_xpad)
    left_pos = poss[*,0]
    left_poss = sgcalcpos(nypanel, position=left_pos, ypans=ypans, ypad=panel_ypads, xchsz=xchsz, ychsz=ychsz)
    right_pos = poss[*,1]
    tx = total(right_pos[[0,2]])*0.5
    ty = total(right_pos[[1,3]])*0.5
    dx = total(right_pos[[0,2]]*[-1,1])*0.4
    dy = total(right_pos[[1,3]]*[-1,1])*0.4
    right_pos = [tx-dx,ty-dy,tx+dx,ty+dy]


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
    foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

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


;---Color bar.
    tpos = left_poss[*,1]
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    tpos = left_poss[*,2]
    cbpos[1] = tpos[1]


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

    data_var_suffix = 'scaled_theta'
    ztitle = 'Scaled detr.tilt '+tex2str('theta')+' (deg)'
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
        zzs = get_var_data(data_var, times=xxs, in=short_time_range)
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

    square_symbol = 6
    cross_symbol = 1
    square_size = label_size
    foreach df, df_list, ii do plots, df.obs_time,df.obs_mlt,/data, color=probe_colors[ii], psym=square_symbol, symsize=square_size
    foreach df, df_list, ii do plots, df.xcor_time,df.obs_mlt,/data, color=probe_colors[ii], psym=cross_symbol, symsize=square_size


;    ; Linear-fit result.
;    hr_per_sec_2_deg_per_min = 15.*60
;    xxs = obs_times
;    yys = obs_mlts
;    fit_result = linfit(xxs, yys, yfit=yfit)
;    ss_res = total((yys-yfit)^2)
;    ss_tot = total((yys-mean(yys))^2)
;    r_square = 1-ss_res/ss_tot
;    omega_fit = fit_result[1]*hr_per_sec_2_deg_per_min
;
;    fit_times = short_time_range
;    fit_mlts = fit_times*fit_result[1]+fit_result[0]
;    oplot, fit_times, fit_mlts, linestyle=constant_linestyle
;    tx = tpos[0]+xchsz*1
;    ty = tpos[3]-ychsz*1.2
;    msg = tex2str('omega')+'!Dfit!N = '+strmid(string(omega_fit,format='(F6.1)'),2)+' deg/min, r!U2!N = '+string(r_square,format='(F4.2)')
;    xyouts, tx,ty,/normal, msg


    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


    ; Color bar.
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'b. UT-'+strmid(ytitle,0,strpos(ytitle,' '))




;---Panel b. UT-R short.
    tpos = left_poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    ytitle = 'Rxy (Re)'
    yminor = 5
    yrange = minmax(make_bins(obs_rxys,yminor))
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0

    data_var_suffix = 'scaled_theta'
    ztitle = 'Scaled detr.tilt '+tex2str('theta')+' (deg)'
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

    pos_var_suffix = 'r_sm'
    flags = list()
    foreach probe, sorted_probes do begin
        prefix = probe+'_'
        pos_var = prefix+pos_var_suffix
        data_var = prefix+data_var_suffix
        zzs = get_var_data(data_var, times=xxs, in=short_time_range)
        rsm = get_var_data(pos_var, at=xxs)
        yys = snorm(rsm[*,0:1])
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

    square_symbol = 6
    square_size = label_size
    foreach df, df_list, ii do plots, df.obs_time,df.obs_rxy,/data, color=probe_colors[ii], psym=square_symbol, symsize=square_size
    foreach df, df_list, ii do plots, df.xcor_time,df.obs_rxy,/data, color=probe_colors[ii], psym=cross_symbol, symsize=square_size


    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'c. UT-'+strmid(ytitle,0,strpos(ytitle,' '))





;---Time lag removed.
    tpos = left_poss[*,3]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    xtitle = ' '
    data_thick = 1
    arrow_hsize = keyword_set(test)? 4:80
    common_scale = 16.  ; deg.
    heights = []
    foreach df, df_list do heights = [heights,df.height]
    if min(heights) le 5 then common_scale *= 0.5

    data_var_suffix = 'theta'
    var_unit = 'deg'
    strlength = nprobe*3
    ty_label = tpos[3]+ychsz
    arrival_time0 = min(obs_times)
    foreach df, df_list, ii do begin
        probe = df.probe
        prefix = probe+'_'
        arrival_time = obs_times[ii]
        time_lag = arrival_time-arrival_time0
        color = probe_colors[ii]

        get_data, prefix+data_var_suffix, xxs, yys
        xxs -= time_lag

        index = lazy_where(xxs, '[]', short_time_range)
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
    xyouts, tx,ty_label,/normal, 'd. Detr.tilt '+tex2str('theta')+', time lag removed'


;---Right panel.
    tpos = right_pos
    if n_elements(xy_range) ne 2 then xy_range = [5,-15]  ; sun to the left.
    panel_xrange = -xy_range
    panel_yrange = (region_name eq 'post_midn')? -xy_range: reverse(xy_range)

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

    ; Add figure label.
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1.2
    label = 'e. SM x-y plane'
    xyouts, tx,ty,/normal,alignment=0, label

    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    plots, circle_x, circle_y
    foreach r, make_bins(minmax(abs(xy_range)),5,/inner) do oplot, circle_x*r, circle_y*r, linestyle=1


    ; Draw triad.
    triad_list = event.triad_list
    ntriad_vertex = 3
    triad_color = sgcolor('wheat')
    if triad_list.length le 3 then begin
        foreach triad, triad_list do begin
            r_sms = triad.r_sms[*,0:1]
            center_r_sm = triad.center_r_sm[0:1]
            for ii=0,ntriad_vertex-1 do begin
                txs = [center_r_sm[0],r_sms[ii,0]]
                tys = [center_r_sm[1],r_sms[ii,1]]
                plots, txs, tys, color=triad_color
            endfor
        endforeach
    endif


    ; Draw v_2d.
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

    ; Calculate mean omega_2d.
    if omegas.length lt 1 then begin
        omega_2d = 0.
        omega_2d_error = 0.
    endif else if omegas.length lt 2 then begin
        omega_2d = omegas[0]
        omeag_2d_error = 0.
    endif else begin
        omegas = omegas.toarray()
        omega_2d = abs(mean(omegas))
        omega_2d_error = stddev(omegas)
    endelse
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*2.4
    msg = '|'+tex2str('omega')+'!D2D!N| = '+strmid(string(omega_2d,format='(F6.1)'),2)+tex2str('pm')+strmid(string(omega_2d_error,format='(F6.1)'),2)+' deg/min'
    xyouts, tx,ty,/normal, alignment=0, msg


    ; Add SC as vertex.
    foreach df, df_list, ii do begin
        rsm = df.obs_r_sm
        probe = df.probe
        color = probe_colors[ii]
        plots, rsm[0],rsm[1],/data, psym=square_symbol, symsize=square_size, color=color
        probe_info = resolve_probe(probe)
        short_name = probe_info.short_name
        tmp = convert_coord(rsm[0:1], /data, /to_normal)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]
        xyouts, tx,ty,/normal, strupcase(short_name), color=color
    endforeach


    if keyword_set(test) then stop
    sgclose

end
