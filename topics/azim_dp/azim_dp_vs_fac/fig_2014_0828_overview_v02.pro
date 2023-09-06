;+
; DP lines, SC distribution on XY plane, DP EWOgram, and AE.
;
; Adopted from fig1_example_event, test_themis_flow_speed.
;-

;---Settings.
test = 0
    long_time_range = time_double(['2014-08-28/08:00','2014-08-28/24:00'])
    short_time_range = time_double(['2014-08-28/09:30','2014-08-28/11:30'])
    vel_time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
    plot_dir = join_path([homedir(),'Dropbox','mypapers','dp_vs_fac','plot'])
    plot_file = join_path([plot_dir,'fig_2014_0828_overview_v02.pdf'])

;---Settings for the start location and slopes (omega).
    the_line_style = 2
    the_line_color = sgcolor('dim_grey')
    the_start = dictionary($
        'times', time_double(['2014-08-28/10:08']), $
        'mlts', [0], $
        'omega1', [-2.1], $
        'omega2', [3] )
    pos_times = time_double('2014-08-28/'+['10:10','10:20'])


;---Find the event.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_find_dfgroup(project=project)
    foreach event, events do if product(event.time_range-short_time_range) lt 0 then break
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
    nprobe = ndf

;---Load ion velocity.
    root_dir = join_path([googledir(),'works','global_efield','data'])
    time_step = 1*60d
    common_times = make_bins(long_time_range, time_step)
    ntime = n_elements(common_times)
    ndim = 3
    foreach probe, ['a','d','e'], probe_id do begin
        prefix = 'th'+probe+'_'

        base = prefix+'ion_vel_'+time_string(long_time_range[0],tformat='YYYY')+'.tplot'
        file = join_path([root_dir,'ion_vel',base])
        tplot_restore, filename=file

        ; Convert data from gsm to sm.
        var = prefix+'u'
        in_var = var+'_gsm'
        out_var = var+'_sm'
        vec = get_var_data(in_var, times=times, in=long_time_range+[0,time_step])
        vec = cotran(vec, times, 'gsm2sm')
        store_data, out_var, times, vec
        lim = {$
            display_type: 'vector', $
            unit: 'km/s', $
            short_name: 'U!S!Uion!N!R', $
            coord: 'SM', $
            coord_labels: constant('xyz'), $
            colors: constant('rgb')}
        add_setting, out_var, /smart, lim

        vec_avg = fltarr(ntime,ndim)
        foreach time, common_times, time_id do begin
            time_index = where_pro(times, '[]', time+[0,time_step])
            for dim_id=0,ndim-1 do vec_avg[time_id,dim_id] = mean(vec[time_index,dim_id],/nan)
            ;for dim_id=0,ndim-1 do vec_avg[time_id,dim_id] = median(vec[time_index,dim_id])
        endforeach
        store_data, out_var, common_times, vec_avg
    endforeach


;---For mapping.
    model = 't01'
    map_time = time_double('2014-08-28/10:20')
    tilt = geopack_recalc(map_time)
    line_info = list()
    up_current = dictionary('lat_range', [60,67], 'color', sgcolor('salmon'))
    down_current = dictionary('lat_range', [68,72], 'color', sgcolor('dodger_blue'))
    rad = constant('rad')
    foreach lat, make_bins(up_current.lat_range,1)*rad do line_info.add, dictionary('lat',lat, 'color', up_current.color)
    foreach lat, make_bins(down_current.lat_range,1)*rad do line_info.add, dictionary('lat',lat, 'color', down_current.color)
    foreach lat, [50,60,70,80,90,-80,-90,100,110,120,130]*rad do line_info.add, dictionary('lat',lat, 'color', sgcolor('black'))
    nline = line_info.length



;---Plot settings.
    ; Cut off data outside the range.
    region_name = (strsplit(event.region,'%',/extract))[0]
    mlt_range = (region_name eq 'post_midn')? [-3,9]: [-9,3]
    probe_xys = [-1,1]
    foreach df, event.df_list do probe_xys = [probe_xys,df.obs_rxy]
;    xy_range = [-1.5,12]
    pos_xrange = [2,-15]
    pos_yrange = [2,-12]
    pos_zrange = [0,6]
    rxy_range = [4.,30]
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 2
    label_yshift = full_ychsz
    label_size = 0.8
    constant_linestyle = 1

    probe_colors = smkarthm(90,240,ndf, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)

    data_thick = 1
    arrow_hsize = keyword_set(test)? 4:80
    arrow_hsize *= 1.5
    arrow_solid = 1


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=xchsz, ychsz=ychsz
    sgclose, /wdelete
    abs_xchsz = xchsz
    abs_ychsz = ychsz

    ; Get the sizes of the long panels.
    panel_x_ratio = 3.     ; inch per hour.
    long_pan_xsize = abs(total(short_time_range*[-1,1]))/3600*panel_x_ratio
    long_ypans = [1,0.6]*0.8
    panel_ypad = 0.4
    nypanel = n_elements(long_ypans)
    long_pan_ysize = total(long_ypans)+(nypanel-1)*panel_ypad


    ; Get the sizes of right and left panels.
    margins = [12.,5,8,2]
    panel_xskip = 8
    panel_yskip = 4

    right_pan_xsize = long_pan_xsize*0.5
    right_pan_ysize = right_pan_xsize/abs(total(pos_xrange*[-1,1]))*(abs(total(pos_yrange*[-1,1]))+abs(total(pos_zrange*[-1,1])))+panel_ypad*abs_ychsz
    left_pan_xsize = long_pan_xsize-right_pan_xsize-panel_xskip*abs_xchsz
    left_pan_ysize = right_pan_ysize
    right_pan_shift = 5



    fig_ysize = right_pan_ysize+total(margins[[1,3]])*abs_ychsz+$
        panel_yskip*abs_ychsz+long_pan_ysize
    fig_xsize = long_pan_xsize+total(margins[[0,2]])*abs_xchsz
    event_id = time_string(mean(event.time_range),tformat='YYYY_MMDD_hhmm_ss')
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz


    pos = margins*[xchsz,ychsz,xchsz,ychsz]
    pos = [pos[0:1],1-pos[2:3]]
    tposs = sgcalcpos(2, position=pos, ypans=[right_pan_ysize,long_pan_ysize], ypad=panel_yskip)
    lower_pos = tposs[*,1]
    long_pan_poss = sgcalcpos(2, position=lower_pos, ypans=long_ypans, ypad=lower_panel_ypad)

    upper_pos = tposs[*,0]
    poss = sgcalcpos(1,2, position=upper_pos, xpans=[left_pan_xsize,right_pan_xsize], xpad=panel_xskip)
    left_pos = poss[*,0]
    right_pos = poss[*,1]
    left_pos[2] += right_pan_shift*xchsz
    right_pos[[0,2]] += right_pan_shift*xchsz

    right_poss = sgcalcpos(2,position=right_pos,ypad=panel_ypad, ypans=abs([total(pos_zrange*[-1,1]),total(pos_yrange*[-1,1])]))



;---Common x-axis.
    xrange = long_time_range
    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    xminor = 6 ; hr.
    xstep = 3600.
    xtickv = make_bins(xrange, xstep, /inner) ; make time line up.
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
    xtickn[1:*:2] = ' '





;---Panel a. Dst/AE.
    tpos = long_pan_poss[*,1]
    xtickformat = ''
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    the_var = 'ae'
    omni_read_index, long_time_range
    ystep = 500.
    constants = [500,1000]
    yminor = 5
    yys = get_var_data(the_var, in=long_time_range, times=xxs)
    yrange = minmax(make_bins([0,minmax(yys)],100))
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    ytitle = '(nT)'

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Add data.
    oplot, xxs, yys, color=ae_color
    foreach constant, constants do oplot, xrange, constant+[0,0], color=ae_color, linestyle=1

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'e. AE'


;---Panel b. UT-MLT long.
    tpos = long_pan_poss[*,0]
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
    ztitle = 'Scaled '+tex2str('theta')+' (deg)'
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

    time_step = 10.
    common_times = make_bins(long_time_range+[0,-1]*time_step, time_step)
    ntime = n_elements(common_times)
    thetas = fltarr(ntime,nprobe)
    mlts = fltarr(ntime,nprobe)
    rxys = fltarr(ntime,nprobe)
    rsms = fltarr(ntime,3,nprobe)
    foreach probe, probes, probe_id do begin
        prefix = probe+'_'
        the_var = prefix+'scaled_theta'
        if tnames(the_var) eq '' then continue
        thetas[*,probe_id] = get_var_data(the_var,at=common_times)
        mlts[*,probe_id] = get_var_data(prefix+'pseudo_mlt',at=common_times)
        rsms[*,*,probe_id] = get_var_data(prefix+'r_sm',at=common_times)
        rxys[*,probe_id] = snorm(reform(rsms[*,0:1,probe_id]))
    endforeach


    foreach time, common_times, time_id do begin
        zzs = thetas[time_id,*]
        yys = mlts[time_id,*]
        rxy = rxys[time_id,*]
        rsm = rsms[time_id,*,*]

        index = where(finite(zzs),count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = azim_df_normalize_theta(zzs[index], zrange=tilt_range, ct=ct, /reverse_ct)
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Filter spatially.
        index = where_pro(yys, '[]', yrange, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = where_pro(rxy, '[]', rxy_range, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = where(check_if_in_magn(rsm, dynamic_pressure=pdyn) eq 1, count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Sort by tilt.
        index = reverse(sort(zzs))
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        for ii=0, count-1 do plots, time, yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach

    oplot, xrange, constant+[0,0], linestyle=constant_linestyle

    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase
    bar_thick = (size(file,/type) eq 7)? 8: 4
    plots, short_time_range, yrange[1]+[0,0], color=sgcolor('black'), thick=bar_thick
    tmp = convert_coord(short_time_range[0], yrange[1], /data, /to_normal)
    ttpos = left_pos[*,0]
    plots, [tmp[0],ttpos[0]], [tmp[1],ttpos[1]], /normal, linestyle=0
    tmp = convert_coord(short_time_range[1], yrange[1], /data, /to_normal)
    ttpos = left_pos[*,0]
    plots, [tmp[0],ttpos[2]], [tmp[1],ttpos[1]], /normal, linestyle=0

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'd. Scaled '+tex2str('theta')


    ; Add "Drift" examples.
    ypos = tpos[1]+ychsz*0.5
    xpos = tpos[0]+xchsz*0.5
    drift_label = '"Drift" in MLT'
    xyouts, xpos,ypos,/normal, drift_label;, charsize=label_size
    drift_times = list()
    drift_times.add, ['10:05','10:45']
    drift_times.add, ['12:35','13:10']
    drift_times.add, ['15:40','16:20']
    drift_times.add, ['16:35','17:10']
    drift_times.add, ['18:40','19:15']
    drift_times.add, ['20:15','21:15']
    date = time_string(short_time_range, tformat='YYYY-MM-DD')
    data_ypos = (convert_coord(xpos,ypos+ychsz*0.8, /normal, /to_data))[1]
    label_ypos = (convert_coord(xpos,ypos, /normal, /to_data))[1]
    foreach drift_time, drift_times, drift_id do begin
        drift_time_range = time_double(date+'/'+drift_time)
        plots, drift_time_range, data_ypos
        xyouts, mean(drift_time_range), label_ypos, $
            /data, alignment=0.5, string(drift_id+1,format='(I0)');, charsize=label_size
    endforeach

    ; Add omega.
    txs = the_start.times
    tys = the_start.mlts
    omega1 = the_start.omega1
    nstart = n_elements(txs)
    for ii=0,nstart-1 do begin
        tmp1 = [txs[ii],short_time_range[1]]
        the_omega = omega1[ii]
        if finite(the_omega) then begin
            slope = the_omega/60./15
            tmp2 = tys[ii]-[0,slope*(tmp1[1]-tmp1[0])]
            oplot, tmp1, tmp2, linestyle=the_line_style, color=the_line_color

            ; Find the position to print omega.
            ty = tpos[3]-ychsz*1.5
            tmp = convert_coord(tpos[0],ty, /normal, /to_data)
            tx = -(tmp[1]-tys[ii])/slope+txs[ii]
            tmp = convert_coord(tx,tmp[1], /data, /to_normal)
            tx = tmp[0]-xchsz*1
            ty = tmp[1]
            msg = tex2str('omega')+' = '+strtrim(string(the_omega,format='(F5.1)'))+' deg/min'
            xyouts, tx,ty,/normal, alignment=1, charsize=label_size, msg

        endif
    endfor

    ; Add color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;---Right panels.
    tpos = right_poss[*,1]
    panel_xrange = pos_xrange
    panel_yrange = pos_yrange

    dis_scale = 4.
    vel_scale = 40
    xtitle = 'SM X (Re)'
    ytitle = 'SM Y (Re)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    plot, panel_xrange, panel_yrange, xstyle=1, ystyle=1, position=tpos, $
        xrange=panel_xrange, yrange=panel_yrange, /nodata, /noerase, /isotropic, $
        xtitle=xtitle, ytitle=ytitle, xticklen=xticklen, yticklen=yticklen

    ; Add lines and circles.
    label_color = sgcolor('silver')
    oplot, panel_xrange, [0,0], linestyle=1, color=label_color
    oplot, [0,0], panel_yrange, linestyle=1, color=label_color

    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    foreach r, make_bins(minmax(abs(pos_xrange)),5,/inner) do oplot, circle_x*r, circle_y*r, linestyle=1, color=label_color


    ; Add figure label.
    tx = tpos[2]-xchsz*1
    ty = tpos[1]+ychsz*0.5
    label = 'c. Equatorial (XY) plane'
    xyouts, tx,ty,/normal,alignment=1, label


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
    ww_angle = abs(omega_2d)*min(df_widths)/60*rad
    ww_angle = 15*constant('rad')


    ; Overall extent.
    mlt_extent_color = sgcolor('blue')
    angles = smkarthm(angle_extent[0]-ww_angle*0.5, angle_extent[1]+ww_angle*0.5,nangle, 'n')
    trs = smkarthm(rxy_extent[1],rxy_extent[0],nangle,'n')
    txs = [rxy_extent[1]*cos(angles),trs*cos(angles[-1]),rxy_extent[0]*reverse(cos(angles)),reverse(trs)*cos(angles[0])]
    tys = [rxy_extent[1]*sin(angles),trs*sin(angles[-1]),rxy_extent[0]*reverse(sin(angles)),reverse(trs)*sin(angles[0])]
    txs = min(panel_xrange)>txs<max(panel_xrange)
    tys = min(panel_yrange)>tys<max(panel_yrange)
    ;polyfill, txs, tys, /data, color=shape_color
    plots, txs,tys, /data, color=mlt_extent_color, linestyle=2

    ; Shape.
    shape_color = sgcolor('purple')
    the_angle = angle_extent[0]+total(angle_extent*[-1,1])*0.6
    angles = smkarthm(the_angle-ww_angle*0.5, the_angle+ww_angle*0.5,nangle, 'n')
    trs = smkarthm(rxy_extent[1],rxy_extent[0],nangle,'n')
    txs = [rxy_extent[1]*cos(angles),trs*cos(angles[-1]),rxy_extent[0]*reverse(cos(angles)),reverse(trs)*cos(angles[0])]
    tys = [rxy_extent[1]*sin(angles),trs*sin(angles[-1]),rxy_extent[0]*reverse(sin(angles)),reverse(trs)*sin(angles[0])]
    txs = min(panel_xrange)>txs<max(panel_xrange)
    tys = min(panel_yrange)>tys<max(panel_yrange)
    plots, txs,tys, /data, color=shape_color, linestyle=2
    rr = 8
    tx = rr*cos(the_angle)
    ty = rr*sin(the_angle)
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]+ychsz*0.2
    msg = 'DP!Cwidth'
    xyouts, tx,ty,normal=1, alignment=0.5, msg, color=shape_color, charsize=label_size


    ; MLT extent.
    rr = 0.8
    the_rxy = rxy_extent[0]*(1-rr)+rxy_extent[1]*rr
    angles = smkarthm(angle_extent[0], angle_extent[1], nangle, 'n')
    txs = the_rxy*cos(angles)
    tys = the_rxy*sin(angles)
    plots, txs, tys, /data, color=mlt_extent_color
    arrow, txs[-2],tys[-2], txs[-1],tys[-1], /data, $
        hsize=arrow_hsize, solid=arrow_solid, color=mlt_extent_color
    tmp = convert_coord(0,ychsz, /normal,/to_data)-convert_coord(0,0, /normal,/to_data)
    drxy = tmp[1]*0.3
    tmp = the_rxy+[-1,1]*drxy
    plots, tmp*cos(angles[0]), tmp*sin(angles[0]), /data, color=mlt_extent_color

    tx = the_rxy*cos(angles[-1])
    ty = the_rxy*sin(angles[-1])
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]-ychsz*1.2
    msg = 'Lower est.!Cpropagation!Cextent'
    xyouts, tx,ty,normal=1, msg, color=mlt_extent_color, charsize=label_size


;---Add earth.
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    polyfill, circle_x>0, circle_y, color=sgcolor('white')
    plots, circle_x, circle_y


;---Draw v_2d and triads.
    ntriad_vertex = 3
    triad_list = event.triad_list
    the_triad_probes = ['g15','thd','the']
    the_key = strjoin(sort_uniq(the_triad_probes),'_')
    triad_color = sgcolor('gray')
    foreach triad, triad_list do begin
        center_r_sm = triad.center_r_sm
        v_mag = triad.vmag_obs_time
        v_hat = triad.vhat_obs_time

        triad_probes = triad.probes
        key = strjoin(sort_uniq(triad_probes),'_')
        if key eq the_key then begin
            r_sms = triad.r_sms
            txs = r_sms[[0,1,2,0],0]
            tys = r_sms[[0,1,2,0],1]
            plots, txs, tys, color=triad_color

            x0 = center_r_sm[0]
            y0 = center_r_sm[1]
            scale = dis_scale/vel_scale*v_mag
            x1 = x0+v_hat[0]*scale
            y1 = y0+v_hat[1]*scale
            arrow, x0,y0,x1,y1,/data, solid=arrow_solid, $
                hsize=arrow_hsize, color=triad_color
        endif
    endforeach

;---Draw flow vel.
    scale = dis_scale/vel_scale
    thm_probes = 'th'+['a','d','e']
    foreach df, df_list, probe_id do begin
        probe = df.probe
        rsm = df.obs_r_sm
        obs_time = df.obs_time
        color = probe_colors[probe_id]
        index = where(thm_probes eq probe, count)
        if count eq 0 then continue

        u_sm = get_var_data(probe+'_u_sm', times=times, in=obs_time+[-1,1]*60*10)

        ux_mean = mean(u_sm[*,0],/nan)
        uy_mean = mean(u_sm[*,1],/nan)
        x0 = rsm[0]
        y0 = rsm[1]
        x1 = x0+ux_mean*scale
        y1 = y0+uy_mean*scale
        arrow, x0,y0, x1,y1, data=1, solid=arrow_solid, $
                hsize=arrow_hsize, color=color
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
        if probe eq 'g15' then begin
            tx = tmp[0]-xchsz*0.5
            ty = tmp[1]
            xyouts, tx,ty,/normal, alignment=1, strupcase(short_name), color=probe_colors[ii]
        endif else begin
            tx = tmp[0]+xchsz*0.5
            ty = tmp[1]
            xyouts, tx,ty,/normal, strupcase(short_name), color=probe_colors[ii]
        endelse
    endforeach
    ; Add LANL.
    _2014_0828_10_load_data
    lanl_color = sgcolor('grey')
    r_gsm = get_var_data('1991-080_r_gsm', at=pos_times)
    r_sm = cotran(r_gsm, pos_times, 'gsm2sm')
    rsm = reform(r_sm[0,*])
    plots, rsm[0],rsm[1],/data, psym=square_symbol, symsize=square_size, color=lanl_color
    short_name = 'lanl'
    tmp = convert_coord(rsm[0:1], /data, /to_normal)
    tx = tmp[0];-xchsz*0.5
    ty = tmp[1]-ychsz*1
    xyouts, tx,ty,/normal, alignment=0.5, strupcase(short_name), color=lanl_color


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


    ; The XZ plane.
    tpos = right_poss[*,0]
    times = map_time+dblarr(nline)
    ndim = 3
    r_sm = dblarr(nline,ndim)
    foreach line, line_info, ii do begin
        lat = line.lat
        r_sm[ii,*] = [-cos(lat),0,sin(lat)]
    endforeach

    h0 = 100.   ; km.
    re = constant('re')
    r0 = 1+h0/re
    r_sm *= r0
    r_gsm = cotran(r_sm, times, 'sm2gsm')
    flines = list()


    xrange = pos_xrange
    yrange = pos_zrange
    xstep = 5
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    ystep = 5
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = ystep

    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ytitle = 'SM Z (Re)'

    ; Set up coord.
    plot, xrange, yrange, nodata=1, noerase=1, $
        xstyle=1, xrange=xrange, xticklen=1, xgridstyle=1, xtickformat='(A1)', $
        xticks=xticks, xtickv=xtickv, xminor=1, $
        ystyle=1, yrange=yrange, yticklen=1, ygridstyle=1, ytickformat='(A1)', $
        yticks=yticks, ytickv=ytickv, yminor=1, $
        position=tpos, /iso, color=sgcolor('silver')


    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)>min(yrange)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    plots, circle_x, circle_y


    ; Prepare model info.
    sgeopack_par, map_time+[-1,1]*600, model
    par = get_var_data(model+'_par', at=map_time)
    routine = (model eq 't04s')? 'ts04': model
    routine = 'geopack_'+routine
    t89 = (model eq 't89')? 1: 0
    t96 = (model eq 't96')? 1: 0
    t01 = (model eq 't01')? 1: 0
    t04s = (model eq 't04s')? 1: 0
    storm = (model eq 't04s')? 1: 0


    min_bmag = 20.
    min_bmag_color = sgcolor('dark_gray')
    min_bmag_color = sgcolor('silver')
    foreach line, line_info, ii do begin
        dir = (r_gsm[ii,2]>0)? 1: -1
        geopack_trace, r_gsm[ii,0], r_gsm[ii,1], r_gsm[ii,2], dir, par, xf, yf, zf, $
            fline=fline, /igrf, r0=r0, t89=t89, t96=t96, t01=t01, ts04=t04s, storm=storm
        ndata = n_elements(fline[*,0])
        fline = cotran(fline, map_time+dblarr(ndata), 'gsm2sm')
        oplot, fline[*,0], fline[*,2], color=line.color

        ndata = n_elements(fline[*,0])
        b_gsm = fltarr(ndata,ndim)
        for jj=0,ndata-1 do begin
            geopack_igrf_gsm, fline[jj,0],fline[jj,1],fline[jj,2], bxp,byp,bzp
            call_procedure, routine, par, fline[jj,0],fline[jj,1],fline[jj,2], tbx,tby,tbz
            b_gsm[jj,*] = [bxp,byp,bzp]+[tbx,tby,tbz]
        endfor
        bmag = snorm(b_gsm)
        index = where(bmag le min_bmag, count)
        if count ne 0 then oplot, fline[index,0], fline[index,2], color=min_bmag_color
    endforeach


    plot, xrange, yrange, /nodata, /noerase, $
        xstyle=1, xrange=xrange, xticklen=xticklen, xtitle=' ', $
        xticks=xticks, xtickv=xtickv, xminor=xminor, xtickformat='(A1)', $
        ystyle=1, yrange=yrange, yticklen=yticklen, ytitle=ytitle, $
        yticks=yticks, ytickv=ytickv, yminor=yminor, $
        position=tpos, /iso


    yspace = 0.9
    msg = 'J Downward!C '+strjoin(string(down_current.lat_range,'(I0)'),'-')+' deg'
    tx = -8
    ty = 4
    xyouts, tx,ty,alignment=0.5, /data, msg, color=sgcolor('blue'), charsize=label_size
    msg = 'J Upward!C'+strjoin(string(up_current.lat_range,'(I0)'),'-')+' deg'
    tx = -5.7
    ty = 2.0
    xyouts, tx,ty,alignment=0.5, /data, msg, color=sgcolor('red'), charsize=label_size
    msg = '|B| < 20 nT!CApporx. PS'
    tx = -12.5
    ty = 1.5
    xyouts, tx,ty,alignment=0.5, /data, msg, color=min_bmag_color, charsize=label_size


    tx = tpos[0]
    ty = tpos[3]+ychsz*0.5
    xyouts, tx,ty, /normal, 'b. XZ plane and '+strupcase(model)+' model B field'


;---Left panel.
    xrange = short_time_range
    xminor = 10 ; min.
    xstep = 30*60
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    secofday = 86400d
    xtickn = time_string(xtickv,tformat='hh:mm')
    xtickn[0] = time_string(xrange[0],tformat='YYYY-MM-DD')+'        '

    poss = sgcalcpos(ndf, position=left_pos)
    foreach probe, sorted_probes, probe_id do begin
        prefix = probe+'_'
        theta = get_var_data(prefix+'theta', times=times)
        index = where_pro(times, '[]', short_time_range)
        yrange = minmax(theta[index])
        yrange = [-1,1]*max(abs(make_bins(yrange,2)))
        yticks = 2
        if max(yrange) le 5 then begin
            yminor = max(yrange)
            ytickv = [-1,0,1]*yminor
        endif else begin
            yminor = 5
            ytickv = [-1,0,1]*floor(max(yrange)/5)*5
        endelse
        xtickformat = (probe_id eq ndf-1)? '': '(A1)'


        tpos = poss[*,probe_id]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        ytickformat = ''
        ytitle = '(deg)'
        the_time_range = short_time_range

        index = where_pro(times, '[]', the_time_range)
        xxs = times[index]
        yys = theta[index]

        the_xtickn = xtickn

        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, /nodata, /noerase

        oplot, xxs, yys;, color=line_colors[jj]
        plots, xrange, [0,0], linestyle=1;, color=line_colors[jj]
        ref_time = obs_times[probe_id]
        plots, ref_time, yrange, linestyle=1;, color=line_colors[jj]

        ; Add ramp time.
        tmp = convert_coord(ref_time,yrange[0], /data, /to_normal)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.25
        xyouts, tx,ty,/normal, time_string(ref_time,tformat='hh:mm:ss'), charsize=label_size


        ; Draw box.
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=the_xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytickformat=ytickformat, ytitle=ytitle, $
            position=tpos, /nodata, /noerase

        ; Add probe info.
        tx = tpos[0]-xchsz*7
        ty = (tpos[1]+tpos[3])*0.5+ychsz*half_ychsz
        probe_info = resolve_probe(probe)
        xyouts, tx,ty,/normal, alignment=1, strupcase(probe_info.short_name), color=probe_colors[probe_id]
        df = df_list[probe_id]
        ref_mlt = df.obs_mlt
        ref_rxy = df.obs_rxy
        ty = (tpos[1]+tpos[3])*0.5-ychsz*half_ychsz
        xyouts, tx,ty,/normal, alignment=1, sgnum2str(ref_mlt,ndec=1)+' MLT!C'+$
            sgnum2str(ref_rxy,ndec=1)+' Re', charsize=label_size

        ; Add title.
        if probe_id eq 0 then begin
            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            title = 'a. Detrended tilt angle '+tex2str('theta')+' of magnetic field'
            xyouts, tx,ty,/normal, title, alignment=0
        endif

        ; Add label.
        label = 'a-'
        label += string(probe_id+1,format='(I0)')+'.'
        ty = tpos[3]-ychsz*1
        tx = tpos[0]+xchsz*0.5
        alignment = 0
        xyouts, tx,ty,/normal, label, alignment=alignment


        ; Add width of DP and dipolarized region.
        dy = ychsz*0.25
        width_color = sgcolor('purple')
        if probe eq 'g15' then begin
            ; Dipolarization.
            dp_range = df.time_range
            txs = dblarr(2)
            foreach tx, dp_range, ii do begin
                tmp = convert_coord(tx,0, data=1, to_normal=1)
                txs[ii] = tmp[0]
            endforeach
            ty = tpos[3]-ychsz*0.4
            plots, txs, ty+[0,0], normal=1, color=width_color
            foreach ii,[0,1] do plots, txs[ii], ty+[-1,1]*0.5*dy, normal=1, color=width_color
            msg = 'dipolarization '+strtrim(string(total(dp_range*[-1,1]/60),format='(F5.1)'),2)+' min'
            xyouts, txs[1]+xchsz*0.5,ty-ychsz*0.3,normal=1, msg, charsize=label_size, color=width_color

            ; Dipolarized region.
            dr_range = df.obs_time+[0.,25*60]
            txs = dblarr(2)
            foreach tx, dr_range, ii do begin
                tmp = convert_coord(tx,0, data=1, to_normal=1)
                txs[ii] = tmp[0]
            endforeach
            ty = 0
            tmp = convert_coord(0,ty, data=1, to_normal=1)
            ty = tmp[1]-ychsz*0.5
            plots, txs, ty+[0,0], normal=1, color=width_color
            foreach ii,[0,1] do plots, txs[ii], ty+[-1,1]*0.5*dy, normal=1, color=width_color
            msg = 'dipolarized region '+strtrim(string(total(dr_range*[-1,1]/60),format='(I0)'),2)+' min'
            xyouts, txs[1]+xchsz*0.5,ty-ychsz*0.3,normal=1, msg, charsize=label_size, color=width_color
        endif
    endforeach


    if keyword_set(test) then stop
    sgclose

end
