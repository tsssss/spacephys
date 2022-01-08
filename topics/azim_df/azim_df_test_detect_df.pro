;+
; Load theta for a given time range and a list of probes.
; Detect dipolarizations and return times at each probe.
;-

pro azim_df_test_detect_df, time_range, probes=probes, project=project

test = 1
    test_time = time_double('2014-08-28/10:00')
;    test_time = time_double('2016-10-13/12:30')
    ;test_time = time_double('2008-01-07/08:00')
;    test_time = time_double('2008-02-29/08:00')
    ;test_time = time_double('2017-03-01/05:00')
    candidates = azim_df_search_candidate(project=project)
    if keyword_set(test_time) then begin
        candidate_ids = list()
        foreach candidate, candidates do begin
            if product(candidate.time_range-test_time) gt 0 then continue else candidate_ids.add, candidate.id
        endforeach
        candidate_ids = candidate_ids.toarray()-1   ; id starts from 1.
        candidates = candidates[candidate_ids]
    endif
    candidate = candidates[0]
    time_range = candidate.time_range
    probes = candidate.all_probes
    probe_infos = project.probe_infos
    search_name = candidate.search_type
    foreach tmp, project.search_settings do if tmp.name eq search_name then search_setting = tmp
    full_time_range = search_setting.time_range



; Test for Runov's events.
    search_name = 'beyond_15Re'
    foreach tmp, project.search_settings do if tmp.name eq search_name then search_setting = tmp
    full_time_range = search_setting.time_range

    ; Fig 1-1 in Runov+2010.
    probes = 'th'+letters('e')
    time_range = time_double(['2009-02-27/07:30','2009-02-27/09:00'])
    region_name = 'post_midn'

    ; Fig 1-2 in Runov+2010.
    probes = 'th'+letters('e')
    time_range = time_double(['2009-03-05/02:30','2009-03-05/04:00'])
    region_name = 'pre_midn'

    ; Fig 1-3 in Runov+2010.
    probes = 'th'+letters('e')
    time_range = time_double(['2009-03-09/08:30','2009-03-09/10:00'])
    region_name = 'post_midn'

    ; Fig 2-1 in Runov+2010.
    probes = 'th'+letters('e')
    time_range = time_double(['2009-03-15/08:30','2009-03-15/09:30'])
    region_name = 'post_midn'

    ; Fig 2-2 in Runov+2010.
    probes = 'th'+letters('e')
    time_range = time_double(['2009-03-19/08:00','2009-03-19/09:00'])
    region_name = 'pre_midn'

    ; Fig 2-3 in Runov+2010.
    probes = 'th'+letters('e')
    time_range = time_double(['2009-03-31/08:00','2009-03-31/09:00'])
    region_name = 'pre_midn'

    ; Fig 1 in Runov+2009.
    probes = 'th'+letters('e')
    time_range = time_double(['2009-02-27/07:30','2009-02-27/09:00'])
    region_name = 'post_midn'

    candidate = dictionary($
        'id', 0, $
        'duration', total(time_range*[-1,1]), $
        'search_type', 'beyond_15Re', $
        'nsection', 1, $
        'probe_list', list('th'+letters('e')), $
        'time_range', time_range, $
        'region', region_name, $
        'all_probes', 'th'+letters('e'), $
        'time_range_list', list(time_range))




    plot_file = join_path([project.plot_dir,'diagnostic_plot','search_event','detect_df','detect_df_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmm'),'_')+'.pdf'])

;---Load data.
    foreach probe, probes do begin
        prefix = probe_infos[probe].prefix
        var_suffixes = ['r_sm','theta','mlt']
        foreach var_suffix, var_suffixes do begin
            the_var = prefix+var_suffix
            if check_if_update(the_var, full_time_range) then azim_df_read_data, var_suffix, time_range=full_time_range, probe=probe, project=project
            if check_if_update(prefix+'x_sm', full_time_range) then begin
                get_data, prefix+'r_sm', times, r_sm
                store_data, prefix+'x_sm', times, r_sm[*,0]
            endif
        endforeach
    endforeach
    if check_if_update('dst', full_time_range) then omni_read_index, full_time_range


;---Loop through.
    probe_df = dictionary()
    time_step = project.time_step
    common_times = make_bins(time_range, time_step)
    section_duration = 120. ; sec.
    section_boundary_times = smkarthm(time_range[0], time_range[1], section_duration, 'dx')
    section_times = (section_boundary_times[1:-1]+section_boundary_times[0:-2])*0.5
    nsection = n_elements(section_times)
    boxcar_width = section_duration/time_step
    foreach probe, probes do begin
        prefix = probe_infos[probe].prefix
        df_list = list()

        ; Make a copy to play with.
        theta_orig_var = prefix+'theta'
        theta = get_var_data(theta_orig_var, at=common_times)
        get_data, theta_orig_var, limits=lim
        theta_var = prefix+'theta_test'
        store_data, theta_var, common_times, theta, limits=lim

        ; Scale by MLT.
        mlt_var = prefix+'mlt'
        mlt = get_var_data(mlt_var, at=common_times)
        theta_scale = theta*sqrt(1+abs(mlt))
        theta_scale = theta
        theta_scale_var = prefix+'theta_test_scale'
        store_data, theta_scale_var, common_times, theta_scale, limits=lim

        section_stddev = fltarr(nsection)
        section_median = fltarr(nsection)
        section_mean = fltarr(nsection)
        for ii=0, nsection-1 do begin
            the_data = theta_scale[lazy_where(common_times, '[]', section_boundary_times[ii:ii+1])]
            if n_elements(the_data) eq 1 then begin
                section_stddev[ii] = the_data
                section_median[ii] = the_data
                section_mean[ii] = the_data
            endif else begin
                section_stddev[ii] = stddev(the_data,/nan)
                section_median[ii] = median(the_data)
                section_mean[ii] = mean(the_data,/nan)
            endelse
        endfor


        theta_stddev = interpol(section_stddev, section_times, common_times)
        theta_median = interpol(section_median, section_times, common_times, /quadratic)
        theta_mean = interpol(section_mean, section_times, common_times, /quadratic)
;
;        section_deriv = deriv(section_median)/section_duration
;        section_deriv = smooth(section_deriv,5,/edge_truncate)
;        theta_deriv = interpol(section_deriv, section_times, common_times)

;        theta_stddev = boxcar(theta_scale, boxcar_width, type='stddev')
;        theta_stddev = boxcar(theta_stddev, boxcar_width, type='mean')
        mean_stddev = mean(theta_stddev)
;        theta_median = boxcar(theta_scale, boxcar_width, type='median')    ; moving boxcar median has wired shape of waveform.
;        theta_mean = boxcar(theta_scale, boxcar_width, type='mean')
        theta_wanted = theta_mean
;        theta_wanted = theta_median
        theta_deriv = deriv(theta_wanted)/section_duration
        theta_deriv = boxcar(theta_deriv, boxcar_width, type='median')

        stddev_var = prefix+'theta_stddev'
        store_data, stddev_var, common_times, theta_stddev, limits=lim
        median_var = prefix+'theta_median'
        store_data, median_var, common_times, theta_median, limits=lim
        mean_var = prefix+'theta_mean'
        store_data, mean_var, common_times, theta_mean, limits=lim
        deriv_var = prefix+'theta_deriv'
        store_data, deriv_var, common_times, theta_deriv, limits={constant:0, ytitle:'d'+tex2str('theta')+'/dt (deg/sec)'}


        test_combo_var = prefix+'theta_test_combo'
        store_data, test_combo_var, common_times, [[theta_scale],[theta_median],[theta_mean],[theta_stddev],[-theta_stddev]], limits=lim
        add_setting, test_combo_var, {$
            labels: ['orig','median','mean','stddev',''], $
            constant: 0, $
            colors: sgcolor(['silver','blue','purple','tan','tan'])}
        if n_elements(mean_stddev) ne 0 then options, test_combo_var, 'constant', [0,[-1,1]*mean_stddev]

        ; Find node of positive slope, when stddev>median.
        prod = [0,theta_wanted[1:-1]*theta_wanted[0:-2]]
        index = where(prod lt 0, nnode)
        node_times = common_times[index]
        ; Select nodes with positive slope.
        slope_flags = bytarr(nnode)
        for ii=0, nnode-1 do slope_flags[ii] = theta_wanted[index[ii]-1] lt 0
        index = where(slope_flags eq 1, nnode)
        node_times = node_times[index]

        ; Remove nodes that just "touch" 0.
        node_boundary_times = [time_range[0],node_times,time_range[1]]
        stddev_flags = bytarr(nnode)
;        mean_stddev = mean(theta_stddev)
        nsigma = 0.5
        for ii=0, nnode-1 do begin
            ; Expect a min-value.
            index = lazy_where(common_times, '[]', node_boundary_times[ii+0:ii+1], count=count)
            if count eq 0 then continue
            the_times = common_times[index]
            the_median = theta_median[index]
            the_stddev = theta_stddev[index]
            threshold = the_stddev*nsigma
;            threshold = mean_stddev*nsigma
            index = where(the_median le -threshold, count)
            if count eq 0 then continue
            min_value = min(the_median[index], min_value_index)
            ;if min_value gt -mean_stddev then continue

            ; Expect a max-value.
            index = lazy_where(common_times, '[]', node_boundary_times[ii+1:ii+2], count=count)
            if count eq 0 then continue
            the_times = common_times[index]
            the_median = theta_median[index]
            the_stddev = theta_stddev[index]
            threshold = the_stddev*nsigma
;            threshold = mean_stddev*nsigma
            index = where(the_median ge  threshold, count)
            if count eq 0 then continue
            max_value = max(the_median[index], max_value_index)
            ;if max_value lt  mean_stddev then continue

            stddev_flags[ii] = 1
        endfor
        index = where(stddev_flags eq 1, nnode)
        if nnode eq 0 then continue
        node_times = node_times[index]

        ; Find min/max around a node.
        node_boundary_times = [time_range[0],node_times,time_range[1]]
        df_list = list()
        deriv_prod = [0,theta_deriv[1:-1]*theta_deriv[0:-2]]
        for ii=0, nnode-1 do begin
            ; Find the min/max value.
            ; Local minimum/maximum do not work, b/c there are "hidden" minimum/maximum in a slope.
            ; Try to use derivative to get the min/max values. Works well!!!
            node_time = node_boundary_times[ii+1]
            node_index = where(common_times eq node_time)
            ; Find max value and its time.
            index = lazy_where(common_times, '[]', node_boundary_times[ii+1:ii+2], count=count)
            if count eq 0 then continue
            the_times = common_times[index]
            the_flags = (deriv_prod[index] le 0) and (theta_median[index] ge  theta_stddev[index]*nsigma)
            index = where(the_flags eq 1, count)
            if count eq 0 then continue
            max_value_index = index[0]
            max_value_time = the_times[max_value_index]
            ; Sometimes the found point is not the minimum.
            index = lazy_where(common_times, '[]', [node_time,max_value_time], count=count)
            if count eq 0 then message, 'Inconsistency, stop here'
            the_times = common_times[index]
            max_value = max(theta_median[index], max_value_index)
            max_value_time = the_times[max_value_index]
            ; Find min value and its time.
            index = lazy_where(common_times, '[]', node_boundary_times[ii+0:ii+1], count=count)
            if count eq 0 then continue
            the_times = common_times[index]
            the_flags = (deriv_prod[index] le 0) and (theta_median[index] le -theta_stddev[index]*nsigma)
            index = where(the_flags eq 1, count)
            if count eq 0 then continue
            min_value_index = index[count-1]
            min_value_time = the_times[min_value_index]
            index = lazy_where(common_times, '[]', [min_value_time,node_time], count=count)
            if count eq 0 then message, 'Inconsistency, stop here'
            the_times = common_times[index]
            min_value = min(theta_median[index], min_value_index)
            min_value_time = the_times[min_value_index]

            ; Find the width and height.
            width = max_value_time-min_value_time
            height = max_value-min_value



;            ; Find the width and height.
;            index = lazy_where(common_times, '[]', [min_value_time,node_time])
;            the_times = common_times[index]
;            the_median = theta_median[index]
;            half_min_value = min_value*0.5
;            index = where(the_median le half_min_value)
;            half_min_value_time = the_times[index[-1]]
;
;            index = lazy_where(common_times, '[]', [node_time,max_value_time])
;            the_times = common_times[index]
;            the_median = theta_median[index]
;            half_max_value = max_value*0.5
;            index = where(the_median ge half_max_value)
;            half_max_value_time = the_times[index[0]]
;            width = half_max_value_time-half_min_value_time
;            height = half_max_value-half_min_value

            df_list.add, dictionary($
                'time', node_time, $
                'min_value', min_value, $
                'max_value', max_value, $
                'min_value_time', min_value_time, $
                'max_value_time', max_value_time, $
;                'half_min_value', half_min_value, $
;                'half_max_value', half_max_value, $
;                'half_min_value_time', half_min_value_time, $
;                'half_max_value_time', half_max_value_time, $
                'width', width, $
                'height', height)
        endfor

        probe_df[probe] = df_list

        ; Test plot.
        tpos = [0.1,0.2,0.9,0.9]
        xrange = time_range
        yrange = sg_autolim(minmax(theta_scale))
        options, test_combo_var, 'yrange', yrange
        options, test_combo_var, 'ytitle', strupcase(probe)+' detr.tilt!C'+tex2str('theta')+' (deg)'
;        options, test_combo_var, 'constant', [[-1,1]*mean_stddev,0]
        test_file = join_path([project.plot_dir,'test_plot','azim_df_test_detect_df_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmm')+'_')+'_'+probe+'.pdf'])
        sgopen, test_file, xsize=10, ysize=5, /inch
        tplot, test_combo_var, trange=time_range, position=tpos
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, /nodata, /noerase
        foreach df, df_list do begin
            xxs = [df.min_value_time,df.max_value_time]
            yys = [df.min_value,df.max_value]
            plots, df.time+[0,0], yrange, color=sgcolor('red')
            plots, mean(xxs), yys, color=sgcolor('green')
            plots, xxs, mean(yys)+[0,0], color=sgcolor('green')
        endforeach
        sgclose
;        stop
    endforeach
    min_height = 0
    max_height = 10000.
    max_value_threshold =  0
    min_value_threshold = -0


;---General settings.
    time_range = candidate.time_range
    all_probes = candidate.all_probes
    region_name = candidate.region
    nsection = candidate.nsection
    the_region = project.search_candidate.search_roi.regions[region_name]
    search_name = candidate.search_type
    search_settings = project.search_settings
    search_setting = (search_settings[0].name eq search_name)? search_settings[0]: search_settings[1]
    time_range_list = candidate.time_range_list
    probe_list = candidate.probe_list
    probe_infos = project.probe_infos
    xticklen_chsz = -0.15   ; in ychsz.
    yticklen_chsz = -0.30   ; in xchsz.
    full_ychsz = constant('full_ychsz')
    half_ychsz = full_ychsz*0.5
    lineskip = constant('lineskip')
    label_size = 0.7
    label_xshift = 10
    label_yshift = full_ychsz
    secofday = constant('secofday')
    bar_thick = keyword_set(test)? 0.5: 4

    rad = constant('rad')
    deg = constant('deg')


;---Size of panel and figure.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ; 3 line panels: Dst/AE, MLT-UT, X-UT.
    line_panel_ysize = 1.5
    ypans = [0.5,1,1]*line_panel_ysize
    nypanel = n_elements(ypans)
    ypads = dblarr(nypanel-1)+0.5
    line_panel_aspect_ratio = 0.6     ; inch per hour.
    duration = total(time_range*[-1,1])/constant('secofhour')
    line_pos_xsize = line_panel_ysize*line_panel_aspect_ratio*duration
    line_pos_ysize = total(ypans)+total(ypads)*abs_ychsz
    ; 1 orbit panel: X-Y plane.
    orbit_pos_xsize = line_pos_ysize
    orbit_pos_ysize = line_pos_ysize
    ; Space b/w line and orbit panels.
    xpads = [15]
    ; Figure size.
    margins = [12,5,10,2]
    nxpanel = 2
    pos_xsize = line_pos_xsize+total(xpads)*abs_xchsz+orbit_pos_xsize
    pos_ysize = line_pos_ysize
    fig_xsize = pos_xsize+total(margins[[0,2]])*abs_xchsz
    fig_ysize = pos_ysize+total(margins[[1,3]])*abs_ychsz
    ; Panel positions.
    pos = [margins[0]*abs_xchsz/fig_xsize,margins[1]*abs_ychsz/fig_ysize, $
        1-margins[2]*abs_xchsz/fig_xsize,1-margins[3]*abs_ychsz/fig_ysize]
    line_pos = [pos[2]-line_pos_xsize/fig_xsize,pos[1],pos[2],pos[3]]
    orbit_pos = [pos[0],pos[1],pos[0]+orbit_pos_xsize/fig_xsize,pos[3]]
    pos_list = list()
    pos_list.add, orbit_pos
    pos_tops = dblarr(nypanel)+line_pos[3]
    for ii=1, nypanel-1 do pos_tops[ii] = pos_tops[ii-1]-(ypans[ii-1]+ypads[ii-1]*abs_ychsz)/fig_ysize
    for ii=0, nypanel-1 do begin
        pos_top = pos_tops[ii]
        tpos = [line_pos[0],pos_top-ypans[ii]/fig_ysize,line_pos[2],pos_top]
        pos_list.add, tpos
    endfor

    if n_elements(plot_file) eq 0 then plot_file = 0
;    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz
;    if keyword_set(test) then foreach tpos, pos_list do plot, [0,1], [0,1], /nodata, /noerase, position=tpos


;---Panel 1: SC orbit in XY plane.
    tpos = pos_list[0]

    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xrange = [10.,-30]
    xstep = 10
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = 5
    xtitle = 'SM X (Re)'
    xticklen = xticklen

    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    yrange = (region_name eq 'pre_midn')? [25.,-15]: [15.,-25]
    ystep = 20
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 5
    ytitle = 'SM Y (Re)'
    yticklen = yticklen

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase


    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    plots, circle_x, circle_y

    ; Add magnetopause.
    magn_test = fltarr(nangle,3)
    magn_test[*,0] = circle_x*30
    magn_test[*,1] = circle_y*30
    tmp = check_if_in_magn(magn_test, magn_pos=magn_pos, dynamic_pressure=10)
    magn_pos = magn_pos[*,[0,2,1]]
    magn_pos = [magn_pos,magn_pos]
    magn_pos[nangle:nangle*2-1,1] *= -1
    oplot, magn_pos[*,0], magn_pos[*,1]

    ; Add ROI.
    mlt_range = the_region.mlt_range
    rxy_range = project.search_candidate.search_roi.rxy_range
;    rxy_range = [4,30]
    foreach tmp, [5,10,20] do oplot, circle_x*tmp, circle_y*tmp, linestyle=1
    foreach tmp, rxy_range do oplot, circle_x*tmp, circle_y*tmp
    foreach tmp, mlt_range do oplot, -rxy_range*cos(tmp*15*rad), -rxy_range*sin(tmp*15*rad)

    ; Draw probe position for the full_time_range.
    foreach probe, all_probes do begin
        prefix = probe_infos[probe].prefix
        color = probe_infos[probe].color
        rsm = get_var_data(prefix+'r_sm', in=time_range)
        xxs = rsm[*,0]
        yys = rsm[*,1]
        oplot, xxs, yys, color=color
    endforeach

    ; Highlight the sections and probes within ROI.
    for ii=0, nsection-1 do begin
        section_time_range = time_range_list[ii]
        section_probes = probe_list[ii]
        foreach probe, section_probes do begin
            prefix = probe_infos[probe].prefix
            color = probe_infos[probe].color
            rsm = get_var_data(prefix+'r_sm', in=time_range)
            xxs = rsm[*,0]
            yys = rsm[*,1]
            index = where(section_probes eq probe, count)
            oplot, xxs, yys, color=color, thick=bar_thick
        endforeach
    endfor

    ; Add notations.
    ty0 = (region_name eq 'post_midn')? tpos[3]: tpos[1]+ychsz*full_ychsz*(3+lineskip)
    tx = tpos[0]+xchsz*1
    step = 4

    ty = ty0-ychsz*full_ychsz*1
    xyouts, tx,ty,/normal, 'Probes of the candidate: ', charsize=label_size
    ttx = tx+15*xchsz*label_size
    foreach probe, all_probes do begin
        probe_info = probe_infos[probe]
        xyouts, ttx,ty,/normal, strupcase(probe_info.short_name), color=probe_info.color, charsize=label_size
        ttx += step*xchsz*label_size
    endforeach

    ty = ty0-ychsz*full_ychsz*2
    xyouts, tx,ty,/normal, 'Probes of the search: ', charsize=label_size
    ttx = tx+15*xchsz*label_size
    foreach probe, search_setting.probes do begin
        probe_info = probe_infos[probe]
        xyouts, ttx,ty,/normal, strupcase(probe_info.short_name), color=probe_info.color, charsize=label_size
        ttx += step*xchsz*label_size
    endforeach

    ty = ty0-ychsz*full_ychsz*3
    duration = total(time_range*[-1,1])/constant('secofhour')
    msg = 'Orbit from '+strjoin(time_string(time_range,tformat='YYYY-MM-DD/hh:mm'), ' to ')+', duration is '+sgnum2str(duration,ndec=1)+' hour(s)'
    xyouts, tx,ty,/normal, msg, charsize=label_size

    ; Add axes.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtitle=xtitle, xticklen=xticklen, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, yticklen=yticklen, $
        position=tpos, /noerase, /nodata

    ; Add labels.
    fig_label = 'a. XY!C    plane'
    tx = tpos[0]-label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, fig_label




;---The common x-axis setting.
    xrange = time_range
    xstep = constant('secofhour')
    xticks = floor(total(xrange*[-1,1])/xstep)
    xtickv = smkarthm(xrange[0], xstep, xticks+1, 'x0')
    xminor = 6
    xtickn = strarr(xticks+1)
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


;---Panel b. Dst/AE.
    tpos = pos_list[1]
    xtickformat = '(A1)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    the_var = 'dst'
    ystep = 50.
    constants = [-50,0]
    yys = get_var_data(the_var, in=time_range, times=xxs)
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
    yys = get_var_data(the_var, in=time_range, times=xxs)
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

    tx = tpos[0]-label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'b. Dst/'
    xyouts, tx,ty,/normal, '           AE', color=ae_color


;---Common z-range, colorbar, symbol.
    theta_var = 'theta_median'
    ztitle = 'Lin.Log['+tex2str('theta')+' detr.tilt (deg)]'
    theta_range = [-1,1]*64
    zrange = theta_range
    ztickv = [-64,-16,-4,0,4,16,64]
    ztickn = string(ztickv,format='(I0)')
    ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
    zrange = alog10(abs(zrange))*[-1,1]
    zticks = n_elements(ztickv)-1
    spec_psym = 8
    spec_symsize = 0.3
    usersym, [1,1,-1,-1,1]*0.5, [-1,1,1,-1,-1]*1, /fill
    spec_ct = 70


;---Common x data.
    time_step = 30.
    xxs = make_bins(time_range, time_step)
    nxx = n_elements(xxs)


;---MLT-UT.
    tpos = pos_list[2]
    pos_var = 'mlt'
    xtickformat = '(A1)'
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    ytitle = 'MLT (hr)'
    yminor = 3
    yrange = (region_name eq 'pre_midn')? mlt_range+[0,1]*yminor: mlt_range+[-1,0]*yminor
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    constants = 0


    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Plot data.
    foreach probe, all_probes do begin
        prefix = probe_infos[probe].prefix
        yys = get_var_data(prefix+pos_var, at=xxs)
        zzs = get_var_data(prefix+theta_var, at=xxs)
        mlt = get_var_data(prefix+'mlt', at=xxs)
        zzs = azim_df_normalize_theta(zzs, mlt, zrange=theta_range, ct=spec_ct, /reverse_ct)

        ; Remove data outside ROI.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        rsm = get_var_data(prefix+'r_sm', at=xxs)
        rxy = snorm(rsm[*,0:1])
        index = lazy_where(rxy, '][', rxy_range, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        index = where(finite(zzs,/nan), count)
        if count ne 0 then yys[index] = !values.f_nan
        index = where(finite(yys), count)
        if count eq 0 then continue

        ; Plot data.
        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=spec_psym, symsize=spec_symsize
    endforeach


    ; Add DF times.
    foreach probe, all_probes do begin
        prefix = probe_infos[probe].prefix

        if ~probe_df.haskey(probe) then continue
        df_list = probe_df[probe]
        df_times = list()
        foreach df, df_list do begin
            if df.max_value lt max_value_threshold and df.min_value gt min_value_threshold then continue
            df_times.add, df.time
        endforeach
        txxs = df_times.toarray()
        if n_elements(txxs) eq 0 then continue
        yys = get_var_data(prefix+pos_var, at=txxs)
        mlt = get_var_data(prefix+'mlt', at=txxs)

        ; Remove data outside ROI.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        rsm = get_var_data(prefix+'r_sm', at=txxs)
        if n_elements(txxs) eq 1 then rxy = snorm(rsm[0:1]) else rxy = snorm(rsm[*,0:1])
        index = lazy_where(rxy, '][', rxy_range, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        index = where(finite(yys), count)
        if count eq 0 then continue

        ; Plot data.
        plots, txxs, yys, psym=1, symsize=spec_symsize, color=sgcolor('black')
    endforeach


    ; Add notations.
    foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Add label.
    tx = tpos[0]-label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'c. UT-MLT'

    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;---X-UT.
    tpos = pos_list[3]
    pos_var = 'x_sm'
    xtickformat = ''
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    ytitle = 'SM X (Re)'
    yrange = (search_name eq 'beyond_15Re')? [-30,10]: [-20,10]
    ystep = (search_name eq 'beyond_15Re')? 20: 10
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 5


    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Plot data.
    foreach probe, all_probes do begin
        prefix = probe_infos[probe].prefix
        rsm = get_var_data(prefix+'r_sm', at=xxs)
        yys = rsm[*,0]
        zzs = get_var_data(prefix+theta_var, at=xxs)
        mlt = get_var_data(prefix+'mlt', at=xxs)
        zzs = azim_df_normalize_theta(zzs, mlt, zrange=theta_range, ct=spec_ct, /reverse_ct)

        ; Remove data outside ROI.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        rsm = get_var_data(prefix+'r_sm', at=xxs)
        index = lazy_where(rxy, '][', rxy_range, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        index = where(finite(zzs,/nan), count)
        if count ne 0 then yys[index] = !values.f_nan
        index = where(finite(yys), count)
        if count eq 0 then continue

        ; Plot data.
        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=spec_psym, symsize=spec_symsize
    endforeach


    ; Add DF times.
    foreach probe, all_probes do begin
        prefix = probe_infos[probe].prefix

        if ~probe_df.haskey(probe) then continue
        df_list = probe_df[probe]
        df_times = list()
        foreach df, df_list do begin
            if df.max_value lt max_value_threshold and df.min_value gt min_value_threshold then continue
            df_times.add, df.time
        endforeach
        txxs = df_times.toarray()
        if n_elements(txxs) eq 0 then continue
        yys = get_var_data(prefix+pos_var, at=txxs)
        mlt = get_var_data(prefix+'mlt', at=txxs)

        ; Remove data outside ROI.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        rsm = get_var_data(prefix+'r_sm', at=txxs)
        if n_elements(txxs) eq 1 then rxy = snorm(rsm[0:1]) else rxy = snorm(rsm[*,0:1])
        index = lazy_where(rxy, '][', rxy_range, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        index = where(finite(yys), count)
        if count eq 0 then continue

        ; Plot data.
        plots, txxs, yys, psym=1, symsize=spec_symsize, color=sgcolor('black')
    endforeach


    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Add label.
    tx = tpos[0]-label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'd. UT-X!USM'

    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;    if keyword_set(test) then stop
    sgclose

end
