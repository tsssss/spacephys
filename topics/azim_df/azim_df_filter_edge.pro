;+
; Check if a given candidate is a DF group.
; Adopted from azim_df_filter_df_group.
;-

function azim_df_filter_edge, candidate, project=project, $
    settings=settings, $
    log_file=log_file

    tab = constant('4space')
    if n_elements(log_file) eq 0 then log_file = -1
    if n_elements(settings) eq 0 then begin
        settings = dictionary($
            'dtime_range', [0.,30]*60, $
            'dtime_min_mean', 8.*60, $
            'min_mean_xcor', 0.5 )
    endif

    dtime_range = settings.dtime_range
    dtime_min_mean = settings.dtime_min_mean
    min_mean_xcor = settings.min_mean_xcor
    if n_elements(candidate) eq 0 then begin
        lprmsg, 'Settings for filtering edge ...', log_file
        lprmsg, 'dtime range (sec): ['+strjoin(string(dtime_range,format='(I0)'),',')+']', log_file
        lprmsg, 'Minimum dtime_mean (sec): '+string(dtime_min_mean,format='(I0)'), log_file
        lprmsg, 'Minimum xcors_mean (#): '+string(min_mean_xcor,format='(F3.1)'), log_file
        lprmsg, '', log_file
        return, settings
    endif

    retval = dictionary()
    if n_elements(project) eq 0 then project = azim_df_load_project()
    azim_df_load_basic_data, project=project, scale_width=scale_width
    lprmsg, '', log_file
    msg = 'Processing candidate: '+strjoin(time_string(candidate.time_range),' to ')
    lprmsg, msg, log_file


;---Check dtimes.
    dtime_max = dtime_range[1]
    dtime_min = 0.
    df_list = candidate.df_list
    obs_times = dblarr(df_list.length)
    foreach df, df_list, df_id do obs_times[df_id] = df.obs_time
    dtimes = obs_times[1:df_id-1]-obs_times[0:df_id-2]
    msg = tab+'dtimes (min): '+strjoin(string(dtimes,format='(F6.1)'),', ')
    lprmsg, msg, log_file
    ; maximum dtime.
    index = where(dtimes ge dtime_max, count)
    if count ne 0 then begin
        msg = tab+'dtime too large, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif
    ; mean dtime.
    index = where(dtimes gt dtime_min, count)
    msg = tab+'# >'+string(dtime_min,format='(I0)')+' sec: '+string(count,format='(I0)')
    mean_dtime = mean(dtimes[index])
    msg = tab+'dtimes_mean (min): '+string(mean_dtime,format='(F6.1)')
    if mean_dtime ge dtime_min_mean then begin
        msg = tab+'dtime mean too large, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif


;---Gatther info.
    probe_list = list()
    df_list = candidate.df_list
    foreach df, df_list do begin
        probe = df.probe
        if probe_list.where(probe) eq !null then probe_list.add, probe
    endforeach
    nprobe = probe_list.length
    probes = probe_list.toarray()
    ndim = 3
    df_r_sms = fltarr(nprobe,ndim)
    foreach df, df_list, ii do df_r_sms[ii,*] = df.obs_r_sm
    df_obs_times = dblarr(nprobe)
    foreach df, df_list, ii do df_obs_times[ii] = df.obs_time
    df_widths = fltarr(nprobe)
    foreach df, df_list, ii do df_widths[ii] = df.width


;---Check edges.
    nedge_vertex = 2
    edge_combos = list()
    for ii=0,nprobe-2 do edge_combos.add, probes[ii:ii+1]
    edge_list = list()
    common_data_rate = project.time_step
    foreach edge_probes, edge_combos do begin
        edge_probes = edge_probes[sort(edge_probes)]

        ; Get info about the vertex.
        probe_index = intarr(nedge_vertex)
        foreach probe, edge_probes, ii do probe_index[ii] = where(probes eq probe)
        vertex_obs_times = df_obs_times[probe_index]
        vertex_widths = df_widths[probe_index,*]

        ; Sort by obs_time.
        index = sort(vertex_obs_times)
        edge_probes = edge_probes[index]
        vertex_obs_times = vertex_obs_times[index]
        vertex_widths = vertex_widths[index]
        ; Must be positive or 0.
        time_lag_df = total(vertex_obs_times*[-1,1])

    ;---xcor.
        ; This setting works.
        the_width = min(vertex_widths)
        the_time_range = the_width*[-1,1]*2.5
        nlag = floor(the_width*2/common_data_rate)
        section_times = make_bins(the_time_range, common_data_rate)
        nsection_time = n_elements(section_times)

        ; Read data.
        nrec = n_elements(section_times)
        data = fltarr(nrec,nedge_vertex)
        for probe_id=0,nedge_vertex-1 do begin
            the_times = vertex_obs_times[probe_id]+section_times
            the_var = edge_probes[probe_id]+'_theta'
            data[*,probe_id]= get_var_data(the_var, at=the_times)
        endfor

        ; Do xcor.
        lags = findgen(nlag)-round(nlag/2)
        xcors = c_correlate(data[*,0], data[*,1], lags)
        xcor_max = max(xcors, index)
        xcor_time_lag = lags[index]*common_data_rate
        ; xcor error.
        xx = total(data,2)/nedge_vertex
        del_x = xx-mean(xx)
        dx_dt = deriv(xx)/common_data_rate
        xcor_err = sqrt(1./(nrec-1)*(1-xcor_max)/xcor_max*2*mean(del_x^2)/mean(dx_dt^2))

        ; Save results.
        time_lag_df = total(vertex_obs_times*[-1,1])
        time_lag_xcor = time_lag_df+xcor_time_lag
        time_lag_diff = abs(time_lag_df-time_lag_xcor)
        time_lag_ratio = time_lag_diff/max(abs([time_lag_df,time_lag_xcor]))

        edge_info = dictionary($
            'probes', edge_probes, $
            'time_lag_df', time_lag_df, $
            'time_lag_xcor', time_lag_xcor, $
            'xcor_max', xcor_max, $
            'xcor_err', xcor_err, $
            'time_lag_diff', time_lag_diff, $
            'time_lag_ratio', time_lag_ratio)
        edge_list.add, edge_info

        ; Check xcor.
        msg = tab+'edge of '+strjoin(edge_probes,'_')+' xcor: '+string(xcor_max,format='(F5.2)')
        lprmsg, msg, log_file
    endforeach

    xcor_maxs = list()
    foreach info, edge_list do xcor_maxs.add, info.xcor_max
    xcor_maxs = xcor_maxs.toarray()
    mean_xcor = mean(xcor_maxs)
    mean_xcor = round(mean_xcor*10)/10.
    msg = tab+'mean xcor: '+string(mean_xcor,format='(F4.1)')
    lprmsg, msg, log_file
    if mean_xcor le min_mean_xcor then begin
        msg = 'Bad coherency, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif


    candidate.edge_list = edge_list
    msg = 'Pass ...'
    lprmsg, msg, log_file
    return, candidate
end
