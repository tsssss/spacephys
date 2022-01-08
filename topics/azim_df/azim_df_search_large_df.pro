;+
; Search large DFs.
;-

function azim_df_filter_large_df_per_event, df_list, project=project, $
    settings=settings, log_file=log_file, $
    test=test


    tab = constant('4space')
    if n_elements(log_file) eq 0 then log_file = -1
    if n_elements(settings) eq 0 then begin
        settings = dictionary($
            'probe_min_count', 4, $
            'clear_window_size', 10.*60 )
    endif

    probe_min_count = settings.probe_min_count
    clear_window_size = settings.clear_window_size
    if n_elements(df_list) eq 0 then begin
        lprmsg, 'Settings for filtering large dipolarizations ...', log_file
        lprmsg, tab+'Minimum # of probes (#): '+string(probe_min_count,format='(I0)'), log_file
        lprmsg, tab+'Clear window size (sec): '+string(clear_window_size,format='(I0)'), log_file
        lprmsg, '', log_file
        return, settings
    endif

    retval = list()
    if n_elements(project) eq 0 then project = azim_df_load_project()
    ndf = df_list.length
    obs_times = dblarr(ndf)
    foreach df, df_list, ii do obs_times[ii] = df.obs_time
    time_range = minmax(obs_times)
    lprmsg, '', log_file
    msg = 'Processing event: '+tab+strjoin(time_string(time_range),' to ')
    lprmsg, msg, log_file
    
    ; Sort by probe.
    df_dict = dictionary()
    foreach df, df_list do begin
        probe = df.probe
        if ~df_dict.haskey(probe) then df_dict[probe] = list()
        df_dict[probe].add, df 
    endforeach
    probes = df_dict.keys()

    ; Loop DFs of each probe.
    large_df_list = list()
    foreach probe, probes do begin
        probe_df_list = df_dict[probe]
        nprobe_df = probe_df_list.length
        foreach df, probe_df_list, df_id do begin
        ;---Check df shape.
            large_df = azim_df_filter_vertex(df, project=project, log_file=log_file)
            if n_elements(large_df) eq 0 then continue

        ;---Check clear_window_size.
            if df_id lt nprobe_df-1 then begin
                dtime = probe_df_list[df_id+1].obs_time-df.obs_time
                lprmsg, tab+'Next DF time (sec): '+string(dtime,format='(I0)'), log_file
                if dtime lt clear_window_size then begin
                    lprmsg, 'Within clear window, skip ...', log_file
                    continue
                endif
            endif

        ;---Pass.
            large_df_list.add, large_df
        endforeach
    endforeach
    

;---Test plot.
    plot_file = join_path([project.plot_dir,'diagnostic_plot','azim_df_search_large_df',$
        'azim_df_ut_mlt_plot_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'])
    if file_test(plot_file) eq 0 then azim_df_ut_mlt_plot, df_list, large_df_list, filename=plot_file


;---Filter.
    ndf = large_df_list.length
    lprmsg, tab+'Found '+string(ndf,format='(I0)')+' DFs ...', log_file
    if ndf lt probe_min_count then begin
        lprmsg, 'Not enough DF, skip ...', log_file
        return, retval
    endif

    ; Update probes.
    probe_list = list()
    foreach df, large_df_list do begin
        probe = df.probe
        if probe_list.where(probe) eq !null then probe_list.add, probe
    endforeach
    nprobe = probe_list.length
    lprmsg, tab+'Found '+string(nprobe,format='(I0)')+'probes ...', log_file
    if nprobe lt probe_min_count then begin
        lprmsg, 'Not enough probe, skip ...', log_file
        return, retval
    endif

;---Sort by obs_time
    obs_times = dblarr(ndf)
    foreach df, large_df_list, ii do obs_times[ii] = df.obs_time
    index = sort(obs_times)
    large_df_list = large_df_list[index]

    lprmsg, 'Pass ...', log_file
    return, large_df_list

end

function azim_df_search_large_df, search_setting, project=project, $
    reset=reset, test_time=test_time

;test_time = time_double('2007-11-16/20:00')
    retval = list()
    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()

;---Settings for selecting large DF.
    the_key = 'search_large_df'
    if ~search_setting.haskey(the_key) then message, 'No settings for large DF ...'
    search_info = search_setting[the_key]
    min_probe_count = 4.


;---Check if file exists, to avoid search again.
    file_suffix = search_info.file_suffix
    out_file = join_path([project.data_dir,file_suffix])
    log_file = join_path([project.data_dir,file_suffix+'.log'])
    if keyword_set(reset) then begin
        lprmsg, 'Resetting large DF search ...'
        file_delete, out_file, /allow_nonexistent
        file_delete, log_file, /allow_nonexistent
    endif
    if file_test(out_file) eq 1 then begin
        if file_lines(out_file) eq 0 then file_delete, out_file
    endif
    if file_test(log_file) eq 1 then begin
        if file_lines(log_file) eq 0 then file_delete, log_file
    endif

    if file_test(out_file) eq 0 or keyword_set(test_time) then begin
        event_id = 0
        if file_test(out_file) eq 0 then ftouch, out_file
        if file_test(log_file) eq 0 then ftouch, log_file
        if keyword_set(test_time) then log_file = -1    ; output to console.


    ;---Search settings.
        tmp = azim_df_filter_large_df_per_event(project=project, log_file=log_file)
        tmp = azim_df_filter_vertex(project=project, log_file=log_file)

    ;---Load candidates.
        candidates = azim_df_search_df(search_setting, project=project)
        df_count = 0
        foreach candidate, candidates do begin
            df_list = azim_df_filter_large_df_per_event(candidate.df_list, $
                project=project, log_file=log_file)
            if n_elements(df_list) eq 0 then continue

            
            probe_list = list()
            foreach df, df_list do begin
                probe = df.probe
                if probe_list.where(probe) eq !null then probe_list.add, probe
            endforeach
            probes = probe_list.toarray()
            probes = probes[sort(probes)]

            ndf = df_list.length
            obs_times = dblarr(ndf)
            foreach df, df_list, ii do obs_times[ii] = df.obs_time


            ; Write to file.
            event_id += 1
            df_count += df_list.length
            df_group = dictionary($
                'id', event_id, $
                'time_range', minmax(obs_times), $
                'search_name', candidate.search_name, $
                'region', candidate.region, $
                'probes', probes, $
                'df_list', df_list, $
                'triad_list', list(), $
                'edge_list', list())

            azim_df_subgroup_write_file, df_group, filename=out_file
        endforeach
        lprmsg, '', log_file
        lprmsg, 'Found large df (#): '+string(df_count,format='(I0)')+' ...', log_file
    endif

    lprmsg, ''
    lprmsg, 'Read large DF events from file ...'
    return, azim_df_subgroup_read_file(out_file)

end

;---Test event times.
    test_event_times = list()
    ;test_event_times.add, time_double(['2014-08-28/09:10:00','2014-08-28/13:10:00'])
    ;test_event_times.add, time_double(['2016-10-13/10:10:00','2016-10-13/13:55:00'])
    test_event_times.add, time_double(['2007-11-20/16:20:00','2007-11-20/20:15:00'])



    if candidates.length eq 0 then begin
        search_step = 'df'
        candidates = list()
        candidates.add, azim_df_search_post_midn_events(project=project, search_step=search_step), /extract
        candidates.add, azim_df_search_pre_midn_events(project=project, search_step=search_step), /extract
    endif

    test_events = list()
    foreach time_range, test_event_times do begin
        foreach candidate, candidates do begin
            if total(candidate.time_range-time_range) eq 0 then begin
                test_events.add, candidate
                break
            endif
        endforeach
    endforeach

    log_file = -1
    foreach candidate, test_events do begin
        event = azim_df_filter_large_df_per_event(candidate.df_list, project=project, log_file=log_file)
        stop
    endforeach

end
