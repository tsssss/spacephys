;+
; Do all timing for a given event.
;-

function azim_df_do_timing_per_event, event, project=project

    retval = dictionary()
    event = azim_df_filter_edge(event, project=project)
    if n_elements(event) eq 0 then return, retval
    event = azim_df_filter_triad(event, project=project)
    return, event

end


;---Add test times.
    test_event_times = time_double([$
        '2017-04-09/12:04:15', $
        '2014-08-02/16:54:55'])

    event_list = list()
    foreach test_event_time, test_event_times do event_list.add, test_event_time


;---Load all candidates.
    search_step = 'uniq_subgroup'
    candidates = azim_df_search_all_events(search_step=search_step, project=project)


;---Select the test events.
    index = list()
    foreach test_event_time, test_event_times do begin
        foreach candidate, candidates, id do begin
            if test_event_time eq candidate.time_range[0] then begin
                index.add, id
                lprmsg, 'Found the test event ...'
                break
            endif
        endforeach
    endforeach
    index = index.toarray()
    test_candidates = candidates[index]

    foreach candidate, test_candidates do begin
        tline = 'rbspb     2017-04-09/12:40:40      780       3.37           8.0    -5.27     5.26    2017-04-09/12:34:10 2017-04-09/12:47:10      -2.28    1.09     -1.00,  5.16, -1.73'
        the_df = azim_df_vertex_read(tline)
        df_list = candidate.df_list
        foreach df, df_list, ii do if df.probe eq the_df.probe then df_list[ii] = the_df
        obs_times = dblarr(df_list.length)
        foreach df, df_list, ii do obs_times[ii] = df.obs_time
        index = sort(obs_times)
        candidate.df_list = df_list[index]
        
        event = azim_df_do_timing_per_event(candidate, project=project)
        stop
    endforeach

end