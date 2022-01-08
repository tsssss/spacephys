;+
; Check edge coherency.
;-


;---Add test times.
    test_event_times = time_double([$
        '2007-11-20/17:18:10', $
        '2008-01-09/11:27:45', $
        '2008-01-19/12:02:55', $
        '2008-02-29/08:26:50', $
        '2014-08-28/10:10:40', $
        '2014-12-26/01:05:25', $
        '2016-10-13/12:22:35', $
        '2016-12-11/09:46:35', $
        '2017-03-28/03:00:40'])
    event_list = list()
    foreach test_event_time, test_event_times do event_list.add, test_event_time


;---Load all candidates.
    search_step = 'uniq_df'
    candidates = list()
    candidates.add, azim_df_search_post_midn_events(search_step=search_step, reset=reset), /extract
    candidates.add, azim_df_search_pre_midn_events(search_step=search_step, reset=reset), /extract


;---Sort by time.
    times = dblarr(candidates.length)
    foreach candidate, candidates, ii do times[ii] = candidate.time_range[0]
    index = sort(times)
    candidates = candidates[index]
    

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


;---Loop through.
;    foreach candidate, test_candidates do begin
     foreach candidate, candidates do begin
        lprmsg, 'Candidate '+strjoin(time_string(candidate.time_range),' to ')+' ...'
        event = azim_df_filter_edge(candidate, project=project)
        if n_elements(event) eq 0 then stop
    endforeach
    stop

end
