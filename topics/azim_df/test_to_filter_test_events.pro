;---Load all df_group.
    project = azim_df_load_project()
    search_step = 'df_group'
    candidates = list()
    candidates.add, azim_df_search_post_midn_events(search_step=search_step), /extract
    candidates.add, azim_df_search_pre_midn_events(search_step=search_step), /extract
;stop
;    dirname = 'df_group'
;    azim_df_subgroup_gen_diagnostic_plot, candidates, dirname=dirname, project=project

;---Sort by time.
    times = list()
    foreach candidate, candidates do times.add, candidate.time_range[0]
    index = sort(times.toarray())
    candidates = candidates[index]
    foreach candidate, candidates, id do candidate.id = id


;---Get overlap_list
    overlap_list = list()
    current_list = list()
    isolated_list = list()
    event_list = list()
    foreach the_group, candidates, id do begin
        str_id = string(id,format='(I0)')
        lprmsg, 'current candidate: '+str_id
        if current_list.length ne 0 then begin
            pre_group = current_list[-1]
            pre_time = pre_group.time_range
            the_time = the_group.time_range
            if min(the_time) ge max(pre_time) then begin
                lprmsg, 'no overlap with previous candidate, pop current_list'
                if current_list.length gt 1 then begin
                    lprmsg, 'add current_list to overlap_list'
                    overlap_list.add, current_list
                    event_list.add, current_list
                endif else begin
                    lprmsg, 'add current_list to isolated_list'
                    isolated_list.add, current_list
                    event_list.add, current_list
                endelse
                lprmsg, 'reset current_list'
                current_list = list()
            endif
        endif
        lprmsg, 'add '+str_id+' to current_list'
        current_list.add, the_group
        if the_group.id eq candidates.length-1 then begin
            if current_list.length gt 1 then begin
                lprmsg, 'add current_list to overlap_list'
                overlap_list.add, current_list
                event_list.add, current_list
            endif else begin
                lprmsg, 'add current_list to isolated_list'
                isolated_list.add, current_list
                event_list.add, current_list
            endelse
        endif
    endforeach

;---Select test_list.
    index = list()
    foreach candidate_list, overlap_list do begin
        start_times = dblarr(event_list.length)
        foreach candidate, candidate_list, ii do start_times[ii] = candidate.time_range[0]
        contain_test_event = 0
        foreach test_event_time, test_event_times do begin
            index = where(start_times eq test_event_time, count)
            the_test_event_time = start_times[index]
            if count ne 0 then begin
                contain_test_event = 1
                break
            endif
        endforeach
        if ~contain_test_event then continue

        selected_event = azim_df_resolve_overlap_subgroup(candidate_list, project=project)
        if selected_event.time_range[0] ne the_test_event_time then begin
            message, 'Did not find the wanted event ...'
        endif else begin
            print, 'Found the wanted event ...'
        endelse
    endforeach
end
