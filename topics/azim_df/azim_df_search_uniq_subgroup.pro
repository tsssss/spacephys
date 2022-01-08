;+
; Resolve overlap events.
;-

function azim_df_search_uniq_subgroup, search_setting, project=project, $
    reset=reset, test_time=test_time

    retval = list()
    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()

;---Settings for selecting large DF.
    the_key = 'search_uniq_subgroup'
    if ~search_setting.haskey(the_key) then message, 'No settings for unique DF ...'
    search_info = search_setting[the_key]

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


    ;---Load candidates.
        candidates = azim_df_search_df_group(search_setting, project=project)

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
            lprmsg, 'current candidate: '+str_id, log_file
            if current_list.length ne 0 then begin
                pre_group = current_list[-1]
                pre_time = pre_group.time_range
                the_time = the_group.time_range
                if min(the_time) ge max(pre_time) then begin
                    lprmsg, 'no overlap with previous candidate, pop current_list', log_file
                    if current_list.length gt 1 then begin
                        lprmsg, 'add current_list to overlap_list'
                        overlap_list.add, current_list
                    endif else begin
                        lprmsg, 'add current_list to isolated_list', log_file
                        isolated_list.add, current_list
                        event_list.add, current_list[0]
                    endelse
                    lprmsg, 'reset current_list', log_file
                    current_list = list()
                endif
            endif
            lprmsg, 'add '+str_id+' to current_list', log_file
            current_list.add, the_group
            if the_group.id eq candidates.length-1 then begin
                if current_list.length gt 1 then begin
                    lprmsg, 'add current_list to overlap_list', log_file
                    overlap_list.add, current_list
                endif else begin
                    lprmsg, 'add current_list to isolated_list', log_file
                    isolated_list.add, current_list
                    event_list.add, current_list[0]
                endelse
            endif
        endforeach

        foreach candidate_list, overlap_list do begin
            event_list.add, azim_df_resolve_overlap_subgroup(candidate_list, project=project)
        endforeach

        times = dblarr(event_list.length)
        foreach event, event_list, ii do times[ii] = event.time_range[0]
        index = sort(times)
        event_list = event_list[index]

        foreach event, event_list, ii do begin
            event.id = ii+1
            azim_df_subgroup_write_file, event, filename=out_file
        endforeach
    endif

    if keyword_set(test_time) then stop
    return, azim_df_subgroup_read_file(out_file)
end
