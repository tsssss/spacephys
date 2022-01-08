;+
; Return a event overlaps with given time range.
;-

function azim_df_find_dfgroup, project=project, reset=reset

    retval = list()
    if n_elements(project) eq 0 then project = azim_df_load_project()
    file_suffix = project.name+'_find_dfgroup.txt'
    out_file = join_path([project.data_dir,file_suffix])
    if keyword_set(reset) then file_delete, out_file, /allow_nonexistent
    if file_test(out_file) eq 1 then begin
        if file_lines(out_file) eq 0 then file_delete, out_file
    endif

    if file_test(out_file) eq 0 then begin
        ftouch, out_file
        
        events = list()
        candidates = azim_df_find_subgroup(project=project)
        foreach candidate, candidates do begin
            dir = azim_df_analyze_direction(candidate)
            if dir eq '' then continue
            candidate.region += '%'+dir
            events.add, candidate
        endforeach

    ;---Output.
        foreach event, events, ii do events[ii].id = ii+1
        foreach event, events do azim_df_subgroup_write_file, event, filename=out_file
    endif

    events = azim_df_subgroup_read_file(out_file)
    ndf = 0
    foreach event, events do ndf += event.df_list.length
    lprmsg, '# of events: '+string(events.length,format='(I0)')
    lprmsg, '# of DFs: '+string(ndf,format='(I0)')

    directions = strarr(events.length)
    foreach event, events, ii do directions[ii] = (strsplit(event.region,'%',/extract))[1]
    all_directions = sort_uniq(directions)
    foreach direction, all_directions do begin
        index = where(directions eq direction, count)
        lprmsg, '# of event along '+direction+': '+string(count,format='(I0)')
    endforeach

    return, events

end

reset = 0
events = azim_df_find_dfgroup(project=project, reset=reset)
end