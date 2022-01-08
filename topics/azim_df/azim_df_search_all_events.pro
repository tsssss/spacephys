
function azim_df_search_all_events, search_step=search_step, project=project, reset=reset

    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(search_step) eq 0 then search_step = 'uniq_subgroup'


    file_suffix = 'azim_df_search_event_'+search_step+'.txt'
    out_file = join_path([project.data_dir,file_suffix])
    if keyword_set(reset) then file_delete, out_file, /allow_nonexistent
    if file_test(out_file) eq 0 then begin
        candidates = list()
        candidates.add, azim_df_search_pre_midn_events(project=project, search_step=search_step), /extract
        candidates.add, azim_df_search_post_midn_events(project=project, search_step=search_step), /extract

        times = dblarr(candidates.length)
        foreach candidate, candidates, ii do times[ii] = candidate.time_range[0]
        index = sort(times)
        candidates = candidates[index]

        foreach candidate, candidates, ii do candidate.id = ii+1

        ftouch, out_file
        foreach candidate, candidates do azim_df_subgroup_write_file, candidate, filename=out_file
    endif

    lprmsg, ''
    lprmsg, 'Read event from file ...'
    return, azim_df_subgroup_read_file(out_file)
end