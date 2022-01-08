;+
; Search DF group containing large DF from different sc.
;-

function azim_df_search_df_group, search_setting, project=project, $
    reset=reset, test_time=test_time

    retval = list()
    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()

;---Settings for selecting subgroup.
    the_key = 'search_df_group'
    if ~search_setting.haskey(the_key) then message, 'No settings for DF group ...'
    search_info = search_setting[the_key]


;---Check if file exists, to avoid search again.
    file_suffix = search_info.file_suffix
    out_file = join_path([project.data_dir,file_suffix])
    log_file = join_path([project.data_dir,file_suffix+'.log'])
    if keyword_set(reset) then begin
        lprmsg, 'Resetting DF group search ...'
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


    ;---Settings.
        tmp = azim_df_filter_edge(project=project, log_file=log_file)


    ;---Load candidates.
        candidates = azim_df_search_subgroup(search_setting, project=project)
        foreach candidate, candidates do begin
            event = azim_df_filter_edge(candidate, project=project, log_file=log_file)
            if n_elements(event) eq 0 then continue

        ;---Output to file.
            event_id += 1
            event.id = event_id
            azim_df_subgroup_write_file, event, filename=out_file
        endforeach
    endif

    lprmsg, ''
    lprmsg, 'Read DF group from file ...'
    return, azim_df_subgroup_read_file(out_file)
end
