;+
; A wrapper to return the final events.
; 
; project=.
; after_analysis=. Set to perform vel analysis.
;-

function azim_df_load_event, project=project, after_analysis=after_analysis
;    events = azim_df_search_df_group_mlt_select(project=project)
    events = azim_df_search_df_group_substorm(project=project)

    if ~keyword_set(after_analysis) then return, events

    data_file_suffix = project.name+'_all_analyzed_events.tplot'
    data_file = join_path([project.data_dir,data_file_suffix])
    the_var = 'all_analyzed_events'
    if file_test(data_file) eq 0 then begin
        foreach event, events do event = azim_df_analysis_event(event, project=project)
        store_data, the_var, 0, events
        tplot_save, the_var, filename=data_file
    endif
    tplot_restore, filename=data_file
    return, get_var_data(the_var)
end