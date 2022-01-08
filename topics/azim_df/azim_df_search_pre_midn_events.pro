
function azim_df_search_pre_midn_events, search_step=search_step, project=project, $
    reset=reset, test_time=test_time

    project = azim_df_load_project()

;---Some common settings.
    region_name = 'pre_midn'
    mlt_limit = 9.  ; hour.
    mlt_range = [-1,0]*mlt_limit
    search_settings = azim_df_search_event_settings()
    foreach search_setting, search_settings do begin
        search_setting.search_roi.mlt_range = mlt_range
        search_name = search_setting.name
        foreach key, search_setting.keys() do begin
            the_setting = search_setting[key]
            if ~isa(the_setting,'dictionary') then continue
            if ~the_setting.haskey('file_suffix') then continue
            the_prefix = 'azim_df_search_event_'+region_name+'_'+search_name+'_'
            the_setting['file_suffix'] = the_prefix+key+'.txt'
        endforeach
    endforeach


    if n_elements(search_step) eq 0 then search_step = 'uniq_subgroup'
    routine = 'azim_df_search_'+search_step
    candidates = list()
    foreach search_setting, search_settings do begin
        candidates.add, call_function(routine, search_setting, project=project, reset=reset, test_time=test_time), /extract
    endforeach

    return, candidates
end

events = azim_df_search_pre_midn_events(search_step=search_step, reset=reset)
end
