
function azim_df_search_around_midn_events, search_step=search_step, project=project, $
    reset=reset, test_time=test_time

    project = azim_df_load_project()
    search_name = 'beyond_15Re'
    region_name = 'around_midn'

    the_prefix = 'azim_df_search_event_'+region_name+'_'+search_name+'_'
    search_setting = dictionary($
        'name', region_name, $
        'probes', 'th'+letters('e'), $
        'time_range', time_double(['2007-03-01','2009-09-01']), $
        'data_file_suffix', 'azim_df_primitive_data_themis.cdf')

;---Search times within ROI.
    mlt_limit = 3.  ; hour.
    mlt_range = [-1,1]*mlt_limit
    rxy_range = [4.,30] ; Re.
    roi_min_count = 4
    pdyn = 10.  ; nPa, ~8 Re at nose.
    search_roi = dictionary($
        'name', search_name, $
        'mlt_range', mlt_range, $
        'rxy_range', rxy_range, $
        'pdyn', pdyn, $
        'roi_min_count', roi_min_count, $
        'roi_min_duration', constant('secofhour'), $
        'roi_probes', !null, $
        'file_suffix', the_prefix+'search_roi.txt' )
    search_setting['search_roi'] = search_roi


;---Search times with enough triads.
    small_angle_limit = 15.
    triad_angle_range = [0,180]+[1,-1]*small_angle_limit
    triad_min_count = 2
    search_triad = dictionary($
        'name', search_name, $
        'small_angle_limit', small_angle_limit, $
        'triad_angle_range', triad_angle_range, $
        'triad_min_count', triad_min_count, $
        'triad_min_duration', constant('secofhour'), $
        'file_suffix', the_prefix+'search_triad.txt')
    search_setting['search_triad'] = search_triad


;---Search DF within the candidate periods.
    search_df = dictionary($
        'name', search_name, $
        'file_suffix', the_prefix+'search_df.txt')
    search_setting['search_df'] = search_df


;---Search large DF within the candidate periods.
    df_height_range = [1.,180]
    df_width_range = [60.,1800]
    large_df_min_height = 10. ; deg.
    large_df_min_count = 4.
    scale_width = 3.    ; MLT.
    search_large_df = dictionary($
        'name', search_name, $
        'large_df_min_count', large_df_min_count, $
        'large_df_min_height', large_df_min_height, $
        'scale_width', scale_width, $
        'df_height_range', df_height_range, $
        'df_width_range', df_width_range, $
        'file_suffix', the_prefix+'search_large_df.txt')
    search_setting['search_large_df'] = search_large_df


;---Select subgroups.
    min_consecutive_df_count = 4
    search_subgroup = dictionary($
        'name', search_name, $
        'min_consecutive_df_count', min_consecutive_df_count, $
        'file_suffix', the_prefix+'search_subgroup.txt')
    search_setting['search_subgroup'] = search_subgroup


;---Select DF groups.
    search_df_group = dictionary($
        'name', search_name, $
        'xcor_section_ratio', [-1,1]*1.5, $
        'xcor_min_value', 0.5, $
        'obs_time_diff_range', [60.,1800], $
        'max_time_lag_ratio', 0.3, $
        'file_suffix', the_prefix+'search_df_group.txt')
    search_setting['search_df_group'] = search_df_group


;---Select uniq DF.
    search_uniq_subgroup = dictionary($
        'name', search_name, $
        'file_suffix', the_prefix+'search_uniq_subgroup.txt')
    search_setting['search_uniq_subgroup'] = search_uniq_subgroup



    if n_elements(search_step) eq 0 then search_step = 'uniq_subgroup'
    routine = 'azim_df_search_'+search_step
    candidates = call_function(routine, search_setting, project=project, reset=reset, test_time=test_time)
    return, candidates
end

reset = 0
search_steps = ['uniq_subgroup']
foreach search_step, search_steps do $
    events = azim_df_search_around_midn_events(search_step=search_step, reset=reset)
end
