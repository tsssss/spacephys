;---Test events.
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



;---Load all df_group.
    project = azim_df_load_project()
    search_step = 'uniq_subgroup'
    candidates = azim_df_search_all_events(project=project, search_step=search_step)



;---Select test_list.
    test_candidates = list()
    foreach test_event_time, test_event_times do begin
        lprmsg, 'Find event '+time_string(test_event_time)+' ...'
        foreach candidate, candidates do begin
            if candidate.time_range[0] eq test_event_time then begin
                lprmsg, 'Found the wanted event ...'
                test_candidates.add, candidate
                break
            endif
        endforeach
    endforeach

    log_file = -1
    vmag_ratio_range = []
    vhat_angle_range = []
    foreach candidate, test_candidates do begin
        event = azim_df_filter_triad(candidate, project=project, log_file=log_file, vmag_ratios=vmag_ratios, vhat_angles=vhat_angles)
        vmag_ratio_range = [vmag_ratio_range,minmax(vmag_ratios)]
        vhat_angle_range = [vhat_angle_range,minmax(vhat_angles)]
        if n_elements(event) eq 0 then stop
    endforeach

    print, minmax(vmag_ratio_range)
    print, minmax(vhat_angle_range)

end
