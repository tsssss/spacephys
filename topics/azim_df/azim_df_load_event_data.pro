;+
; Load data based on given time range.
; Will search over all events to find the event that matches the given time range.
;-

pro azim_df_load_event_data, time_range, project=project, data_file=data_file, event=event

    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(event) eq 0 then begin
        event = azim_df_find_event(time_range, project=project)
        if isa(event, 'list') then event = events[0]
    endif

    data_time_range = time_range
    data_file_suffix = project.name+'_event_data_'+$
        strjoin(time_string(event.time_range+[1,-1]*constant('secofday'),tformat='YYYY_MMDD_hhmmss'),'_')+'_v01.tplot'
    if n_elements(data_file) eq 0 then data_file = join_path([project.data_dir,'event_data',data_file_suffix])
    if file_test(data_file) eq 0 then begin
        lprmsg, 'Load data for time range: '+strjoin(time_string(data_time_range),' to ')+' ...'
        del_data, '*'

        data_vars = ['theta','r_sm','b_sm']
        save_vars = list()
        time_step = project.time_step
        probe_infos = project.probe_infos
        probes = event.probes
        foreach probe, probes do begin
            lprmsg, 'Processing probe: '+probe+' ...'
            prefix = probe_infos[probe].prefix
            foreach var, data_vars do begin
                azim_df_read_data, var, time_range=data_time_range, probe=probe, project=project
                the_var = prefix+var
                uniform_time, the_var, time_step
            endforeach
            save_vars.add, prefix+data_vars, /extract

            mlt_var = prefix+'mlt'
            r_sms = get_var_data(prefix+'r_sm', times=times)
            mlt = azim_df_calc_pseudo_mlt(r_sms)
            store_data, mlt_var, times, mlt
            add_setting, mlt_var, /smart, {$
                unit: 'hr', $
                display_type: 'scalar', $
                short_name: 'MLT'}
            save_vars.add, mlt_var
        endforeach

        ; AE and Dst.
        omni_read_index, data_time_range
        save_vars.add, ['ae','dst'], /extract

        save_vars = save_vars.toarray()
        tplot_save, save_vars, filename=data_file
        lprmsg, 'Data saved to '+data_file+' ...'
    endif

    if file_test(data_file) then tplot_restore, filename=data_file

end
