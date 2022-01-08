;+
; To distinguish pflux and dipolarization.
;-

    test = 1
    ; Angelopoulos event.
    short_time_range = time_double(['2007-03-23/11:00','2007-03-23/12:00'])
    probes = 'th'+letters('e')
    ; Ogasawara Figure 3 event.
    short_time_range = time_double(['2008-02-26/03:50','2008-02-26/04:50'])
    probes = 'th'+['a','d','e']
    ; Ogasawara Figure 6 event.
;    short_time_range = time_double(['2008-03-03/07:50','2008-03-03/08:40'])
;    probes = 'th'+['d','e']

;---Check inputs.
    if n_elements(short_time_range) ne 2 then begin
        errmsg = handle_error('No input short_time_range ...')
        return
    endif

    if nprobe eq 0 then begin
        errmsg = handle_error('No input probes ...')
        return
    endif

;---Load data.
    mean_time = mean(short_time_range)
    event_id = time_string(mean_time,tformat='YYYY_MMDD_hhmm')
    common_data_rate = 10.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    local_root = project.data_dir
    data_file = join_path([local_root,'data','manual_event','fig_df_vs_pflux_'+event_id+'.tplot'])

    if file_test(data_file) eq 0 then begin
        vars = ['r_sm','theta']
        save_vars = list()
        foreach probe, probes do begin
            prefix = probe+'_'
            data_time_range = mean_time+[-1,1]*3600*5
            foreach var, vars do begin
                azim_df_read_data, var, time_range=data_time_range, probe=probe, project=project
                the_var = prefix+var
                save_vars.add, the_var
            endforeach
            
            
            the_probe = strmid(probe,0,1,/reverse)
            themis_read_bfield, data_time_range, probe=the_probe
            themis_read_efield, data_time_range, probe=the_probe
            stop

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
