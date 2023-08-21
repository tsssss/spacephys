;+
; Load the E and B fields to check if a "dipolarization" is due to Alfven wave.
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

    probes = 'th'+['d','e']
    short_time_range = time_double(['2008-03-28/06:25','2008-03-28/06:45'])

;---Check inputs.
    if n_elements(center_time) ne 0 then begin
        if size(center_time,/type) eq 7 then center_time = time_double(center_time)
        pad_time = constant('secofhour')*0.5
        short_time_range = center_time+[-1,1]*pad_time
    endif

    if n_elements(short_time_range) ne 2 then begin
        errmsg = handle_error('No input short_time_range ...')
        return
    endif

    nprobe = n_elements(probes)
    if nprobe eq 0 then begin
        errmsg = handle_error('No input probes ...')
        return
    endif

;---Load data.
    center_time = mean(short_time_range)
    event_id = time_string(center_time,tformat='YYYY_MMDD_hhmm')
    common_data_rate = 10.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    local_root = project.data_dir
    data_file = join_path([local_root,'data','manual_event','fig_df_vs_pflux_'+event_id+'.tplot'])

    if file_test(data_file) eq 0 then begin
        vars = ['r_sm','theta']
        save_vars = list()
        data_pad_time = constant('secofhour')*2
        time_step = 3.
        min_data_ratio = 0.8
        b0_window_size = 600./time_step     ; 10 min.
        ndim = 3
        foreach probe, probes do begin
            prefix = probe+'_'
            data_time_range = center_time+[-1,1]*data_pad_time

        ;---Load E and B field.
            the_probe = strmid(probe,0,1,/reverse)
            ntime = total(data_time_range*[-1,1])/time_step
            b_var = themis_read_bfield(data_time_range, probe=the_probe, coord='gsm')
            data = get_var_data(b_var)
            if n_elements(data)/3 le ntime*min_data_ratio then begin
                errmsg = handle_error('Not enough B data ...')
                return
            endif
            themis_read_efield, data_time_range, probe=the_probe
            e_var = prefix+'e_gsm'
            data = get_var_data(e_var)
            if n_elements(data)/3 le ntime*min_data_ratio then begin
                errmsg = handle_error('Not enough E data ...')
                return
            endif
            foreach var, [e_var,b_var] do begin
                uniform_time, var, time_step
                save_vars.add, var
            endforeach

        ;---Load position and theta.
            foreach var, vars do begin
                azim_df_read_data, var, time_range=data_time_range, probe=probe, project=project
                the_var = prefix+var
                save_vars.add, the_var
            endforeach
            mlt_var = prefix+'mlt'
            r_sms = get_var_data(prefix+'r_sm', times=times)
            mlt = azim_df_calc_pseudo_mlt(r_sms)
            store_data, mlt_var, times, mlt
            add_setting, mlt_var, /smart, {$
                unit: 'hr', $
                display_type: 'scalar', $
                short_name: 'MLT'}
            save_vars.add, mlt_var

        ;---Separate B0 and dB.
            b_gsm = get_var_data(b_var, times=times)
            db_gsm = b_gsm
            for ii=0, ndim-1 do begin
                tdat = b0_gsm[*,ii]
                db_gsm[*,ii] = tdat-smooth(tdat,b0_window_size, /edge_truncate, nan=1)
            endfor
            b0_gsm = b_gsm-db_gsm
            b0_var = prefix+'b0_gsm'
            store_data, b0_var, times, b0_gsm
            add_setting, b0_var, /smart, dictionary($
                'display_type', 'vector', $
                'short_name', 'B0', $
                'unit', 'nT', $
                'coord', 'GSM', $
                'coord_labels', ['x','y','z'] )

            db_var = prefix+'db_gsm'
            store_data, db_var, times, db_gsm
            add_setting, db_var, /smart, {$
                display_type: 'vector', $
                short_name: 'dB', $
                unit: 'nT', $
                coord: 'GSM', $
                coord_labels: ['x','y','z'] }

        ;---Convert dB and E to FAC.
            get_data, r_var, times, r_sms
            r_gsms = cotran(r_sms, times, 'sm2gsm')
            r_gsm_var = prefix+'r_gsm'
            store_data, r_gsm_var, times, r_gsms
            add_setting, r_gsm_var, smart=1, dictionary($
                'display_type', 'vector', $
                'short_name', 'R', $
                'unit', 'Re', $
                'coord', 'GSM', $
                'coord_labels', ['x','y','z'] )
            define_fac, b0_var, r_gsm_var
            e_fac_var = prefix+'e_fac'
            to_fac, e_var, to=e_fac_var
            db_fac_var = prefix+'db_fac'
            to_fac, db_var, to=db_fac_var
            

            pf_var = prefix+'pf_fac'
            stplot_calc_pflux_mor, e_fac_var, db_fac_var, pf_var
            add_setting, pf_var, smart=1, dictionary($
                'display_type', 'vector', $
                'short_name', 'S', $
                'unit', 'mW/m!U2!N', $
                'coord', 'FAC', $
                'coord_labels', ['b','w','o'] )
            save_vars.add, pf_var
            stop


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
