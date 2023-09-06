;+
; The E dot0 correction coef around dawn is larger than dusk. Check why.
;-


;---Settings.
    probe = 'a'
    time_range = time_double(['2012-10-01','2015-10-01'])
    ;time_range = time_double(['2012-10-01','2013-10-01'])
    prefix = 'rbsp'+probe+'_'

    ; Read orbit.
    r_gse_var = prefix+'r_gse'
    if check_if_update(r_gse_var, time_range) then begin
        rbsp_read_orbit, time_range, probe=probe
    endif
    apogee_mlt_var = prefix+'apogee_mlt'
    if check_if_update(apogee_mlt_var) then begin
        r_gse = get_var_data(r_gse_var, times=times)
        dis = snorm(r_gse)
        index = where(dis ge max(dis)*0.8)
        apogee_time_ranges = times[time_to_range(index,time_step=1)]
        napogee = n_elements(apogee_time_ranges)*0.5
        apogee_times = dblarr(napogee)
        for ii=0,napogee-1 do begin
            apogee_times[ii] = mean(apogee_time_ranges[ii,*])
        endfor
        apogee_r_gse = sinterpol(r_gse, times, apogee_times)
        apogee_mlt = pseudo_mlt(apogee_r_gse)
        store_data, apogee_mlt_var, apogee_times, apogee_mlt
    endif
    
    test_mlt_range = [5,6]
    apogee_mlt = get_var_data(apogee_mlt_var, times=apogee_times)
    index = where_pro(apogee_mlt, '[]', test_mlt_range)
    test_time_ranges = apogee_times[time_to_range(index,time_step=1)]
    ntest_section = n_elements(test_time_ranges)*0.5
    for section_id=0,ntest_section-1 do begin
        the_time_range = test_time_ranges[section_id,*]
        
        e_var = prefix+'e_mgse'
        flag_var = prefix+'efw_flags'
        ew_var = prefix+'ew_dot0'
        br_var = prefix+'bw_ratio'
        if check_if_update(e_var, the_time_range) then pflux_grant_read_e_mgse, the_time_range, probe=probe
        if check_if_update(flag_var, the_time_range) then rbsp_efw_read_flags, the_time_range, probe=probe
        if check_if_update(ew_var, the_time_range) then pflux_grant_read_ew_dot0, the_time_range, probe=probe
        
        get_data, e_var, common_times, e_mgse, limits=lim
        ew = get_var_data(ew_var)
        e_mgse[*,0] = ew
        edot0_var = prefix+'edot0_mgse'
        store_data, edot0_var, common_times, e_mgse, limits=lim

        bad_index = []

        ; B ratio for edot0        
        min_br = 0.2
        get_data, br_var, common_times, br
        bad_index = [bad_index, where(abs(br) le min_br)]
        
        ; Flags for general bad data.
        flags = get_var_data(flag_var, times=times, limits=lim)
        boom_strs = ['1','2','3','4']
        wanted_flags = ['eclipse','maneuver','efw_sweep',$
            'v'+boom_strs+'_saturation','boomflag'+boom_strs]
        nwanted_flag = n_elements(wanted_flags)
        wanted_flag_index = fltarr(nwanted_flag)
        for ii=0,nwanted_flag-1 do wanted_flag_index[ii] = where(lim.labels eq wanted_flags[ii])
        the_flags = flags[*,wanted_flag_index]
        overall_flags = total(the_flags, 2) ne 0
        ; Extend by 10 min.
        pad_time = 600.
        time_step = sdatarate(times)
        pad_index = pad_time/time_step
        ntime = n_elements(times)
        index = where(overall_flags eq 1, count)
        bad_sections = times[time_to_range(index,time_step=1)]
        the_bad_index = (bad_sections-times[0])/time_step
        the_bad_index[*,0] = (the_bad_index[*,0]-pad_index)>0
        the_bad_index[*,1] = (the_bad_index[*,1]+pad_index)<(ntime-1)
        nbad_section = n_elements(bad_sections)*0.5
        for ii=0,nbad_section-1 do begin
            i0 = the_bad_index[ii,0]
            i1 = the_bad_index[ii,1]
            overall_flags[i0:i1] = 1
        endfor
        ; Interpolate to common_times.
        overall_flags = sinterpol(overall_flags, times, common_times)
        bad_index = [bad_index, where(overall_flags eq 1)]
        ;bad_index = where(overall_flags eq 1)
        store_data, prefix+'flag', times, the_flags, limits={$
            yrange:[-0.2,1.2], labels:wanted_flags}
        
        ; Apply.
        bad_index = sort_uniq(bad_index)
        foreach var, [e_var,ew_var,edot0_var] do begin
            get_data, var, times, data
            data[bad_index,*] = !values.f_nan
            store_data, var+'_corr', times, data
        endforeach
        
        
        ; Rotate into FAC.
        b0_gsm_var = prefix+'b0_gsm'
        if check_if_update(b0_gsm_var, the_time_range) then pflux_grant_read_bfield, the_time_range, probe=probe
        rbsp_read_orbit, time_range, probe=probe
        get_data, b0_gsm_var, times
        r_gse = get_var_data(r_gse_var, at=times)
        r_gsm = cotran(r_gse, times, 'gse2gsm')
        store_data, prefix+'r_gsm', times, r_gsm
        add_setting, prefix+'r_gsm', /smart, dictionary($
            'diaplay_type', 'vector', $
            'short_name', 'R', $
            'unit', 'Re', $
            'coord', 'GSM', $
            'coord_labels', xyz )        
        define_fac, b0_gsm_var, prefix+'r_gsm', time_var=prefix+'r_gsm'
        q_var = prefix+'q_gsm2fac'
        get_data, q_var, times, q_gsm2fac
        store_data, q_var+'_lowres', times, q_gsm2fac
        q_gsm2fac = qslerp(q_gsm2fac, times, common_times)
        store_data, q_var, common_times, q_gsm2fac

        foreach var, prefix+['e','edot0'] do begin
            vec = get_var_data(var+'_mgse_corr')
            vec = cotran(vec, common_times, 'mgse2gsm', probe=probe)
            store_data, var+'_gsm', common_times, vec
            add_setting, var+'_gsm', /smart, dictionary($
                'display_type', 'vector', $
                'short_name', 'E', $
                'unit', 'mV/m', $
                'coord', 'GSM', $
                'coord_labels', xyz )
        endforeach
        stop
    endfor
stplot2cdf
end
