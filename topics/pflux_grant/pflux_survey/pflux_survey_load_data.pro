;+
; Load data from prepared data.
;-

pro pflux_survey_load_mlt, probe=probe, data_file=data_file

    var_type = 'r_sm'
    pflux_survey_load_data, var_type, probe=probe

    prefix = 'rbsp'+probe+'_'
    var_name = prefix+var_type
    r_sm = get_var_data(var_name, times=times)
    data = pseudo_mlt(r_sm)
    the_var = prefix+'mlt'
    store_data, the_var, times, data
    add_setting, the_var, /smart, dictionary($
        'display_type', 'scalar', $
        'unit', 'hr', $
        'short_name', 'MLT')

end


pro pflux_survey_load_mlat, probe=probe, data_file=data_file

    var_type = 'r_sm'
    pflux_survey_load_data, var_type, probe=probe

    prefix = 'rbsp'+probe+'_'
    var_name = prefix+var_type
    r_sm = get_var_data(var_name, times=times)
    data = asin(r_sm[*,2]/snorm(r_sm))*constant('deg')
    the_var = prefix+'mlat'
    store_data, the_var, times, data
    add_setting, the_var, /smart, dictionary($
        'display_type', 'scalar', $
        'unit', 'deg', $
        'short_name', 'MLat')

end


pro pflux_survey_load_dis, probe=probe, data_file=data_file

    var_type = 'r_sm'
    pflux_survey_load_data, var_type, probe=probe

    prefix = 'rbsp'+probe+'_'
    var_name = prefix+var_type
    r_sm = get_var_data(var_name, times=times)
    data = snorm(r_sm)
    the_var = prefix+'dis'
    store_data, the_var, times, data
    add_setting, the_var, /smart, dictionary($
        'display_type', 'scalar', $
        'unit', 'Re', $
        'short_name', 'R')

end


pro pflux_survey_load_pf_sm_norm, probe=probe, data_file=data_file

    var_types = ['pf_fac_norm','q_sm2fac']
    foreach var_type, var_types do pflux_survey_load_data, var_type, probe=probe

    prefix = 'rbsp'+probe+'_'
    the_var = prefix+'pf_sm_norm'
    from_fac, prefix+'pf_fac_norm', to=the_var, q_var=prefix+'q_sm2fac'
    add_setting, the_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'mW/m!U2!N', $
        'short_name', 'S', $
        'coord', 'SM', $
        'coord_labels', constant('xyz'))

end


pro pflux_survey_load_pfdot0_sm_norm, probe=probe, data_file=data_file

    var_types = ['pfdot0_fac_norm','q_sm2fac']
    foreach var_type, var_types do pflux_survey_load_data, var_type, probe=probe

    prefix = 'rbsp'+probe+'_'
    the_var = prefix+'pfdot0_sm_norm'
    from_fac, prefix+'pfdot0_fac_norm', to=the_var, q_var=prefix+'q_sm2fac'
    add_setting, the_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'mW/m!U2!N', $
        'short_name', 'S!S!Udot0!R!N', $
        'coord', 'SM', $
        'coord_labels', constant('xyz'))

end


pro pflux_survey_load_de_fac, probe=probe, data_file=data_file

    var_types = ['de_sm','q_sm2fac']
    foreach var_type, var_types do pflux_survey_load_data, var_type, probe=probe

    prefix = 'rbsp'+probe+'_'
    the_var = prefix+'de_fac'
    q_var = prefix+'q_sm2fac'
    to_fac, prefix+'de_sm', to=the_var, q_var=q_var
    add_setting, the_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'mV/m', $
        'short_name', 'E', $
        'coord', 'FAC', $
        'coord_labels', get_setting(q_var, 'out_coord_labels'))

end



pro pflux_survey_load_data, var_type, probe=probe

    if n_elements(var_type) eq 0 then message, 'No var_type ...'

    project = pflux_survey_load_project()
    data_file = join_path([project.data_dir,project.data_file_suffix])
    pflux_survey_prepare_data, data_file, probe=probe

    ; Default vars are: rbspx_ut,q_sm2fac,[r,b0,de,dedot0]_sm,[pf,pfdot0]_fac_norm
    the_var = (n_elements(probe) ne 0)? 'rbsp'+probe+'_'+var_type: var_type
    if cdf_has_var(the_var, filename=data_file) then begin
        cdf_load_var, the_var, filename=data_file
    endif else begin
        call_procedure, 'pflux_survey_load_'+var_type, probe=probe, data_file=data_file
    endelse

;
;;---Remove when dot0 is invalid.
;    fillval = !values.f_nan
;    vars = ['dedot0_sm','pfdot0_sm_norm']
;    index = where(vars eq var_type, count)
;    if count ne 0 then begin
;        min_bw_ratio = 0.2
;        bw_ratio_var = prefix+'bw_ratio'
;        cdf_load_var, bw_ratio_var, filename=data_file
;        bw_ratio = get_var_data(bw_ratio_var, times=times)
;
;        bad_flag = abs(bw_ratio) le min_bw_ratio
;        index = where(bad_flag eq 1, count)
;        time_step = total(times[0:1]*[-1,1])
;        pad_time = 1800d
;        pad_nrec = pad_time/time_step
;        ntime = n_elements(times)
;        for ii=0, count-1 do begin
;            i0 = (index[ii]-pad_nrec)>0
;            i1 = (index[ii]+pad_nrec)<(ntime-1)
;            bad_flag[i0:i1] = 1
;        endfor
;
;        bad_index = where(bad_flag eq 1, count)
;        if count ne 0 then begin
;            the_var = prefix+var_type
;            get_data, the_var, times, data
;            data[bad_index,*] = fillval
;            store_data, the_var, times, data
;        endif
;    endif
;
;
;;---Remove times of bad B0 data.
;    vars = ['b0_sm','pf_fac_norm','pfdot0_fac_norm','pf_sm_norm','pfdot0_sm_norm']
;    index = where(vars eq var_type, count)
;    if count ne 0 and probe eq 'b' then begin
;        ; RBSP-B.
;        bad_time_ranges = list()
;        bad_time_ranges.add, time_double(['2014-08-20','2014-08-28'])
;        bad_time_ranges.add, time_double(['2015-09-14','2015-09-17'])
;        the_var = prefix+var_type
;        get_data, the_var, times, data
;        foreach bad_time_range, bad_time_ranges do begin
;            index = lazy_where(times, '[]', bad_time_range)
;            data[index,*] = fillval
;        endforeach
;        store_data, the_var, times, data
;    endif


end
