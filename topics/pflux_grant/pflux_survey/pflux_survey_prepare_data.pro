;+
; Prepare data to be used for binning.
;-

pro pflux_survey_prepare_data_time, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range

    if cdf_has_var(time_var, filename=cdf_file) then return

    common_time_step = 60. ; sec.
    project = pflux_survey_load_project()
    common_times = make_bins(mission_time_range, common_time_step)
    settings = dictionary($
        'unit', 'sec', $
        'time_step', common_time_step, $
        'mission_time_range', mission_time_range, $
        'time_var_type', 'unix' )
    cdf_save_var, time_var, value=common_times, filename=cdf_file
    cdf_save_setting, settings, filename=cdf_file, varname=time_var

end


pro pflux_survey_prepare_data_omni, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range, var_name=the_vars

    has_var = 1
    foreach the_var, the_vars do if ~cdf_has_var(the_var, filename=cdf_file) then has_var = 0
    if has_var then return

    omni_read_index, mission_time_range, resolution='1min'
    pflux_survey_prepare_data_time, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range
    common_times = cdf_read_var(time_var)
    foreach the_var, the_vars do begin
        interp_time, the_var, common_times

        data = float(get_var_data(the_var))
        settings = dictionary($
            'depend_0', time_var, $
            'display_type', 'scalar', $
            'short_name', strupcase(the_var), $
            'unit', 'nT')
        cdf_save_var, the_var, value=data, filename=cdf_file
        cdf_save_setting, settings, filename=cdf_file, varname=the_var
    endforeach

end


pro pflux_survey_prepare_data_r_gse, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range, var_name=the_var

    if cdf_has_var(the_var, filename=cdf_file) then return

    rbsp_read_orbit, mission_time_range, probe=probe
    pflux_survey_prepare_data_time, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range
    common_times = cdf_read_var(time_var)
    interp_time, the_var, common_times

    data = float(get_var_data(the_var))
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'R', $
        'unit', 'Re', $
        'coord', 'GSE', $
        'coord_labels', ['x','y','z'])
    cdf_save_var, the_var, value=data, filename=cdf_file
    cdf_save_setting, settings, filename=cdf_file, varname=the_var

end


pro pflux_survey_prepare_data_r_sm, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range, var_name=the_var

    if cdf_has_var(the_var, filename=cdf_file) then return

    pflux_survey_prepare_data_time, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range
    common_times = cdf_read_var(time_var)
    
    r_gse_var = 'rbsp'+probe+'_r_gse'
    pflux_survey_prepare_data_r_gse, cdf_file, time_var=time_var, time_range=mission_time_range, probe=probe, var_name=r_gse_var
    r_gse = cdf_read_var(r_gse_var)
    
    data = float(cotran(r_gse, common_times, 'gse2sm'))
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'R', $
        'unit', 'Re', $
        'coord', 'SM', $
        'coord_labels', ['x','y','z'])
    cdf_save_var, the_var, value=data, filename=cdf_file
    cdf_save_setting, settings, filename=cdf_file, varname=the_var

end


pro pflux_survey_prepare_data_b_gsm, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range, var_name=the_vars

    has_var = 1
    foreach the_var, the_vars do if ~cdf_has_var(the_var, filename=cdf_file) then has_var = 0
    if has_var then return

    pflux_survey_prepare_data_time, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range
    common_times = cdf_read_var(time_var)
    ncommon_time = n_elements(common_times)
    ndim = 3
    secofday = constant('secofday')
    days = make_bins(mission_time_range+[0,-1]*secofday, secofday)
    foreach the_var, the_vars do begin
        data_type = (strsplit(the_var,'_',/extract))[1]
        b_gsm = fltarr(ncommon_time,ndim)
        foreach day, days do begin
            day_time_range = day+[0,secofday]
            pflux_grant_read_bfield, day_time_range, probe=probe, id=data_type
            time_index = where_pro(common_times, '[]', day_time_range)
            times = common_times[time_index]
            interp_time, the_var, times
            b_gsm[time_index,*] = get_var_data(the_var)
        endforeach

        data = float(temporary(b_gsm))
        settings = dictionary($
            'depend_0', time_var, $
            'display_type', 'vector', $
            'short_name', strupcase(data_type), $
            'unit', 'nT', $
            'coord', 'GSM', $
            'coord_labels', ['x','y','z'])
        cdf_save_var, the_var, value=data, filename=cdf_file
        cdf_save_setting, settings, filename=cdf_file, varname=the_var
    endforeach

end


pro pflux_survey_prepare_data_b_sm, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range, var_name=the_vars

    has_var = 1
    foreach the_var, the_vars do if ~cdf_has_var(the_var, filename=cdf_file) then has_var = 0
    if has_var then return

    pflux_survey_prepare_data_time, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range
    common_times = cdf_read_var(time_var)
    b_gsm_vars = 'rbsp'+probe+'_b_gsm'
    pflux_survey_prepare_data_b_gsm, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range, var_name=b_gsm_vars
    foreach the_var, the_vars, ii do begin
        data_type = (strsplit(the_var,'_',/extract))[1]
        b_gsm = get_var_data(b_gsm_vars[ii])
        
        data = float(cotran(b_gsm, common_times, 'gsm2sm'))
        settings = dictionary($
            'depend_0', time_var, $
            'display_type', 'vector', $
            'short_name', strupcase(data_type), $
            'unit', 'nT', $
            'coord', 'SM', $
            'coord_labels', ['x','y','z'])
        cdf_save_var, the_var, value=data, filename=cdf_file
        cdf_save_setting, settings, filename=cdf_file, varname=the_var
    endforeach

end


pro pflux_survey_prepare_ion_t, cdf_file, probe=probe, time_var=time_var, time_range=mission_time_range, var_name=the_vars

    ;;---Prepare ion_t.
    ;ion_t_var = prefix+'ion_t'
    ;if ~cdf_has_var(ion_t_var, filename=cdf_file) then begin
    ;the_var = ion_t_var
    ;rbsp_read_ion_temp, mission_time_range, probe=probe
    ;ion_t = get_var_data(ion_t_var, at=common_times)
    ;get_data, ion_t_var, the_times
    ;data_time_range = minmax(the_times)
    ;data_time_range = data_time_range-(data_time_range mod common_time_step)+[common_time_step,0]
    ;index = where_pro(common_times, '][', data_time_range, count=count)
    ;if count ne 0 then ion_t[index] = !values.f_nan
    ;index = where(ion_t lt 0, count)
    ;if count ne 0 then ion_t[index] = !values.f_nan
    ;data = float(temporary(ion_t))
    ;settings = dictionary($
    ;'depend_0', time_var, $
    ;'display_type', 'scalar', $
    ;'short_name', 'T!Dp!N', $
    ;'unit', 'eV')
    ;cdf_save_var, the_var, value=data, filename=cdf_file
    ;cdf_save_setting, settings, filename=cdf_file, varname=the_var
    ;endif
end



pro pflux_survey_prepare_data, cdf_file, probe=probe

    project = pflux_survey_load_project()
    probe_infos = project.probe_infos

    if n_elements(cdf_file) eq 0 then $
        cdf_file = join_path([project.data_dir,project.data_file_suffix])
    if file_test(cdf_file) eq 0 then begin
        path = fgetpath(cdf_file)
        if ~file_test(path,/directory) then file_mkdir, path
        cdf_id = cdf_create(cdf_file)
        cdf_close, cdf_id
    endif


    ; Settings.
    mission_time_range = project.time_range

    ; Time.
    time_var = 'time'
    pflux_survey_prepare_data_time, cdf_file, time_var=time_var, time_range=mission_time_range, probe=probe

    ; Dst and AE.
    dst_var = 'dst'
    ae_var = 'ae'
    omni_vars = [dst_var,ae_var]
    pflux_survey_prepare_data_omni, cdf_file, time_var=time_var, time_range=mission_time_range, probe=probe, var_name=omni_vars

    ; Orbit.
    prefix = 'rbsp'+probe+'_'
    r_sm_var = prefix+'r_sm'
    pflux_survey_prepare_data_r_sm, cdf_file, time_var=time_var, time_range=mission_time_range, probe=probe, var_name=r_sm_var


    ; B0 and B1 GSM.
    b_sm_vars = prefix+['b0','b1']+'_sm'
    pflux_survey_prepare_data_b_sm, cdf_file, time_var=time_var, time_range=mission_time_range, probe=probe, var_name=b_sm_vars


    ; T ion.
;
;;---Prepare q_sm2fac.
    ;q_var = prefix+'q_sm2fac'
    ;if ~cdf_has_var(q_var, filename=cdf_file) then begin
        ;the_var = q_var
        ;cdf_load_var, r_sm_var, filename=cdf_file
        ;cdf_load_var, b1_sm_var, filename=cdf_file
        ;define_fac, b1_sm_var, r_sm_var
        ;data = get_var_data(the_var)
        ;settings = dictionary($
            ;'depend_0', time_var, $
            ;'in_coord', coord, $
            ;'in_coord_labels', xyz, $
            ;'out_coord', get_setting(the_var, 'out_coord'), $
            ;'out_coord_labels', get_setting(the_var, 'out_coord_labels'))
        ;cdf_save_var, the_var, value=data, filename=cdf_file
        ;cdf_save_setting, settings, filename=cdf_file, varname=the_var
    ;endif
;
;;---Prepare pflux, pflux_dot0.
    ;pf_vars = prefix+['pf_fac_norm','pfdot0_fac_norm']
    ;ndim = 3
    ;ncommon_time = n_elements(common_times)
    ;if ~cdf_has_var(pf_vars[0], filename=cdf_file) then begin
    ;;---Break down to smaller time range, because the pflux sample rate is high.
        ;time_range_cadence = 'day'
        ;time_ranges = break_down_times(mission_time_range, time_range_cadence)
        ;time_ranges = sort_uniq([time_ranges,mission_time_range])
        ;ntime_range = n_elements(time_ranges)-1
;
        ;data_ptr = ptrarr(n_elements(pf_vars))
        ;foreach the_var, pf_vars, var_id do data_ptr[var_id] = ptr_new(fltarr(ncommon_time,ndim))
;
        ;for ii=0, ntime_range-1 do begin
            ;current_time_range = time_ranges[ii:ii+1]
            ;time_index = where_pro(common_times, '[)', current_time_range, count=ntime)
            ;if ntime eq 0 then continue
;
            ;pflux_grant_read_preprocessed_pflux, current_time_range, probe=probe
;
            ;foreach the_var, pf_vars, var_id do begin
                ;; Downsample: Tried interpolation and average.
                ;; The latter is better, b/c averaging preserves energy.
                ;;the_data = get_var_data(the_var, at=common_times[time_index])
                ;the_data = fltarr(ntime,ndim)
                ;get_data, the_var, full_times, data
                ;times = common_times[time_index]
                ;for jj=0, ntime-1 do begin
                    ;index = where_pro(full_times, '[)', times[jj]+[0,common_time_step], count=count)
                    ;the_data[jj,*] = total(data[index,*],1, /nan)/count
                ;endfor
;
                ;; Change earthward to parallel.
                ;if check_if_update(r_sm_var) then cdf_load_var, r_sm_var, filename=cdf_file
                ;r_sm = get_var_data(r_sm_var, in=current_time_range)
                ;index = where(r_sm[*,2] lt 0, count)
                ;if count ne 0 then the_data[index,0] *= -1
;
                ;(*data_ptr[var_id])[time_index,*] = the_data
            ;endforeach
        ;endfor
;
        ;foreach the_var, pf_vars, var_id do begin
            ;data = temporary(*data_ptr[var_id])
            ;ptr_free, data_ptr[var_id]
            ;if n_elements(data)/ndim ne n_elements(common_times) then message, 'Inconsistent data ...'
            ;settings = dictionary($
                ;'depend_0', time_var, $
                ;'display_type', 'vector', $
                ;'short_name', 'S', $
                ;'unit', 'mW/m!U2!N', $
                ;'coord', 'FAC', $
                ;'coord_labels', ['b','w','o'])
            ;cdf_save_var, the_var, value=data, filename=cdf_file
            ;cdf_save_setting, settings, filename=cdf_file, varname=the_var
        ;endforeach
    ;endif
;
;;---Prepare bw_ratio.
    ;bw_var = prefix+'bw_ratio'
    ;if ~cdf_has_var(bw_var, filename=cdf_file) then begin
    ;;---Break down to smaller time range, because the pflux sample rate is high.
        ;time_range_cadence = 'day'
        ;time_ranges = break_down_times(mission_time_range, time_range_cadence)
        ;time_ranges = sort_uniq([time_ranges,mission_time_range])
        ;ntime_range = n_elements(time_ranges)-1
;
        ;the_var = bw_var
        ;data = fltarr(ncommon_time)
;
        ;for ii=0, ntime_range-1 do begin
            ;current_time_range = time_ranges[ii:ii+1]
            ;time_index = where_pro(common_times, '[)', current_time_range, count=ntime)
            ;if ntime eq 0 then continue
;
            ;pflux_grant_read_preprocessed_ebfield, current_time_range, probe=probe, id='bw_ratio'
;
            ;; Downsample. Interpolate is fine here, do not need to be very accurate.
            ;data[time_index] = get_var_data(the_var, at=common_times[time_index])
        ;endfor
;
        ;settings = dictionary($
            ;'depend_0', time_var, $
            ;'display_type', 'scalar', $
            ;'short_name', 'Bw ratio', $
            ;'unit', '#')
        ;cdf_save_var, the_var, value=data, filename=cdf_file
        ;cdf_save_setting, settings, filename=cdf_file, varname=the_var
    ;endif
;
;;---Prepare de, de_dot0.
    ;de_vars = prefix+['de_sm','dedot0_sm']
    ;de_gsm_vars = prefix+['e_gsm','edot0_gsm']
    ;if ~cdf_has_var(de_vars[0], filename=cdf_file) then begin
    ;;---Break down to smaller time range, because the pflux sample rate is high.
        ;time_range_cadence = 'day'
        ;time_ranges = break_down_times(mission_time_range, time_range_cadence)
        ;time_ranges = sort_uniq([time_ranges,mission_time_range])
        ;ntime_range = n_elements(time_ranges)-1
;
        ;data_ptr = ptrarr(n_elements(de_vars))
        ;foreach the_var, de_vars, var_id do data_ptr[var_id] = ptr_new(fltarr(ncommon_time,ndim))
;
        ;for ii=0, ntime_range-1 do begin
            ;current_time_range = time_ranges[ii:ii+1]
            ;time_index = where_pro(common_times, '[)', current_time_range, count=ntime)
            ;if ntime eq 0 then continue
;
            ;pflux_grant_read_preprocessed_ebfield, current_time_range, probe=probe, coord='gsm'
;
            ;foreach the_var, de_gsm_vars, var_id do begin
                ;; Downsample.
                ;the_data = fltarr(ntime,ndim)
                ;get_data, the_var, full_times, data
                ;times = common_times[time_index]
                ;for jj=0, ntime-1 do begin
                    ;index = where_pro(full_times, '[)', times[jj]+[0,common_time_step], count=count)
                    ;the_data[jj,*] = total(data[index,*],1, /nan)/count
                ;endfor
;
                ;(*data_ptr[var_id])[time_index,*] = the_data
            ;endforeach
        ;endfor
;
        ;foreach the_var, de_vars, var_id do begin
            ;data = temporary(*data_ptr[var_id])
            ;ptr_free, data_ptr[var_id]
            ;if n_elements(data)/ndim ne n_elements(common_times) then message, 'Inconsistent data ...'
            ;data = float(cotran(data, common_times, 'gsm2sm'))
            ;settings = dictionary($
                ;'depend_0', time_var, $
                ;'display_type', 'vector', $
                ;'short_name', 'E', $
                ;'unit', 'mV/m', $
                ;'coord', coord, $
                ;'coord_labels', xyz)
            ;cdf_save_var, the_var, value=data, filename=cdf_file
            ;cdf_save_setting, settings, filename=cdf_file, varname=the_var
        ;endforeach
    ;endif


end

foreach probe, ['a','b'] do pflux_survey_prepare_data, probe=probe
end
