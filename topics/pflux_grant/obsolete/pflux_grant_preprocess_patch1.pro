;+
; Calculate Ew using the E dot B = 0 assumption.
; Save B0, which is B-B1.
;-

pro pflux_grant_preprocess_patch1, probe=probe, project=project

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project.pflux_calc_setting

    if n_elements(probe) eq 0 then return
    if (where(probe eq ['a','b']))[0] eq -1 then return
    rbspx = 'rbsp'+probe
    mission_time_range = pflux_calc_setting[rbspx].time_range
    secofday = constant('secofday')
    nday = total(mission_time_range*[-1,1])/secofday
    dates = mission_time_range[0]+findgen(nday)*secofday
    dates = (probe eq 'a')? time_double(['2012-09-30','2015-10-01']): time_double(['2018-04-01'])

;---Some settings.
    common_time_step = project.common_time_step
    orbit_time_step = 60.
    xyz = constant('xyz')
    prefix = rbspx+'_'

    local_root = join_path([default_local_root(),'sdata','rbsp'])
    paths = [rbspx,'ebfield']
    time_var = 'ut_sec'
    ; R related var.
    r_gsm_var = prefix+'r_gsm'
    orbit_time_var = 'ut_orbit'
    q_var = prefix+'q_uvw2gsm'
    orbit_time_setting = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix')
    ; B0 related var.
    b_gsm_var = prefix+'b_gsm'
    b1_gsm_var = prefix+'b1_gsm'
    db_gsm_var = prefix+'b1_gsm'
    b0_gsm_var = prefix+'b0_gsm'
    b0_setting = dictionary($
        'depend_0', orbit_time_var, $
        'display_type', 'vector', $
        'short_name', 'B0', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz)
    ; Ew related var.
    e_uv_var = prefix+'e_uv'
    ew_var = prefix+'ew'
    ew_setting = dictionary($
        'depend_0', time_var, $
        'display_type', 'scalar', $
        'unit', 'mV/m', $
        'short_name', 'E!Dw!N')
    bw_ratio_var = prefix+'bw_ratio'
    bw_ratio_setting = dictionary($
        'depend_0', time_var, $
        'display_type', 'scalar', $
        'unit', '#', $
        'short_name', 'B!Dw!N/|B|')
    foreach date, dates do begin
        time_range = date+[0,secofday]
        base_name = rbspx+'_ebfield_'+time_string(date,tformat='YYYY_MMDD')+'_v01.cdf'
        year = time_string(date,tformat='YYYY')
        data_file = join_path([local_root,paths,year,base_name])
        if file_test(data_file) eq 0 then continue
        in_id = cdf_open(data_file)

    ;---Orbit time.
        cdf_load_var, b1_gsm_var, filename=in_id
        b1_gsm = get_var_data(b1_gsm_var, times=common_times)
        norbit_time = secofday/orbit_time_step
        orbit_time_index = smkarthm(0,orbit_time_step/common_time_step,norbit_time,'x0')
        orbit_times = common_times[orbit_time_index]
        if ~cdf_has_var(orbit_time_var, filename=in_id) then begin
            cdf_save_var, orbit_time_var, value=orbit_times, filename=in_id
            cdf_save_setting, orbit_time_setting, filename=in_id, varname=orbit_time_var
        endif

    ;---B0.
        if ~cdf_has_var(b0_gsm_var, filename=in_id) then begin
            rbsp_read_bfield, time_range, probe=probe, resolution='4sec'
            b_gsm = get_var_data(b_gsm_var, at=common_times)
            b0_gsm = b_gsm-b1_gsm
            b0_gsm = b0_gsm[orbit_time_index,*]
            b0_gsm = float(b0_gsm)
            cdf_save_var, b0_gsm_var, value=b0_gsm, filename=in_id
            cdf_save_setting, b0_setting, filename=in_id, varname=b0_gsm_var
        endif

    ;---Ew and Bw/|B|. Saved separately, to cut bad data with different thresholds.
        if ~cdf_has_var(ew_var, filename=in_id) or ~cdf_has_var(bw_ratio_var, filename=in_id) then begin
            rbsp_read_quaternion, time_range, probe=probe
            cdf_load_var, b0_gsm_var, filename=in_id
            interp_time, b0_gsm_var, common_times
            b0_uvw_var = prefix+'b0_uvw'
            rbsp_gsm2uvw, b0_gsm_var, b0_uvw_var, quaternion=q_var, probe=probe
            b0_uvw = get_var_data(b0_uvw_var)
            cdf_load_var, e_uv_var, filename=in_id
            e_uv = get_var_data(e_uv_var)
            ew = -(e_uv[*,0]*b0_uvw[*,0]+e_uv[*,1]+b0_uvw[*,1])/b0_uvw[*,2]
            ew = float(ew)
            bw_ratio = float(b0_uvw[*,2]/snorm(b0_uvw))
            if ~cdf_has_var(ew_var, filename=in_id) then begin
                cdf_save_var, ew_var, value=ew, filename=in_id
                cdf_save_setting, ew_setting, filename=in_id, varname=ew_var
            endif
            if ~cdf_has_var(bw_ratio_var) then begin
                cdf_save_var, bw_ratio_var, value=bw_ratio, filename=in_id
                cdf_save_setting, bw_ratio_setting, filename=in_id, varname=bw_ratio_var
            endif
        endif
        cdf_close, in_id
    endforeach

end

foreach probe, ['b'] do pflux_grant_preprocess_patch1, probe=probe, project=project
end