;+
; Read rbsp orbit vars: 
; [x,y,z]_sm, lshell, mlat, mlt
;-
pro ml_rbsp_load_orbit_var_gen_file, input_date, probe=probe, filename=out_file, resolution=time_step, downsample_method=downsample_method, _extra=extra

test = 0

    if n_elements(input_date) eq 0 then begin
        errmsg = 'No input time ...'
        return
    endif
    if n_elements(probe) eq 0 then begin
        errmsg = 'No input probe ...'
        return
    endif
    if n_elements(time_step) eq 0 then time_step = ml_time_step()
    ; can be mean, median, log_mean.
    if n_elements(downsample_method) eq 0 then downsample_method = 'interp'

    year_str = time_string(input_date[0],tformat='YYYY')
    year = float(year_str)
    time_range = time_double([year_str,string(year+1,format='(I04)')])
    secofday = constant('secofday')
if keyword_set(test) then time_range = time_range[0]+[0,30*secofday]
    days = make_bins(time_range+[0,-1]*secofday,secofday)
    common_times = make_bins(time_range+[0,-1]*time_step,time_step)+time_step*0.5

    prefix = 'rbsp'+probe+'_'

;---Init file.
    if file_test(out_file) eq 0 then cdf_touch, out_file

;---Add time.
    time_var = 'unix_time'
    if ~cdf_has_var(time_var, filename=out_file) then begin
        common_times = make_bins(time_range+[0,-1]*time_step, time_step)+time_step*0.5

        cdf_save_var, time_var, value=common_times, filename=out_file
        settings = dictionary($
            'FIELDNAM', 'Time', $
            'UNITS', 'sec', $
            'TEXT', 'center time of the sampling interval', $
            'VAR_TYPE', 'support_data' )
        cdf_save_setting, settings, var=time_var, filename=out_file
    endif
    common_times = cdf_read_var(time_var, filename=out_file)
    ntime = n_elements(common_times)

;---Add r_sm, lshell, mlt, mlat.
    fillval = !values.f_nan
    files = rbsp_load_spice(time_range, probe=probe)

    var_list = list()
    in_vars = prefix+['r_gse','mlt','mlat','lshell']
    var_list.add, dictionary($
        'in_vars', in_vars, $
        'time_var_name', 'Epoch', $
        'time_var_type', 'epoch' )
    read_vars, time_range, files=files, var_list=var_list
    foreach var, in_vars do interp_time, var, common_times

    ; r_sm.
    r_gse_var = prefix+'r_gse'
    r_gse = get_var_data(r_gse_var)
    r_sm = cotran(r_gse, common_times, 'gse2sm')
    settings = dictionary($
        'FIELDNAM', 'Spacecraft position in SM coord in Re', $
        'UNITS', 'Re', $
        'VAR_TYPE', 'data', $
        'DEPEND_0', time_var )
    r_sm_var = prefix+'r_sm'
    cdf_save_var, r_sm_var, value=r_sm, filename=out_file
    cdf_save_setting, settings, var=r_sm_var, filename=out_file

    ; mlt.
    the_var = prefix+'mlt'
    val = get_var_data(the_var)
    settings = dictionary($
        'FIELDNAM', 'Spacecraft MLT', $
        'UNITS', 'hr', $
        'VAR_TYPE', 'data', $
        'DEPEND_0', time_var )
    cdf_save_var, the_var, value=val, filename=out_file
    cdf_save_setting, settings, var=the_var, filename=out_file

    ; mlat.
    the_var = prefix+'mlat'
    val = get_var_data(the_var)
    settings = dictionary($
        'FIELDNAM', 'Spacecraft MLat', $
        'UNITS', 'deg', $
        'VAR_TYPE', 'data', $
        'DEPEND_0', time_var )
    cdf_save_var, the_var, value=val, filename=out_file
    cdf_save_setting, settings, var=the_var, filename=out_file

    ; lshell.
    the_var = prefix+'lshell'
    val = get_var_data(the_var)
    settings = dictionary($
        'FIELDNAM', 'Spacecraft L-shell', $
        'UNITS', '#', $
        'VAR_TYPE', 'data', $
        'DEPEND_0', time_var )
    cdf_save_var, the_var, value=val, filename=out_file
    cdf_save_setting, settings, var=the_var, filename=out_file


end

date = '2015-03-17'
probe = 'a'
file = join_path([homedir(),'test_ml_rbsp_orbit.cdf'])
file_delete, file, allow_nonexistent=1
ml_rbsp_load_orbit_var_gen_file, date, probe=probe, filename=file


end