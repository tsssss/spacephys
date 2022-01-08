;+
; Preprocess B1 field to:
;   1. Rotate dB field to FAC.
;-

pro pflux_grant_read_b1_fac_gen_file, time, probe=probe, filename=file

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    if n_elements(file) eq 0 then begin
        errmsg = handle_error('No output file ...')
        return
    endif

    if n_elements(probe) eq 0 then begin
        errmsg = handle_error('No input probe ...')
        return
    endif

    if n_elements(time) eq 0 then begin
        errmsg = handle_error('No input time ...')
        return
    endif


;---Settings.
    secofday = constant('secofday')
    date = time[0]-(time[0] mod secofday)
    day_time_range = date+[0,1]*secofday
    time_range = day_time_range
    project = pflux_grant_load_project()
    common_time_step = project.common_time_step
    common_times = make_bins(time_range, common_time_step)
    common_times = common_times[0:-2]
    ncommon_time = n_elements(common_times)
    ndim = 3
    prefix = 'rbsp'+probe+'_'
    xyz = constant('xyz')
    fac = ['earthward','westward','outward']
    fillval = !values.f_nan


;---Load E field.
    pflux_grant_read_bfield, time_range, probe=probe
    pflux_grant_read_q_gsm2fac, time_range, probe=probe

    q_var = prefix+'q_gsm2fac'
    get_data, q_var, times, q_gsm2fac
    q_gsm2fac = qslerp(q_gsm2fac, times, common_times)
    store_data, q_var, common_times, q_gsm2fac

    ; Convert to FAC.
    gsm_var = prefix+'b1_gsm'
    add_setting, gsm_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'dB', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz )
    to_fac, gsm_var, to=prefix+'b1_fac', q_var=q_var


;---Save data to file.
    cdf_id = (file_test(file))? cdf_open(file): cdf_create(file)

    time_var = 'ut_sec'
    times = common_times
    time_var_settings = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix')
    cdf_save_var, time_var, value=times, filename=cdf_id
    cdf_save_setting, time_var_settings, filename=cdf_id, varname=time_var

    var = prefix+'b1_fac'
    data = float(get_var_data(var))
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'dB', $
        'unit', 'nT', $
        'coord', 'FAC', $
        'coord_labels', fac)
    cdf_save_var, var, value=data, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=var

    cdf_close, cdf_id

end


probe = 'b'
time = time_double('2013-06-07')
file = join_path([homedir(),'test.cdf'])
if file_test(file) eq 1 then file_delete, file
tic
pflux_grant_read_b1_fac_gen_file, time, probe=probe, filename=file
toc
end
