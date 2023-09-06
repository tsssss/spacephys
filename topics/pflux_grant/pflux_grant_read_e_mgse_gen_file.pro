;+
; Preprocess E field for the whole mission, to:
;   1. remove perigee field.
;-

pro pflux_grant_read_e_mgse_gen_file, time, probe=probe, filename=file;, pad_time=pad_time

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
    pad_time = 30.  ; for removing DC offset.
    time_range = day_time_range+[-1,1]*pad_time
;time_range = date+[-1,1]*secofday
    project = pflux_grant_load_project()
    common_time_step = project.common_time_step
    common_times = make_bins(time_range, common_time_step)
    ncommon_time = n_elements(common_times)
    ndim = 3
    prefix = 'rbsp'+probe+'_'
    xyz = constant('xyz')
    fillval = !values.f_nan
    min_ntime = 1800d*64    ; 0.5 hour of data.

;---Load E field.
    e_mgse_var = prefix+'e_mgse'
    rbsp_efw_read_e_mgse, time_range, probe=probe
    interp_time, e_mgse_var, common_times
    add_setting, e_mgse_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'mV/m', $
        'short_name', 'E', $
        'coord', 'MGSE', $
        'coord_labels', constant('xyz') )
    e_mgse = get_var_data(e_mgse_var)


;---Calculate E-E model.
    rbsp_read_e_model, time_range, probe=probe
    emod_var = prefix+'emod_mgse'
    interp_time, emod_var, common_times, /quadratic
    add_setting, emod_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'mV/m', $
        'short_name', 'Model E', $
        'coord', 'MGSE', $
        'coord_labels', xyz)
    emod_mgse = get_var_data(emod_var)

    de_mgse = e_mgse-emod_mgse
    de_mgse[*,0] = 0
    de_mgse_var = prefix+'de_mgse'
    store_data, de_mgse_var, common_times, de_mgse
    add_setting, de_mgse_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'mV/m', $
        'short_name', 'Raw dE', $
        'coord', 'MGSE', $
        'coord_labels', xyz)



;---Save data to file.
    index = where_pro(common_times, '[)', day_time_range)
    times = common_times[index]
    de_mgse = get_var_data(de_mgse_var)
    de_mgse = float(de_mgse[index,*])



    cdf_id = (file_test(file))? cdf_open(file): cdf_create(file)

    time_var = 'ut_sec'
    time_var_settings = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix')
    cdf_save_var, time_var, value=times, filename=cdf_id
    cdf_save_setting, time_var_settings, filename=cdf_id, varname=time_var


    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'E', $
        'unit', 'mV/m', $
        'coord', 'MGSE', $
        'coord_labels', xyz)
    cdf_save_var, e_mgse_var, value=de_mgse, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=e_mgse_var

    cdf_close, cdf_id

end

probe = 'b'
;time = time_double('2012-12-04')    ; gap.
;time = time_double('2012-10-02')    ; spikes.
time = time_double('2014-08-28')    ; spin tone.
time = time_double('2014-06-15')
;time = time_double('2014-06-14')

probe = 'b'
time = time_double('2014-06-17')
time = time_double('2015-05-15')    ; flags on.

file = join_path([homedir(),'test.cdf'])
if file_test(file) then file_delete, file
pflux_grant_read_e_mgse_gen_file, time, probe=probe, filename=file


end
