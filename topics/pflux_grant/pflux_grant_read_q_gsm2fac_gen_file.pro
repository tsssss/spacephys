;+
; Preprocess the quaternion to rotate data from GSM to FAC.
;-
pro pflux_grant_read_q_gsm2fac_gen_file, time, probe=probe, filename=file, local_root=local_root

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
    time_range = date+[0,secofday]
    prefix = 'rbsp'+probe+'_'
    xyz = constant('xyz')
    fillval = !values.f_nan

;---Read B0, R GSM.
    b0_gsm_var = prefix+'b0_gsm'
    pflux_grant_read_bfield, time_range, probe=probe
    rbsp_read_orbit, time_range, probe=probe
    r_var = prefix+'r_gse'
    r_gse = get_var_data(r_var, times=orbit_times)
    r_gsm = cotran(r_gse, orbit_times, 'gse2gsm')
    store_data, prefix+'r_gsm', orbit_times, r_gsm
    add_setting, prefix+'r_gsm', /smart, dictionary($
        'diaplay_type', 'vector', $
        'short_name', 'R', $
        'unit', 'Re', $
        'coord', 'GSM', $
        'coord_labels', xyz )

;---Define q_gsm2fac.
    define_fac, prefix+'b0_gsm', prefix+'r_gsm', time_var=prefix+'r_gsm'


;---Save data to file.
    cdf_id = (file_test(file))? cdf_open(file): cdf_create(file)

    ; Orbit time.
    time_var = 'ut_orbit'
    time_var_settings = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix')
    cdf_save_var, time_var, value=orbit_times, filename=cdf_id
    cdf_save_setting, time_var_settings, filename=cdf_id, varname=time_var

    var = prefix+'q_gsm2fac'
    data = get_var_data(var)
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'Q', $
        'unit', '#', $
        'coord', 'GSM2FAC', $
        'coord_labels', ['a','b','c','d'], $
        'colors', sgcolor(['red','green','blue','black']) )
    cdf_save_var, var, value=data, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=var

    cdf_close, cdf_id

end

probe = 'b'
time = time_double('2013-06-07')
file = join_path([homedir(),'test.cdf'])
if file_test(file) eq 1 then file_delete, file
tic
pflux_grant_read_q_gsm2fac_gen_file, time, probe=probe, filename=file
toc
end