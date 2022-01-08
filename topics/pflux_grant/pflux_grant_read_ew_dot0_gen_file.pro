;+
; Get Ew and Bw_ratio in MGSE.
;-

pro pflux_grant_read_ew_dot0_gen_file, time, probe=probe, filename=file, $
    local_root=local_root

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

    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','rbsp'])


;---Settings.
    secofday = constant('secofday')
    date = time[0]-(time[0] mod secofday)
    day_time_range = date+[0,secofday]
    time_range = day_time_range

    ndim = 3
    prefix = 'rbsp'+probe+'_'
    xyz = constant('xyz')
    fillval = !values.f_nan


;---Load data and get Ew and Bw_ratio.
    version = 'v01'
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_bfield_l2_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'bfield','l2','%Y']
    bfield_file = join_path([local_path,base_name])
    bfield_file = apply_time_to_pattern(bfield_file, date)
    b0_gsm_var = prefix+'b0_gsm'
    cdf_load_var, b0_gsm_var, filename=bfield_file

;    base_name = rbspx+'_efield_l1_%Y_%m%d_'+version+'.cdf'
;    local_path = [local_root,rbspx,'efield','l1','%Y']
;    efield_file = join_path([local_path,base_name])
;    efield_file = apply_time_to_pattern(efield_file, date)
;    cdf_load_var, prefix+'e_mgse', filename=efield_file
    pflux_grant_read_e_mgse, time_range, probe=probe

    b_gsm = get_var_data(b0_gsm_var, times=common_times)
    b_mgse = cotran(b_gsm, common_times, 'gsm2mgse', probe=probe)
    bw_ratio = b_mgse[*,0]/snorm(b_mgse)
    bw_ratio_var = prefix+'bw_ratio'
    store_data, bw_ratio_var, common_times, bw_ratio

    e_mgse = get_var_data(prefix+'e_mgse', limits=lim)
    edot0_mgse = e_mgse
    ew_dot0 = -(e_mgse[*,1]*b_mgse[*,1]+e_mgse[*,2]*b_mgse[*,2])/b_mgse[*,0]
    ew_dot0_var = prefix+'ew_dot0'
    store_data, ew_dot0_var, common_times, ew_dot0


;---Trim to day_time_range.
    times = common_times
    bw_ratio = float(get_var_data(bw_ratio_var))
    ew_dot0 = float(get_var_data(ew_dot0_var))


;---Save data to file.
    path = file_dirname(file)
    if file_test(path) eq 0 then file_mkdir, path
    cdf_id = (file_test(file))? cdf_open(file): cdf_create(file)

    time_var = 'ut_sec'
    time_var_settings = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix')
    cdf_save_var, time_var, value=times, filename=cdf_id
    cdf_save_setting, time_var_settings, filename=cdf_id, varname=time_var

    ; Bw ratio.
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'scalar', $
        'short_name', 'B!Dw!N/|B|', $
        'unit', '#' )
    cdf_save_var, bw_ratio_var, value=bw_ratio, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=bw_ratio_var

    ; Ew dot0.
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'scalar', $
        'short_name', 'E!Dw!N', $
        'unit', 'mV/m' )
    cdf_save_var, ew_dot0_var, value=ew_dot0, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=ew_dot0_var

    cdf_close, cdf_id

end



;probe = 'b'
;time = time_double('2012-12-04')    ; gap.
;time = time_double('2013-06-07')    ; Good case.
;
;file = join_path([homedir(),'test.cdf'])
;pflux_grant_read_ew_dot0_gen_file, time, probe=probe, filename=file
;stop

valid_range = time_double(['2012-09-30','2015-10-02'])
print, time_string(valid_range)
secofday = constant('secofday')
days = make_bins(valid_range+[0,-1], secofday, /inner)

local_root = join_path([default_local_root(),'sdata','rbsp'])
version = 'v01'

probes = ['b']
foreach probe, probes do begin
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_efield_l2_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'efield','l2','%Y']
    foreach day, days do begin
        print, time_string(day)
        file = join_path([local_path,base_name])
        file = apply_time_to_pattern(file, day)
        pflux_grant_read_ew_dot0_gen_file, day, probe=probe, filename=file
    endforeach
endforeach


end
