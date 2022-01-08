;+
; Level 2 data includes:
;   b0_gsm, b1_gsm, which are the background and perturbation magnetic field.
;   e_gsm, edot0_gsm, which are the electric field with spin=0 and e_dot_b=0 correction.
;-

pro pflux_grant_read_level2_bfield, time, probe=probe, $
    print_data_type=print_data_type, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, remote_root=remote_root


    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 0
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','rbsp'])
    if n_elements(version) eq 0 then version = 'v01'

;---Init settings.
    type_dispatch = hash()
    valid_range = rbsp_info('spice_data_range', probe=probe)
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_bfield_l2_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'ebfield','l2','bfield','%Y']

    request = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', rbspx+'_'+['b0_gsm','b1_gsm'], $
                'time_var_name', 'ut_field', $
                'time_var_type', 'unix')))

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)
    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            pflux_grant_gen_level2_bfield, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif

;---Read data from files and save to memory.
    read_files, time, files=files, request=request
    
    vars = rbspx+'_'+['b0_gsm','b1_gsm']
    foreach var, vars do begin
        add_setting, var, /smart, dictionary($
            'display_type', 'vector', $
            'unit', 'nT', $
            'short_name', 'B', $
            'coord', 'GSM', $
            'coord_labels', constant('xyz') )
    endforeach

end


;probe = 'b'
;secofday = constant('secofday')
;day = time_double('2015-02-22')
;time_range = day+[0,secofday]
;pflux_grant_read_level2_bfield, time_range, probe=probe
;stop


project = pflux_grant_load_project()
pflux_setting = project.pflux_calc_setting
probes = ['a']
secofday = constant('secofday')
foreach probe, probes do begin
    rbspx = 'rbsp'+probe
    valid_range = pflux_setting[rbspx].time_range
    ;valid_range = valid_range<time_double('2015-10-01')
    days = make_bins(valid_range+[0,-1], secofday, /inner)
    foreach day, days do begin
        time_range = day+[0,secofday]
        pflux_grant_read_level2_bfield, time_range, probe=probe
    endforeach
endforeach
end