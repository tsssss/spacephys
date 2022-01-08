;+
; Level 1 data includes bmod_gsm, b_gsm, e_mgse.
; All data are on uniform time at 16 Samples/sec.
;-

pro pflux_grant_read_level1_data, time, id=data_type, probe=probe, $
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
    base_name = rbspx+'_ebfield_l1_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'ebfield','l1','%Y']
    ; E field.
    type_dispatch['efield'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', rbspx+'_'+['e_uvw'], $
                'time_var_name', 'ut_field', $
                'time_var_type', 'unix')))
    ; B field.
    type_dispatch['bfield'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', rbspx+'_'+['bmod_gsm','b_gsm'], $
                'time_var_name', 'ut_field', $
                'time_var_type', 'unix')))
    ; E and B field.
    type_dispatch['ebfield'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', rbspx+'_'+['bmod_gsm','b_gsm','e_uvw'], $
                'time_var_name', 'ut_field', $
                'time_var_type', 'unix')))

    if keyword_set(print_data_type) then begin
        print, 'Suported data type: '
        ids = type_dispatch.keys()
        foreach id, ids do print, '  * '+id
        return
    endif


;---Dispatch patterns.
    if n_elements(data_type) eq 0 then begin
        errmsg = handle_error('No input data_type ...')
        return
    endif
    if not type_dispatch.haskey(data_type) then begin
        errmsg = handle_error('Do not support type '+data_type+' yet ...')
        return
    endif
    request = type_dispatch[data_type]

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)
    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            pflux_grant_gen_level1_data, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif

;---Read data from files and save to memory.
    read_files, time, files=files, request=request


    prefix = rbspx+'_'
    the_var = prefix+'b_gsm'
    if tnames(the_var) ne '' then begin
        add_setting, the_var, /smart, {$
            display_type: 'vector', $
            unit: 'nT', $
            short_name: 'B', $
            coord: 'GSM', $
            coord_labels: ['x','y','z'] }
    endif

    the_var = prefix+'bmod_gsm'
    if tnames(the_var) ne '' then begin
        add_setting, the_var, /smart, {$
            display_type: 'vector', $
            unit: 'nT', $
            short_name: 'T89 B', $
            coord: 'GSM', $
            coord_labels: ['x','y','z'] }
    endif

    the_var = prefix+'e_uvw'
    if tnames(the_var) ne '' then begin
        add_setting, the_var, /smart, {$
            display_type: 'vector', $
            unit: 'mV/m', $
            short_name: 'E', $
            coord: 'UVW', $
            coord_labels: ['u','v','w'] }
    endif

end


project = pflux_grant_load_project()
pflux_setting = project.pflux_calc_setting
probes = ['a']
secofday = constant('secofday')
foreach probe, probes do begin
    rbspx = 'rbsp'+probe
    valid_range = pflux_setting[rbspx].time_range+[-1,1]*secofday
    print, time_string(valid_range)
    days = make_bins(valid_range+[0,-1], secofday, /inner)
    foreach day, days do begin
        time_range = day+[0,secofday]
        pflux_grant_read_level1_data, time_range, probe=probe, id='ebfield'
    endforeach
endforeach
end
