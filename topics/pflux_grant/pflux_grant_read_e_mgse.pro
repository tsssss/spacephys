;+
; Read E MGSE.
;-

pro pflux_grant_read_e_mgse, time, probe=probe, id=data_type, $
    edot0=edot0, $
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
    data_type = 'e_mgse'

;---Init settings.
    type_dispatch = hash()
    valid_range = time_double(['2012-09-30','2015-10-02'])
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_efield_l1_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'efield','l1','%Y']
    ; B field.
    type_dispatch['e_mgse'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', rbspx+'_'+['e']+'_mgse', $
                'time_var_name', 'ut_sec', $
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
            pflux_grant_read_e_mgse_gen_file, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif

;---Read data from files and save to memory.
    read_files, time, files=files, request=request


    prefix = rbspx+'_'
    e_var = prefix+'e_mgse'

;---Edot0.
    if keyword_set(edot0) then begin
        pflux_grant_read_ew_dot0, time, probe=probe
        e_mgse = get_var_data(e_var, times=common_times)
        ew = get_var_data(prefix+'ew_dot0')
        e_mgse[*,0] = ew
        edot0_var = prefix+'edot0_mgse'
        store_data, edot0_var, common_times, e_mgse
        add_setting, edot0_var, /smart, {$
            display_type: 'vector', $
            unit: 'mV/m', $
            short_name: 'dot0 E', $
            coord: 'MGSE', $
            coord_labels: ['x','y','z'] }
    endif else begin
        if tnames(e_var) ne '' then begin
            add_setting, e_var, /smart, {$
                display_type: 'vector', $
                unit: 'mV/m', $
                short_name: 'E', $
                coord: 'MGSE', $
                coord_labels: ['x','y','z'] }
        endif
    endelse

end



probes = ['b']
secofday = constant('secofday')
foreach probe, probes do begin
    rbspx = 'rbsp'+probe
    valid_range = time_double(['2012-10-01','2015-10-01'])
    print, time_string(valid_range)
    days = make_bins(valid_range+[0,-1], secofday, /inner)
    foreach day, days do begin
        time_range = day+[0,secofday]
        pflux_grant_read_e_mgse, time_range, probe=probe
    endforeach
endforeach
end
