;+
; Read Ew dot0.
; Ew is nan when B ratio < 0.2.
;-

pro pflux_grant_read_ew_dot0, time, probe=probe, id=data_type, $
    print_data_type=print_data_type, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, remote_root=remote_root, $
    min_bw_ratio=min_bw_ratio

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 0
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','rbsp'])
    if n_elements(version) eq 0 then version = 'v01'
    data_type = 'ew_dot0'
    if n_elements(min_bw_ratio) eq 0 then min_bw_ratio = 0.2

;---Init settings.
    type_dispatch = hash()
    valid_range = time_double(['2012-09-30','2015-10-02'])
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_efield_l2_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'efield','l2','%Y']
    ; Ew.
    type_dispatch['ew_dot0'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', rbspx+'_'+['ew_dot0','bw_ratio'], $
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
            pflux_grant_read_ew_dot0_gen_file, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif

;---Read data from files and save to memory.
    read_files, time, files=files, request=request


    prefix = rbspx+'_'
    the_var = prefix+'ew_dot0'
    if tnames(the_var) ne '' then begin
        add_setting, the_var, /smart, dictionary($
            'display_type', 'scalar', $
            'unit', 'mV/m', $
            'short_name', 'E!Dw!N' )
        get_data, prefix+'bw_ratio', times, bw_ratio
        get_data, the_var, times, ew_dot0
        index = where(abs(bw_ratio) le min_bw_ratio, count)
        if count ne 0 then begin
            ew_dot0[index] = !values.f_nan
            store_data, the_var, times, ew_dot0
        endif
    endif

end



probes = ['a','b']
secofday = constant('secofday')
foreach probe, probes do begin
    rbspx = 'rbsp'+probe
    valid_range = time_double(['2012-09-30','2015-10-02'])
    print, time_string(valid_range)
    days = make_bins(valid_range+[0,-1], secofday, /inner)
    foreach day, days do begin
        time_range = day+[0,secofday]
        pflux_grant_read_ew_dot0, time_range, probe=probe
    endforeach
endforeach
end
