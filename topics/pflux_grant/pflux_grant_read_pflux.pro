;+
; Pflux MGSE.
;-

pro pflux_grant_read_pflux, time, probe=probe, id=data_type, $
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
    data_type = 'pflux_fac'


;---Init settings.
    type_dispatch = hash()
    valid_range = time_double(['2012-10-01','2015-10-01'])
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_flux_l1_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'pflux','l1','%Y']
    ; B field.
    type_dispatch['pflux_fac'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', rbspx+'_'+['pf','pfdot0']+'_fac_norm', $
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
            pflux_grant_read_pflux_gen_file, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif

;---Read data from files and save to memory.
    read_files, time, files=files, request=request


    prefix = rbspx+'_'
    fac = ['earthward','westward','outward']
    foreach the_var, prefix+['pf','pfdot0']+'_fac_norm' do begin
        if tnames(the_var) ne '' then begin
            add_setting, the_var, /smart, {$
                display_type: 'vector', $
                unit: 'mW/m!U2!N', $
                short_name: 'S', $
                coord: 'FAC', $
                coord_labels: fac }
        endif
    endforeach

end

probes = ['a']
secofday = constant('secofday')
foreach probe, probes do begin
    rbspx = 'rbsp'+probe
    valid_range = time_double(['2013-10-16','2013-12-01'])
    ;valid_range = time_double(['2014-03-09','2014-10-01'])
    ;valid_range = time_double(['2015-09-23','2015-10-01'])
    print, time_string(valid_range)
    days = make_bins(valid_range+[0,-1], secofday, /inner)
    stop
    foreach day, days do begin
        time_range = day+[0,secofday]
        pflux_grant_read_pflux, time_range, probe=probe
    endforeach
endforeach
end
