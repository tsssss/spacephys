;+
;-

pro rbsp_read_ebfield, time, id=datatype, probe=probe, $
    print_datatype=print_datatype, errmsg=errmsg, $
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
    valid_range = rbsp_info('emfisis_l3_data_range', probe=probe)
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_spice_products_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'spice_product','%Y']
    ; Orbit variables.
    type_dispatch['orbit'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', ['pos_gsm'], $
                'time_var_name', 'ut_pos', $
                'time_var_type', 'unix')))
