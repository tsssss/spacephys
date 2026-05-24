function tracers_load_ephemeris, input_time_range, id=datatype, probe=probe, $
    print_datatype=print_datatype, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, remote_root=remote_root
    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','tracers'])
    if n_elements(remote_root) eq 0 then remote_root = tracers_get_remote_root()
    if n_elements(version) eq 0 then version = 'v[0-9.]+'

    if size(input_time_range[0],type=1) eq 7 then begin
        time_range = time_double(input_time_range)
    endif else begin
        time_range = input_time_range
    endelse

;---Init settings.
  

end
