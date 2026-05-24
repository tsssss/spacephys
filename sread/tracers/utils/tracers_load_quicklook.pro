;+
; Download all quicklook plots.
;-
function tracers_load_quicklook, input_time_range, probe=probe, id=data_type, $
    local_root=local_root, remote_root=remote_root, errmsg=errmsg

    compile_opt idl2
    errmsg = ''
    retval = !null
  
    if n_elements(probe) eq 0 then begin
        errmsg = 'No input probe ...'
        return, retval
    endif
    if not is_valid_probe(tracers_get_probes(),probe) then begin
        errmsg = 'Invalid probe ...'
        return, retval
    endif

;---Check inputs.
    if n_elements(local_root) eq 0 then local_root = tracers_get_local_root() 
    if n_elements(remote_root) eq 0 then remote_root = tracers_get_remote_root()

    probe_dir = strupcase('ts'+probe)
    quicklook_remote_path = join_path([remote_root,'SOC',probe_dir,'QL'])
    quicklook_local_path = join_path([local_root,'SOC',probe_dir,'QL'])
    if n_elements(data_type) eq 0 then data_type = 'roi'

    if n_elements(input_time_range) eq 0 then begin
        input_time_range = [time_double('2025'),systime(seconds=1)]
    endif
    if size(input_time_range[0],type=1) eq 7 then begin
        time_range = time_double(input_time_range)
    endif else begin
        time_range = input_time_range
    endelse

;---Obtain path information.
    year_range = float(sort_uniq(time_string(time_range,tformat='YYYY')))
    years = make_bins(year_range,1)
    year_strs = string(years,format='(I4)')
    local_paths = list()
    remote_paths = list()

    if data_type eq 'roi' then begin
        foreach year_str, year_strs do begin
            local_paths.add, join_path([quicklook_local_path,'roi',year_str])
            remote_paths.add, join_path([quicklook_remote_path,'roi',year_str])
        endforeach
    endif else if data_type eq 'coverage' then begin
        foreach year_str, year_strs do begin
            local_paths.add, join_path([quicklook_local_path,'coverage',year_str])
            remote_paths.add, join_path([quicklook_remote_path,'coverage',year_str])
        endforeach
    endif else begin
        foreach year_str, year_strs do begin
            local_paths.add, join_path([quicklook_local_path,year_str])
            remote_paths.add, join_path([quicklook_remote_path,year_str])
        endforeach
    endelse

;---Download files.
    username = tracers_get_remote_root('username')
    password = tracers_get_remote_root('password')
    npath = n_elements(local_paths)
    for ii=0,npath-1 do begin
        local_path = local_paths[ii]+'/'
        remote_path = remote_paths[ii]+'/'
        tmp = download_recursively(local_path, remote_path, $
            username=username, password=password)
    endfor

    return, 1

end

foreach probe, tracers_get_probes() do begin
    print, tracers_load_quicklook(probe=probe, id='roi')
endforeach

end