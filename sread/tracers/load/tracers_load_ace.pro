
function tracers_load_ace, input_time_range, id=datatype, probe=probe, $
    print_datatype=print_datatype, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, remote_root=remote_root
    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'tracers'])
    if n_elements(remote_root) eq 0 then remote_root = tracers_get_remote_root()
    if n_elements(version) eq 0 then version = 'v[0-9.]+'
    if n_elements(datatype) eq 0 then datatype = 'iowa%l3'


    if size(input_time_range[0],type=1) eq 7 then begin
        time_range = time_double(input_time_range)
    endif else begin
        time_range = input_time_range
    endelse

;---Init settings.
    type_dispatch = hash()
    instr_str = 'ace'
    valid_range = tracers_get_valid_range(probe=probe, id=instr_str)
    probe_str = 'ts'+probe
    prefix = 'ts'+probe+'_'

    ; L3.
    key = 'iowa%l3'
    if datatype eq key then remote_root = 'https://tracers-portal.physics.uiowa.edu'
    remote_path = join_path([remote_root,'L3',strupcase(probe_str),'%Y','%%m','%d'])
    local_path = join_path([local_root,'L3',strupcase(probe_str),'%Y','%%m','%d'])
    base = prefix+'l3_ace_pitch-angle-dist_%Y%m%d_'+version+'.cdf'
    type_dispatch[key] = dictionary($
        'pattern', dictionary($
            'remote_file', join_path([remote_path,base]), $
            'remote_index_file', join_path([remote_path,'']), $
            'local_file', join_path([local_path,base]), $
            'local_index_file', join_path([local_path,default_index_file(/sync)])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base) )


    ; L2.
    key = 'cdaweb%l2'
    if datatype eq key then remote_root = 'https://cdaweb.gsfc.nasa.gov/pub/data/tracers'
    remote_path = join_path([remote_root,'tracers'+probe,'ace_cusp-electrons','l2','def_diff-en-flux','%Y'])
    local_path = join_path([local_root,'cdaweb','tracers'+probe,'ace_cusp-electrons','l2','def_diff-en-flux','%Y'])
    base = prefix+'l2_ace_def_%Y%m%d_'+version+'.cdf'
    type_dispatch[key] = dictionary($
        'pattern', dictionary($
            'remote_file', join_path([remote_path,base]), $
            'remote_index_file', join_path([remote_path,'']), $
            'local_file', join_path([local_path,base]), $
            'local_index_file', join_path([local_path,default_index_file(/sync)])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base) )


    ; L1.
    type_strs = ['epd_x2b2','epd_x2b3']
    foreach type_str, type_strs do begin
        foreach level_str, ['l1a','l1b'] do begin
            key = level_str+'%'+type_str
            remote_path = join_path([remote_root,'SOC',strupcase(probe_str),strupcase(level_str),strupcase(instr_str),'%Y','%%m','%d'])
            local_path = join_path([local_root,'SOC',strupcase(probe_str),strupcase(level_str),strupcase(instr_str),'%Y','%%m','%d'])
            base = probe_str+'_'+level_str+'_'+instr_str+'_'+type_str+'_%Y%m%d_'+version+'.cdf'
            type_dispatch[key] = dictionary($
                'username', tracers_get_remote_root('username'), $
                'password', tracers_get_remote_root('password'), $
                'pattern', dictionary($
                    'remote_file', join_path([remote_path,base]), $
                    'remote_index_file', join_path([remote_path,'']), $
                    'local_file', join_path([local_path,base]), $
                    'local_index_file', join_path([local_path,default_index_file(/sync)])), $
                'valid_range', time_double(valid_range), $
                'cadence', 'day', $
                'extension', fgetext(base) )
        endforeach
    endforeach


;---Dispatch patterns.
    if n_elements(datatype) eq 0 then begin
        errmsg = handle_error('No input datatype ...')
        return, ''
    endif
    if not type_dispatch.haskey(datatype) then begin
        errmsg = handle_error('Do not support type '+datatype+' yet ...')
        return, ''
    endif
    request = type_dispatch[datatype]

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time_range)
    
    if n_elements(files) eq 0 then return, '' else return, files
    
end

compile_opt idl2

; url = 'https://tracers-portal.physics.uiowa.edu/teams/flight/SOC/TS1/L1B/ACE/2025/11/22/'
url = 'https://tracers-portal.physics.uiowa.edu/teams/flight/SOC/TS2/L1B/ACE/2026/02/16/ts2_l1b_ace_epd_x2b2_20260216_v0.7.0.cdf'
; local_file = join_path([homedir(),'Downloads','tmp','tracers_url_test.html'])
; username = tracers_get_remote_root('username')
; password = tracers_get_remote_root('password')
; download_file, local_file, url, username=username, password=password
; stop

time_range = ['2025-11-22','2025-11-23']
probe = '2'
files = tracers_load_ace(time_range, probe=probe, id='iowa%l3')
print, files
stop
files = tracers_load_ace(time_range, probe=probe, id='l1b%ipd_x292')
print, files
files = tracers_load_ace(time_range, probe=probe, id='cdaweb%l2')
print, files
end