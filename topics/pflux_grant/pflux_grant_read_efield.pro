;+
; Read E or E_dot0 in MGSE.
;-

pro pflux_grant_read_efield, time, probe=probe, id=data_type, $
    print_data_type=print_data_type, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, remote_root=remote_root, $
    min_bw_ratio = min_bw_ratio, coord=coord

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 0
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','rbsp'])
    if n_elements(version) eq 0 then version = 'v01'
    if n_elements(data_type) eq 0 then data_type = 'e'
    if n_elements(coord) eq 0 then coord = 'mgse'
    fillval = !values.f_nan


    prefix = 'rbsp'+probe+'_'
    e_mgse_var = prefix+'e_mgse'
    pflux_grant_read_e_mgse, time, probe=probe
    if data_type eq 'edot0' then begin
        edot0_mgse_var = prefix+'edot0_mgse'
        pflux_grant_read_ew_dot0, time, probe=probe, min_bw_ratio=min_bw_ratio
        ew = get_var_data(prefix+'ew_dot0')
        get_data, e_mgse_var, times, e_mgse
        e_mgse[*,0] = ew
        index = where(finite(ew,/nan), count)
        if count ne 0 then e_mgse[index,*] = fillval
        store_data, edot0_mgse_var, times, e_mgse
    endif
    
    if coord ne 'mgse' then begin
        e_var = prefix+data_type
        get_data, e_var+'_mgse', times, e_mgse
        e_coord = cotran(e_mgse, times, 'mgse2'+coord, probe=probe)
        store_data, e_var+'_'+coord, times, e_coord
    endif
   
    
    the_var = prefix+data_type+'_'+coord
    short_name = (data_type eq 'e')? 'E': 'E!S!UDot0!R!N'
    add_setting, the_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', short_name, $
        'unit', 'mV/m', $
        'coord', strupcase(coord), $
        'coord_labels', constant('xyz') )
        
        
    rbsp_efw_read_boom_flag, time, probe=probe
    flags = total(get_var_data(prefix+'boom_flag', times=uts),2) lt 4
    index = where(flags eq 1, count)
    if count ne 0 then begin
        nan_times = uts[time_to_range(index,time_step=1)]
        nnan_time = n_elements(nan_times)*0.5
        get_data, the_var, times, e_mgse
        for ii=0,nnan_time-1 do begin
            index = where_pro(times, '[]', nan_times[ii,*], count=count)
            if count eq 0 then continue
            e_mgse[index,*] = fillval
        endfor
        store_data, the_var, times, e_mgse
    endif
    
end

time_range = time_double(['2013-06-07','2013-06-08'])
probe = 'a'
pflux_grant_read_efield, time_range, probe=probe, id='edot0'
end