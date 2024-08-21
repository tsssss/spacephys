;+
; Read B0 and B1, which is B-B_model plus smoohting over a specific window.
;-

function lets_decompose_bfield, var_info=var_info, $
    b0_window=b0_window, b_var=b_var, bmod_var=bmod_var, $
    update=update, get_name=get_name, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var, $
    coord=coord, resolution=resolution, _extra=ex

    errmsg = ''
    retval = !null

    ; Need b_var and bmod_var.
    if n_elements(b_var) eq 0 then begin
        errmsg = 'No b_var ...'
        return, retval
    endif
    
    if n_elements(bmod_var) eq 0 then begin
        errmsg = 'No bmod_var ...'
        return, retval
    endif

    ; Need mission_probe
    mission_probe = get_var_setting(b_var, 'mission_probe')
    probe_info = resolve_probe(mission_probe)
    prefix = probe_info['prefix']
    time_range = get_var_setting(b_var, 'requested_time_range')

    ; Get var_info.
    if n_elements(b0_window) eq 0 then b0_window = 20*60d   ; sec.
    if n_elements(coord) eq 0 then coord = 'gsm'
    if n_elements(suffix) eq 0 then suffix = ''
    if n_elements(var_info) eq 0 then begin
        var_info = dictionary($
            'b0', prefix+'b0_'+coord+suffix, $
            'b1', prefix+'b1_'+coord+suffix )
    endif
    if keyword_set(get_name) then return, var_info

    ; Check if update in memory.
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info

    ; Read var from routine.
    if n_elements(bmod_var) eq 0 then begin
        bmod_var = lets_read_geopack_bfield(time_range, probe=mission_probe, coord=coord, $
        external_model='t89', internal_model='igrf', $
        errmsg=errmsg, _extra=ex)
        if errmsg ne '' then return, retval
    endif
    b_vec = get_var_data(b_var, times=times)
    bmod_vec = get_var_data(bmod_var, at=times)
    b1_vec = b_vec-bmod_vec
    time_step = sdatarate(times)
    width = b0_window/time_step
    ndim = 3
    for ii=0,ndim-1 do begin
        b1_vec[*,ii] -= smooth(b1_vec[*,ii],width, nan=1, edge_mirror=1)
    endfor
    b0_vec = b_vec-b1_vec
    settings = get_var_setting(b_var)
    settings['short_name'] = 'B0'
    settings['b0_window'] = b0_window
    settings['external_model'] = external_model
    settings['internal_model'] = internal_model
    var_info['b0'] = save_data_to_memory(var_info['b0'], times, b0_vec, settings=settings)
    settings['short_name'] = 'B1'
    var_info['b1'] = save_data_to_memory(var_info['b1'], times, b1_vec, settings=settings)

    
    ; Save to file.
    if n_elements(data_file) ne 0 then begin
        is_success = save_var_to_file(out_vars, file=data_file, time_var=time_var)
    endif
    
    return, var_info


end


time_range = ['2015-03-17','2015-03-18']
probes = ['rbspa','rbspb']
data_file = join_path([homedir(),'test.cdf'])
;if file_test(data_file) eq 1 then file_delete, data_file
foreach probe, probes do begin
    prefix = probe+'_'
    r_var = lets_read_orbit(time_range, probe=probe)
    bmod_var = lets_read_geopack_bfield(orbit_var=r_var, external_model='t89', internal_model='igrf')
    b_var = lets_read_bfield(time_range, probe=probe, resolution='hires')
    b_vars = lets_decompose_bfield(b_var=b_var, bmod_var=bmod_var)
    stop
endforeach
end
