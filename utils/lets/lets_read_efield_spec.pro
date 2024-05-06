function lets_read_efield_spec, time_range, probe=mission_probe, $
    update=update, get_name=get_name, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var, $
    _extra=ex

    ; Need prefix, probe, routine_name, short_name
    if size(mission_probe,type=1) eq 7 then begin
        probe_info = resolve_probe(mission_probe)
    endif else begin
        probe_info = mission_probe
    endelse
    routine_name = probe_info['routine_name']
    probe = probe_info['probe']

    ; Get var_info.
    routine = routine_name+'_read_efield_spec'
    var_info = call_function(routine, time_range, probe=probe, $
        update=update, get_name=1, suffix=suffix, errmsg=errmsg, _extra=ex)
    
    ; Check if update in memory.
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info


    ; Read var from routine.
    var_info = call_function(routine, time_range, probe=probe, $
        update=update, get_name=0, suffix=suffix, errmsg=errmsg, _extra=ex)
    is_success = save_setting_to_memory(var_info, dictionary('mission_probe',probe_info['mission_probe']))
    if n_elements(time_var) ne 0 then is_success = interp_var_to_time(var_info, time_var=time_var)

    
    ; Save to file.
    if n_elements(data_file) ne 0 then begin
        is_success = save_var_to_file(var_info, file=data_file, time_var=time_var)
    endif
    
    return, var_info


end