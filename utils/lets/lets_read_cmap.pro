function lets_read_cmap, var_info=var_info, $
    b_var=b_var, fb_var=fb_var, $
    update=update, get_name=get_name, prefix=prefix, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var, $
    _extra=ex

    errmsg = ''
    retval = !null

    ; Need prefix
    if n_elements(b_var) eq 0 then begin
        errmsg = 'No input b_var ...'
        return, retval
    endif
    if n_elements(fb_var) eq 0 then begin
        errmsg = 'No input b_var ...'
        return, retval
    endif

    if n_elements(prefix) eq 0 then begin
        mission_probe = get_var_setting(b_var, 'mission_probe')
        probe_info = resolve_probe(mission_probe)
        prefix = probe_info['prefix']
    endif

    ; Get var_info.
    if n_elements(suffix) eq 0 then suffix = ''
    if n_elements(var_info) eq 0 then var_info = prefix+'cmap'
    if keyword_set(get_name) then return, var_info

    ; Check if update in memory.
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info

    ; Read from routine.
    bmag = snorm(get_var_data(b_var, times=times, settings=b_settings))
    bfmag = snorm(get_var_data(fb_var, at=times, settings=fb_settings))
    cmap = bfmag/bmag
    settings = inherit_setting(fb_settings)
    settings['short_name'] = 'cmap'
    settings['unit'] = '#'
    foreach key, ['external_model','internal_model'] do begin
        if ~fb_settings.haskey(key) then continue
        settings[key] = fb_settings[key]
    endforeach
    return, save_data_to_memory(var_info, times, cmap, settings=settings)


end