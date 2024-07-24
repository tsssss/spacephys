;+
; Read MLat, MLon, MLT, L-Shell, |R| bundle from orbit_var.
;-

function lets_read_mlat_vars, var_info=var_info, orbit_var=orbit_var, $
    update=update, get_name=get_name, prefix=prefix, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var, $
    _extra=ex

    errmsg = ''
    retval = !null

    ; Need prefix
    if n_elements(orbit_var) eq 0 then begin
        errmsg = 'No input orbit_var ...'
        return, retval
    endif

    if n_elements(prefix) eq 0 then begin
        prefix = get_prefix(orbit_var)
        mission_probe = get_var_setting(orbit_var, 'mission_probe')
        probe_info = resolve_probe(mission_probe)
    endif
    time_range = get_var_setting(orbit_var, 'requested_time_range')

    ; Get var_info.
    if n_elements(suffix) eq 0 then suffix = ''
    if n_elements(var_info) eq 0 then var_info = dictionary($
        'dis', prefix+'dis'+suffix, $
        'lshell', prefix+'lshell'+suffix, $
        'mlat', prefix+'mlat'+suffix, $
        'mlt', prefix+'mlt'+suffix, $
        'mlon', prefix+'mlon'+suffix )
    if keyword_set(get_name) then return, var_info

    ; Check if update in memory.
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info

    ; Read from routine.
    r_vec = get_var_data(orbit_var, times=times, settings=settings)
    r_coord = strlowcase(settings['coord'])
    default_coord = 'mag'
    if r_coord ne default_coord then r_vec = cotran_pro(r_vec, times, coord_msg=[r_coord,default_coord])
    scalar_settings = inherit_setting(settings)
    ; MLat.
    scalar_settings['unit'] = 'deg'
    scalar_settings['short_name'] = 'MLat'
    mlats = calc_mlat(r_vec)
    var_info['mlat'] = save_data_to_memory(var_info['mlat'], times, mlats, settings=scalar_settings)

    scalar_settings['unit'] = 'deg'
    scalar_settings['short_name'] = 'MLon'
    mlons = calc_mlon(r_vec)
    var_info['mlon'] = save_data_to_memory(var_info['mlon'], times, mlons, settings=scalar_settings)

    scalar_settings['unit'] = 'h'
    scalar_settings['short_name'] = 'MLT'
    mlts = calc_mlt(times, mlons)
    var_info['mlt'] = save_data_to_memory(var_info['mlt'], times, mlts, settings=scalar_settings)

    scalar_settings['unit'] = 'Re'
    scalar_settings['short_name'] = '|R|'
    diss = snorm(r_vec)
    if var_info.haskey('dis') then begin
        var_info['dis'] = save_data_to_memory(var_info['dis'], times, diss, settings=scalar_settings)
    endif

    scalar_settings['unit'] = '#'
    scalar_settings['short_name'] = 'L-Shell'
    lshells = calc_lshell(mlats, diss)
    if var_info.haskey('lshell') then begin
        var_info['lshell'] = save_data_to_memory(var_info['lshell'], times, lshells, settings=scalar_settings)
    endif

    ; Convert to the wanted coord.
    if n_elements(time_var) ne 0 then is_success = interp_var_to_time(var_info, time_var=time_var)
    
    ; Save to file.
    if n_elements(data_file) ne 0 then begin
        is_success = save_var_to_file(var_info, file=data_file, time_var=time_var)
    endif
    
    return, var_info

end
