;+
; Calculate elevation angle of a vector (B, E, etc).
; vec_var.
; coord=.
;-
function lets_calc_vec_elev, vec_var, coord=coord_need, $
    var_info=var_info, $
    update=update, get_name=get_name, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var

    short_name = get_var_setting(vec_var, 'short_name')
    if n_elements(var_info) eq 0 then begin
        var_info = strlowcase(get_prefix(vec_var)+short_name+'_elev')
    endif

    ; Check if update in memory.
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info

    the_vec = get_var_data(vec_var, times=times, settings=settings)
    coord = settings['coord']
    if n_elements(coord_need) eq 0 then coord_need = coord
    if coord ne coord_need then begin
        if settings.haskey('probe') then probe = settings['probe'] else probe = !null
        the_vec = cotran_pro(the_vec, times, coord_msg=[coord,coord_need], probe=probe)
    endif
    vec_elev = asin(the_vec[*,2]/snorm(the_vec))*constant('deg')

    store_data, var_info, times, vec_elev
    add_setting, var_info, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', short_name+'!D'+tex2str('theta'), $
        'unit', 'deg' )
    foreach key, ['requested_time_range','probe','mission'] do begin
        if settings.haskey(key) then options, var_info, key, settings[key]
    endforeach

    ; Save to file.
    if n_elements(data_file) ne 0 then begin
        is_success = save_var_to_file(var_info, file=data_file, time_var=time_var)
    endif
    
    return, var_info
end