
function lets_calc_edotb0, e_var=e_var, b_var=b_var, $
    var_info=var_info, $
    update=update, get_name=get_name, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var

    if n_elements(var_info) eq 0 then var_info = streplace(e_var, 'e', 'edot0')

    ; Check if update in memory.
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info


    e_vec = get_var_data(e_var, times=times, settings=e_setting)
    e_coord = e_setting['coord']
    probe = e_setting['probe']

    b_vec = get_var_data(b_var, at=times, settings=b_setting)
    b_coord = strlowcase(b_setting['coord'])
    if b_coord ne e_coord then begin
        b_vec = cotran_pro(b_vec, times, coord_msg=[b_coord,e_coord], probe=probe)
    endif

    sa_component = e_setting['spin_axis']
    coord_labels = e_setting['coord_labels']
    sa_index = where(coord_labels eq sa_component, complement=sp_index)
    e_vec[*,sa_index] = -total(e_vec[*,sp_index]*b_vec[*,sp_index],2)/b_vec[*,sa_index]
    store_data, var_info, times, e_vec, limits=e_setting.tostruct()

    ; B angle.
    b_angle = asin(b_vec[*,sa_index]/snorm(b_vec))*constant('deg')
    add_setting, var_info, dictionary('b_angle', b_angle)

    ; Save to file.
    if n_elements(data_file) ne 0 then begin
        is_success = save_var_to_file(var_info, file=data_file, time_var=time_var)
    endif
    
    return, var_info

end