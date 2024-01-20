
function event_study_read_weygand, event_info, time_var=time_var, get_name=get_name, update=update

    types = ['hor','ver']
    j_vars = 'thg_j_'+types
    if keyword_set(get_name) then return, j_vars
    time_range = event_info['time_range']
    load_data = 0
    foreach var, j_vars do begin
        if check_if_update(var, time_range) then begin
            load_data = 1
            break
        endif
    endforeach

    if keyword_set(update) then load_data = 1
    if load_data eq 0 then return, j_vars

    data_file = event_info['data_file']
    if keyword_set(update) then foreach var, j_vars do cdf_del_var, var, filename=data_file
    if n_elements(time_var) eq 0 then time_var = 'weygand_time'

    foreach type, types do begin
        var = 'thg_j_'+type
        if ~cdf_has_var(var, filename=data_file) then begin
            var = themis_read_weygand_j(time_range, id='j_'+type)

            if ~cdf_has_var(time_var, filename=data_file) then begin
                times = get_var_time(var)
                cdf_save_var, time_var, value=times, filename=data_file
                time_step = 3d
                settings = dictionary($
                    'time_step', time_step, $
                    'data_time_range', time_range )
                cdf_save_setting, settings, filename=data_file, varname=time_var
            endif

            data = get_var_data(var, limits=limits)
            cdf_save_var, var, value=data, filename=data_file
            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=var
        endif

        if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    endforeach

    return, j_vars


end