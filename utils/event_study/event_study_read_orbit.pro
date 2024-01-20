function event_study_read_orbit, event_info, time_var=time_var, coord=coord, update=update

    mission = event_info['mission']
    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    time_range = event_info['orbit_time_range']
    if n_elements(coord) eq 0 then coord = 'gsm'
    routine_name = mission+'_read_orbit'
    var = call_function(routine_name, time_range, probe=probe, coord=coord, get_name=1)
    if keyword_set(get_name) then return, var
    if keyword_set(update) then del_data, var
    if ~check_if_update(var, time_range) then return, var


    if keyword_set(update) then cdf_del_var, var, filename=data_file
    if ~cdf_has_var(var, filename=data_file) then begin
        if n_elements(time_var) eq 0 then time_var = event_study_get_orbit_time(get_name=1)
        times = event_study_get_orbit_time(event_info, time_var=time_var)
        data_tr = time_range+[-1,1]*60
        var = call_function(routine_name, data_tr, probe=probe, coord=coord)
        interp_time, var, times
        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    
    return, var

end