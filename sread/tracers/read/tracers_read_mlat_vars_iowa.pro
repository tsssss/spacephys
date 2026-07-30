
function tracers_read_mlat_vars_iowa, input_time_range, probe=probe, $
    errmsg=errmsg, update=update, get_name=get_name, suffix=suffix, _extra=extra
    compile_opt idl2

    errmsg = ''
    retval = !null

    if n_elements(suffix) eq 0 then suffix = '_iowa'
    r_var = tracers_read_orbit(get_name=1, probe=probe)
    var_info = lets_read_mlat_vars(get_name=1, orbit_var=r_var, probe=probe, suffix=suffix)
    if keyword_set(get_name) then return, var_info

    if keyword_set(update) then tmp = delete_var_from_file(var_info)
    time_range = time_double(input_time_range)
    if ~check_if_update_memory(var_info, time_range) then return, var_info


    r_var = tracers_read_orbit(time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval
    mlat_vars = lets_read_mlat_vars(orbit_var=r_var, errmsg=errmsg)

    return, mlat_vars

end


compile_opt idl2
time_range = ['2026-02-16/04:00','2026-02-16/04:05']
probe = '2'
mlat_vars = tracers_read_mlat_vars_iowa(time_range, probe=probe)
tplot, (mlat_vars.values()).toarray(), trange=time_range
end