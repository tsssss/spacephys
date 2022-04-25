;+
; Read RBSP orbit vars (r_coord, lshell, mlat, mlt)
;-

function ml_rbsp_read_orbit_var, input_time_range, probe=probe, coord=coord, errmsg=errmsg, get_name=get_name, vars=vars

    errmsg = ''
    retval = ''

    if n_elements(vars) eq 0 then vars = ['dis','lshell','mlt','mlat']

    pos_vars = list()
    pos_vars.add, ml_rbsp_read_pos(input_time_range, probe=probe, coord=coord, errmsg=errmsg)
    if errmsg ne '' then return, retval

    routines = 'ml_rbsp_read_'+vars
    foreach routine, routines do begin
        pos_vars.add, call_function(routine, input_time_range, probe=probe, errmsg=errmsg)
        if errmsg ne '' then return, retval
    endforeach
    return, pos_vars.toarray()

end

time_range = ['2013-06-07','2013-06-08']
probe = 'a'
vars = ml_rbsp_read_orbit_var(time_range, probe=probe)
end