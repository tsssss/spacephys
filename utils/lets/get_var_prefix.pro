;+
;-

function get_var_prefix, var

    prefix = get_var_setting(var, 'prefix', exist)
    if exist then return, prefix

    mission_probe = get_var_setting(var, 'mission_probe', exist)
    if exist then begin
        probe_info = resolve_probe(mission_probe)
        if n_elements(probe_info) ne 0 then return, probe_info['prefix']
    endif

    return, get_prefix(var)

end