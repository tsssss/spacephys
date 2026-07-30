
function tracers_read_mlat_vars, input_time_range, probe=probe, $
    update=update, get_name=get_name, suffix=suffix, _extra=extra

    sources = ['iowa','cdaweb']
    foreach source, sources do begin
        func_name = 'tracers_read_mlat_vars_'+source
        retval = call_function(func_name, input_time_range, probe=probe, update=update, get_name=get_name, suffix=suffix, _extra=extra)
        if n_elements(retval) ne 0 then return, retval
    endforeach

    return, []

end