;+
; Add data of var1 and var2 and save the results to var3
;-

function lets_add_vars, var1, var2, save_to=var3, base_on=time_var

    if n_elements(time_var) eq 0 then begin
        dr1 = sdatarate(get_var_time(var1))
        dr2 = sdatarate(get_var_time(var2))
        if dr1 le dr2 then time_var = var1 else time_var = var2
    endif
    times = get_var_time(time_var)
    data3 = get_var_data(var1, at=times, limts=lim)+get_var_data(var2, at=times)
    return, save_data_to_memory(var3, times, data3, limits=lim)

end