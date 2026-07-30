;+
; 
;-

function tracers_read_ele_en_spec_iowa, input_time_range, probe=probe, $
    update=update, get_name=get_name, suffix=suffix, _extra=extra
    compile_opt idl2

    errmsg = ''
    retval = !null
    prefix = 'ts'+probe+'_'

    if n_elements(suffix) eq 0 then suffix = '_iowa'
    var_info = prefix+'e_en_spec'+suffix
    if keyword_set(get_name) then return, var_info
    if keyword_set(update) then del_data, var_info
    time_range = time_double(input_time_range)
    if ~check_if_update(var_info, time_range) then return, var_info

    files = tracers_load_ace(time_range, probe=probe, errmsg=errmsg, id='l1b%epd_x2b3')
    if errmsg ne '' then return, retval

    var_list = list()
    in_vars = prefix+'l1b_ace_flux'
    time_var = 'Epoch'
    var_list.add, dictionary($
        'in_vars', in_vars, $
        'time_var_name', time_var, $
        'time_var_type', 'tt2000' )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    fluxs = var_get_data(in_vars, times=times)
    energy_var = prefix+'l1b_ace_energy'
    energy_bins = cdf_read_var(energy_var, filename=files[0])
    unit = '1/cm!U2!N-s-sr'
    en_specs = total(fluxs, 2, nan=1)
    var_info = var_store(var_info, en_specs, times, energy_bins)
    add_setting, var_info, smart=1, dictionary($
        'display_type', 'spec', $
        'unit', unit, $
        'ylog', 1, $
        'ytitle', 'Energy!C(eV)', $
        'zlog', 1, $
        'short_name', '' )
    return, var_info

end

compile_opt idl2
time_range = ['2025-11-22','2025-11-23']
time_range = ['2025-12-22','2025-12-23']
time_range = ['2026-02-16/04:00','2026-02-16/04:05']
probe = '2'
e_var = tracers_read_ele_en_spec_iowa(time_range, probe=probe)
i_var = tracers_read_ion_en_spec_iowa(time_range, probe=probe)
plot_vars = [e_var, i_var]
tplot, plot_vars, trange=time_range
end
