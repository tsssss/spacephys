
probes = ['a','b']
time_range = time_double(['2018-08-26','2018-09-01'])

foreach probe, probes do begin
    rbsp_read_e_model, time_range, probe=probe
    rbsp_efw_prepare_residue_removal, time_range, probe=probe

    prefix = 'rbsp'+probe+'_'
    emod_mgse = get_var_data(prefix+'emod_mgse', times=times)
    e_mgse = get_var_data(prefix+'e_mgse', at=times, limits=lim)
    de_mgse = e_mgse-emod_mgse

    ; Remove spin-axis data.
    de_mgse[*,0] = !values.f_nan

    ; Remove apogee data.
    dis = snorm(get_var_data(prefix+'r_mgse', at=times))
    index = where(dis ge 2)
    de_mgse[index,*] = !values.f_nan

    ; Remove other invalide data.
    rbsp_efw_read_boom_flag, time_range, probe=probe
    flags = total(get_var_data(prefix+'boom_flag', times=uts),2)
    index = where(flags ne 4, count)
    if count ne 0 then begin
        bad_time_ranges = uts[time_to_range(index,time_step=1)]
        nbad_time_range = n_elements(bad_time_ranges)*0.5
        for ii=0,nbad_time_range-1 do begin
            index = where_pro(times, '[]', bad_time_ranges[ii,*]+[-1,1]*600, count=count)
            if count eq 0 then continue
            de_mgse[index,*] = !values.f_nan
        endfor
    endif

    store_data, prefix+'de_mgse', times, de_mgse, limits=lim

    save_vars = ['dst',prefix+['de','e','emod','r']+'_mgse']
    data_file = join_path([homedir(),prefix+'perigee_e_mgse_zhenpeng_su_request.tplot'])
    tplot_save, save_vars, filename=data_file
endforeach


end
