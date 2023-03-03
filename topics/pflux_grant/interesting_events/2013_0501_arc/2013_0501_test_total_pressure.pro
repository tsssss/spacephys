
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()

    prefix = event_info['prefix']
    probe = event_info['probe']
    time_range = event_info['time_range']
    snapshot_time = event_info['snapshot_time']

    dens = get_var_data(prefix+'e_density', times=times, limits=e_lim)*1
    p_dens = get_var_data(prefix+'p_density', at=times, limits=p_lim)
    o_dens = get_var_data(prefix+'o_density', at=times, limits=o_lim)
    store_data, prefix+'p_full_density', times, dens*p_dens/(p_dens+o_dens), limits=p_lim
    store_data, prefix+'o_full_density', times, dens*o_dens/(p_dens+o_dens), limits=o_lim
    store_data, prefix+'e_full_density', times, dens, limits=e_lim
    
    coefs = make_bins([1,3],0.1)
    targets = []
    
    foreach coef, coefs do begin

        tmp = get_var_data(prefix+'p_density', times=times, in=time_range)
        ntime = n_elements(times)
        foreach species, ['e','p','o'] do begin
            dens = get_var_data(prefix+species+'_density', at=times)
            tavg = get_var_data(prefix+species+'_temp', at=times)
            p_thermal = dens*tavg*1e6*1.6e-19*1e9 ; nPa.
            
            full_dens = get_var_data(prefix+species+'_full_density')
            tavg = (species eq 'e')? 200: 30
            dp = (full_dens-dens)*tavg*1e6*1.6e-19*1e9 ; nPa.
            p_thermal += dp
    
            p_thermal *= coef
    
            store_data, prefix+species+'_p_thermal', times, p_thermal, limits={$
                ytitle: '(nPa)', labels:species+' P thermal'}
        endforeach
    
        b_gsm = get_var_data(prefix+'b_gsm', at=times)
        bmag = snorm(b_gsm)
        p_mag = (bmag*1e-9)^2*0.5/(4*!dpi*1e-7)*1e9 ; nPa.
        store_data, prefix+'p_mag', times, p_mag, limits={$
            ytitle: '(nPa)', labels:'P mag'}
    
    
        p_total = fltarr(ntime)
        foreach var, prefix+['p_mag',['e','p','o']+'_p_thermal'] do begin
    ;        interp_time, var, to=var0
            p_total += get_var_data(var)
        endforeach
        store_data, prefix+'p_total', times, p_total, limits={ytitle:'(nPa)', labels:'P total'}
    
        target = stddev(p_total,nan=1)
        targets = [targets,target]
    endforeach
    

    tmp = min(targets, index)
    the_coef = coefs[index]
    the_coef = 1.5
    
    tmp = get_var_data(prefix+'p_density', times=times, in=time_range)
    ntime = n_elements(times)
    foreach species, ['e','p','o'] do begin
        dens = get_var_data(prefix+species+'_density', at=times)
        tavg = get_var_data(prefix+species+'_temp', at=times)
        p_thermal = dens*tavg*1e6*1.6e-19*1e9 ; nPa.

        full_dens = get_var_data(prefix+species+'_full_density')
        tavg = (species eq 'e')? 200: 30
        dp = (full_dens-dens)*tavg*1e6*1.6e-19*1e9 ; nPa.
        p_thermal += dp

        p_thermal *= the_coef

        store_data, prefix+species+'_p_thermal', times, p_thermal, limits={$
            ytitle: '(nPa)', labels:species+' P thermal'}
    endforeach

    b_gsm = get_var_data(prefix+'b_gsm', at=times)
    bmag = snorm(b_gsm)
    p_mag = (bmag*1e-9)^2*0.5/(4*!dpi*1e-7)*1e9 ; nPa.
    store_data, prefix+'p_mag', times, p_mag, limits={$
        ytitle: '(nPa)', labels:'P mag'}


    p_total = fltarr(ntime)
    foreach var, prefix+['p_mag',['e','p','o']+'_p_thermal'] do begin
        ;        interp_time, var, to=var0
        p_total += get_var_data(var)
    endforeach
    store_data, prefix+'p_total', times, p_total, limits={ytitle:'(nPa)', labels:'P total'}
    
    
    
    
;---Plot the result.
    var = prefix+'combo_pressure'
    the_vars = prefix+['p_mag',['e','p','o']+'_p_thermal','p_total']
    var = stplot_merge(the_vars, output=var, $
        labels=['P!DB','P '+['e-','H+','O+'], 'P Total'], $
        colors=sgcolor(['red','blue','green','purple','black']))
    options, var, 'ytitle', '(nPa)'
    options, var, 'yrange', [-1,20]
    options, var, 'ytickv', [0,10,20]
    options, var, 'yticks', 2
    options, var, 'yminor', 5
    options, var, 'ystyle', 1
    options, var, 'labflag', -1
    
    
    get_data, var, times, p_combo
    p_total = p_combo[*,-1]
    b_lobe = sqrt(p_total*((4*!dpi*1e-7)/1e9)*2)*1e9
    store_data, prefix+'b_lobe', times, b_lobe
    tplot, prefix+'combo_pressure'
    stop

end