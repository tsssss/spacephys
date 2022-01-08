;+
; Test pressure of e, H+, and O+. Also check B pressure.
;-


;---Load data and settings.
    event_info = _2015_0218_02_load_data()
    test = 0

    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']
    
    ; Calc full density.
    var0 = prefix+'e_density'
    get_data, var0, times
    ntime = n_elements(times)
    
    dens = get_var_data(prefix+'efw_density', at=times, limits=lim)
    store_data, prefix+'e_full_density', times, dens, limits=lim
    p_dens = get_var_data(prefix+'p_density', at=times, limits=p_lim)
    o_dens = get_var_data(prefix+'o_density', at=times, limits=o_lim)
    store_data, prefix+'p_full_density', times, dens*p_dens/(p_dens+o_dens), limits=p_lim
    store_data, prefix+'o_full_density', times, dens*o_dens/(p_dens+o_dens), limits=o_lim

    foreach species, ['e','p','o'] do begin
        rbsp_read_hope_moments, time_range, probe=probe, species=species
        dens = get_var_data(prefix+species+'_density', at=times)
        tavg = get_var_data(prefix+species+'_t_avg', at=times)
        p_thermal = dens*tavg*1e6*1.6e-19*1e9 ; nPa.
        
        full_dens = get_var_data(prefix+species+'_full_density')
        tavg = (species eq 'e')? 200: 50
        dp = (full_dens-dens)*tavg*1e6*1.6e-19*1e9 ; nPa.
        print, dp
        p_thermal += dp
        
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
    

        
    vars = prefix+['p_mag',['e','p','o']+'_p_thermal','p_total']
    var = prefix+'pressure'
    stplot_merge, vars, newname=var, $
        labels=['P!DB',['e-','H+','O+']+' P!Dth', 'Total P'], $
        colors=sgcolor(['red','blue','green','purple','black'])
    options, var, 'ytitle', '(nPa)'
    
    
    ; Beta.
    p_mag = get_var_data(prefix+'p_mag', at=times)
    pressure = get_var_data(prefix+'pressure')
    ndim = n_elements(pressure[0,*])
    beta = fltarr(ntime,ndim)
    for ii=0,ndim-1 do begin
        beta[*,ii] = pressure[*,ii]/p_mag
        if ii eq ndim-1 then beta[*,ii] -= 1
    endfor
    var = prefix+'beta'
    str_beta = tex2str('beta')
    store_data, var, times, beta[*,1:*], limits={ytitle:'(#)', constant:1}
    options, var, 'labels', str_beta+' '+['e-','H+','O+','Total']
    options, var, 'colors', sgcolor(['blue','green','purple','black'])
    
    tplot_options, 'labflag', -1

end
