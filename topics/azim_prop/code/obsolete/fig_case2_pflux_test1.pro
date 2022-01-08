;+
; Plot Poynting flux stuff.
;-

;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 1

;---Settings.
    einfo = get_var_data(info_var)
    time_range = einfo.time_range_plot1
    scale_info = {s0:12d, s1:1000, dj:1d/8, ns:0d}

    fac = ['||','east','north']
    rgb = sgcolor(['red','green','blue'])
    probes = ['a','d','e']
    pres = 'th'+probes+'_'
    pf_unit = 'mW/m!U2!N'

;---Calculate pflux.
    foreach pre0, pres do begin
        ; Uniform time.
        vars = pre0+['b0_gsm','db_gsm','r_gsm']
        foreach var, vars do interp_time, var, to=pre0+'e_gsm'

        ; Convert E/B fields to FAC.
        define_fac, pre0+'b0_gsm', pre0+'r_gsm'
        vars = pre0+['db_gsm','e_gsm']
        foreach var, vars do to_fac, var
        
        ; Calculate pflux.
        stplot_calc_pflux_mor, pre0+'e_fac', pre0+'db_fac', pre0+'pf_fac', scaleinfo=scale_info, /power
        
        get_data, pre0+'pf_fac', times, pf_fac
        store_data, pre0+'pf_para', times, pf_fac[*,0]
        get_data, pre0+'cmap', uts, cmap
        cmap = sinterpol(cmap, uts, times)
        store_data, pre0+'pf_para_map', times, pf_fac[*,0]*cmap
        add_setting, pre0+'pf_para_map', /smart, {$
            display_type: 'scalar', $
            unit: pf_unit, $
            short_name: 'S!D||!N'}
    endforeach
    
    vars = pres+'pf_fac'
    options, vars, 'colors', rgb
    options, vars, 'labels', 'FAC S!D'+fac+'!N'
    
    sc1s = 'th'+['d','e','a']
    vars = sc1s+'_pf_para_map'
    options, vars, 'constant', 0
    vars = [vars,sc1s+'_dbmag']
    tplot, vars, trange=time_range
    stop
    
    sc2s = 'th'+['d','e']
    nsc2 = n_elements(sc2s)
    pf_times = time_double(['2014-08-28/10:09:06','2014-08-28/10:10:51'])
    pf_mlons = fltarr(nsc2)
    foreach sc, sc2s, ii do pf_mlons[ii] = get_var_data(sc+'_fmlon_t89', at=pf_times[ii])
    print, pf_mlons
    print, (pf_mlons[1]-pf_mlons[0])/(pf_times[1]-pf_times[0])
end
