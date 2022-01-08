

    event_info = _2015_0218_02_load_data()
    test = 1


    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']

    rgb = constant('rgb')
    xyz = constant('xyz')
    fac_labels = ['para','west','out']
    label_size = 0.8


;---Density.
    efw_density = get_var_data(prefix+'efw_density', times=times, limits=limits)
    ;hope_density = get_var_data(prefix+'ele_n', at=times)
    hope_density = get_var_data(prefix+'e_density', at=times)
    tvar = prefix+'combo_density'
    store_data, tvar, times, [[efw_density],[hope_density]], limits=limits
    options, tvar, 'labels', ['UH line','HOPE']
    options, tvar, 'colors', sgcolor(['black','silver'])
    options, tvar, 'yrange', [1,20]
    options, tvar, 'ylog', 1



;---Bulk flow.
    o_dens = get_var_data(prefix+'o_density', in=time_range)
    p_dens = get_var_data(prefix+'p_density', in=time_range)
    o_num_dens_ratio = o_dens/(o_dens+p_dens)
    p_num_dens_ratio = 1-o_num_dens_ratio
    num_dens = get_var_data(prefix+'efw_density', in=time_range)
    p_v_gsm = get_var_data(prefix+'p_v_gsm', times=times)
    o_v_gsm = get_var_data(prefix+'o_v_gsm', times=times)
    p_mass = 1d
    o_mass = 16d
    avg_mass = p_num_dens_ratio*p_mass+o_num_dens_ratio*o_mass
    v_gsm = p_v_gsm
    for ii=0,2 do begin
        v_gsm[*,ii] = (p_num_dens_ratio*p_mass*p_v_gsm[*,ii]+$
            o_num_dens_ratio*o_mass*o_v_gsm[*,ii])/avg_mass
    endfor
    mhd_v_gsm_var = prefix+'mhd_v_gsm'
    store_data, mhd_v_gsm_var, times, v_gsm
    add_setting, mhd_v_gsm_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'MHD V', $
        'unit', 'km/s', $
        'coord', 'GSM', $
        'coord_labels', xyz )
    mhd_v_fac_var = prefix+'mhd_v_fac'
    to_fac, mhd_v_gsm_var, to=mhd_v_fac_var, q_var=prefix+'q_gsm2fac'

    yrange = [-30,90]
    ystep = 50
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 5
    options, mhd_v_fac_var, 'constant', 0
    options, mhd_v_fac_var, 'yrange', yrange
    options, mhd_v_fac_var, 'ytickv', ytickv
    options, mhd_v_fac_var, 'yticks', yticks
    options, mhd_v_fac_var, 'yminor', yminor
    options, mhd_v_fac_var, 'labels', fac_labels



;---B field related.
    b_gsm = get_var_data(prefix+'b_gsm', times=times)
    b_sm = cotran(b_gsm, times, 'gsm2sm')

    tilt = azim_df_calc_tilt(b_sm)
    tilt_var = prefix+'tilt'
    store_data, tilt_var, times, tilt

    yrange = [40,60]
    ystep = 10.
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 5
    options, tilt_var, 'yrange', yrange
    options, tilt_var, 'ytickv', ytickv
    options, tilt_var, 'yticks', yticks
    options, tilt_var, 'yminor', yminor
    options, tilt_var, 'ytitle', '(deg)'
    options, tilt_var, 'labels', 'B tilt'

    ; |B|.
    bmag = snorm(b_gsm)
    bmag_var = prefix+'bmag'
    store_data, bmag_var, times, bmag
    yrange = [55,85]
    ystep = 20.
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 4
    options, bmag_var, 'yrange', yrange
    options, bmag_var, 'ytickv', ytickv
    options, bmag_var, 'yticks', yticks
    options, bmag_var, 'yminor', yminor
    options, bmag_var, 'ytitle', '(nT)'
    options, bmag_var, 'labels', '|B|'


    ; Temperature.
    var = prefix+'e_temp'
    options, var, 'yrange', [5e2,2e3]
    options, var, 'labels', 'T ele'
    
    var = prefix+'p_temp'
    options, var, 'yrange', [2e3,2e4]

    var = prefix+'o_temp'
    options, var, 'yrange', [2e3,2e4]
    
    
    
    ; Total magnetic pressure.
    bmag = get_var_data(bmag_var, in=time_range, times=times)
    mu0 = 4*!dpi*1e-7
    pmag = (bmag*1e-9)^2/(2*mu0)*1e9
    store_data, prefix+'p_mag', times, pmag, $
        limits = {ytitle:'(nPa)', labels: 'P mag', yrange: [0,5] }
    
    
    ; Total thermal pressure.
    density = get_var_data(prefix+'efw_density', at=times)
    p_density = get_var_data(prefix+'p_density', at=times)
    o_density = get_var_data(prefix+'o_density', at=times)
    p_ratio = p_density/(p_density+o_density)
    o_ratio = 1-p_ratio
    temp = get_var_data(prefix+'e_temp', at=times)+$
        get_var_data(prefix+'p_temp', at=times)*p_ratio+$
        get_var_data(prefix+'o_temp', at=times)*o_ratio
    store_data, prefix+'t_total', times, temp, limits={$
        ytitle:'(eV)', $
        ylog: 1 }
    pthermal = calc_thermal_pressure(density, temp)


    pthermal = calc_thermal_pressure($
        get_var_data(prefix+'efw_density', at=times), $
        get_var_data(prefix+'e_temp', at=times))

    ;pthermal = calc_thermal_pressure($
    ;    get_var_data(prefix+'e_density', at=times), $
    ;    get_var_data(prefix+'e_temp', at=times))
        
    foreach species, ['p','o'] do begin
        pthermal += calc_thermal_pressure($
            get_var_data(prefix+species+'_density', at=times), $
            get_var_data(prefix+species+'_temp', at=times))
    endforeach
    store_data, prefix+'p_thermal', times, pthermal, $
        limits = {ytitle:'(nPa)', labels: 'P thermal', yrange:[0,5] }
    
    
    ; Total pressure.
    store_data, prefix+'p_total', times, pmag+pthermal, $
        limits = {ytitle:'(nPa)', labels:'P total', yrange: [0,5] }
    
    beta = pthermal/pmag
    store_data, prefix+'beta', times, beta, $
        limits= {ytitle:'(#)', labels: 'beta' }
    


    plot_file = join_path([homedir(),'2015_0218_fig_overview3.pdf'])
    if test eq 1 then plot_file = 0
    sgopen, plot_file, xsize=4, ysize=6, xchsz=xchsz, ychsz=ychsz

    vars = prefix+['combo_density','tilt','bmag','e_temp','p_total']
    nvar = n_elements(vars)
    fig_labels = letters(nvar)+'. '+['Density','Tilt','|B|','T ele','P total']
    margins = [8,4,7,1]
    poss = sgcalcpos(nvar, margins=margins)


    xticklen_chsz = -0.2
    yticklen_chsz = -0.3
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 4
    tplot_options, 'labflag', -1
    options, vars, 'xticklen', xticklen
    options, vars, 'yticklen', yticklen


    tplot, vars, trange=time_range, position=poss, /novtitle
    for pan_id=0,nvar-1 do begin
        tpos = poss[*,pan_id]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*0.9
        xyouts, tx,ty,/normal, fig_labels[pan_id]
    endfor


    if keyword_set(test) then stop
    sgclose
end
