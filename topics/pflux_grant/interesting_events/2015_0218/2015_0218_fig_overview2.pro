
    event_info = _2015_0218_02_load_data()
    test = 1


    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']

    rgb = constant('rgb')
    xyz = constant('xyz')
    fac_labels = ['para','west','out']
    label_size = 0.8


;---pflux.
    tvar = prefix+'pf_fac_norm'
    get_data, prefix+'pfdot0_fac', times, data
    get_data, prefix+'cmap', tuts, cmap
    cmap0 = mean(cmap[where(tuts ge time_range[0] and tuts le time_range[1])])
    store_data, tvar, times, data*cmap0, limits={$
        ytitle:'(mW/m!U2!N)', colors:rgb}

    options, tvar, 'ystyle', 1
    options, tvar, 'yrange', [-15,5]
    options, tvar, 'ytickv', [-3,-1,1]*5
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'constant', 0
    options, tvar, 'labels', fac_labels
    options, tvar, 'labflag', -1


;---Density.
    efw_density = get_var_data(prefix+'efw_density', times=times, limits=limits)
    hope_density = get_var_data(prefix+'ele_n', at=times)
    tvar = prefix+'combo_density'
    store_data, tvar, times, [[efw_density],[hope_density]], limits=limits
    options, tvar, 'labels', ['UH line','HOPE']
    options, tvar, 'colors', sgcolor(['black','silver'])
    options, tvar, 'yrange', [1,20]
    options, tvar, 'ylog', 1



;---Bulk flow.
    o_dens = get_var_data(prefix+'o_density', in=time_range, times=times)
    p_dens = get_var_data(prefix+'p_density', in=time_range)
    o_num_dens_ratio = o_dens/(o_dens+p_dens)
    p_num_dens_ratio = 1-o_num_dens_ratio
    num_dens = get_var_data(prefix+'efw_density', in=time_range)
    p_v_gsm = get_var_data(prefix+'p_v_gsm', at=times)
    o_v_gsm = get_var_data(prefix+'o_v_gsm', at=times)
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

    yrange = [-30,30]
    ystep = 30
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 3
    options, mhd_v_fac_var, 'constant', 0
    options, mhd_v_fac_var, 'yrange', yrange
    options, mhd_v_fac_var, 'ytickv', ytickv
    options, mhd_v_fac_var, 'yticks', yticks
    options, mhd_v_fac_var, 'yminor', yminor
    options, mhd_v_fac_var, 'labels', fac_labels



;---EWOgram.
    mlt_range = [-2.5,-0.5]
    mlat_range = [61.5,64.5]
    zrange = [60,180]

    ewo_var = 'thg_asf_ewo'
    yrange = mlt_range
    ystep = 1
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 5
    options, ewo_var, 'yrange', yrange
    options, ewo_var, 'yticks', yticks
    options, ewo_var, 'ytickv', ytickv
    options, ewo_var, 'yminor', yminor
    options, ewo_var, 'zrange', zrange
    options, ewo_var, 'ytitle', 'MLT (deg)'


    keo_var = 'thg_asf_keo'
    yrange = mlat_range
    ystep = 1
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 4
    options, keo_var, 'yrange', yrange
    options, keo_var, 'yticks', yticks
    options, keo_var, 'ytickv', ytickv
    options, keo_var, 'yminor', yminor
    options, keo_var, 'zrange', zrange
    options, keo_var, 'ytitle', 'MLat (deg)'


    plot_file = join_path([homedir(),'2015_0218_fig_overview2.pdf'])
    if test eq 1 then plot_file = 0
    sgopen, plot_file, xsize=4, ysize=6, xchsz=xchsz, ychsz=ychsz

    vars = [prefix+['e_en_spec','e_pa_spec','combo_density','mhd_v_fac',$
        'pf_fac_norm'], 'thg_asf_'+['ewo','keo']]
    nvar = n_elements(vars)
    fig_labels = letters(nvar)+'. '+[$
        'e EN','e PA','Density','MHD Velocity','S norm. 100 km','EWOgram','KEOgram']
    margins = [8,4,7,1]
    poss = sgcalcpos(nvar, margins=margins)

    asf_vars = [ewo_var,keo_var]
    options, asf_vars, 'color_table', 40
    cbpos = poss[*,-1]
    cbpos[0] = poss[2,-1]+xchsz*0.8
    cbpos[2] = cbpos[0]+xchsz*0.8
    cbpos[3] = poss[3,-2]
    options, asf_vars[1], 'zposition', cbpos
    options, asf_vars[0], 'no_color_scale', 1
    options, asf_vars, 'zcharsize', label_size


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
        if pan_id eq 0 or pan_id eq  4 or pan_id eq 5 then ty = tpos[1]+ychsz*0.3
        if pan_id eq nvar-1 or pan_id eq nvar-2 then begin
            color = sgcolor('white')
        endif else color = sgcolor('black')
        xyouts, tx,ty,/normal, fig_labels[pan_id], color=color
    endfor

    foreach var, vars[-2:-1] do begin
        get_data, var, limits=lim
        xrange = time_range
        yrange = lim.yrange
        tpos = reform(poss[*,where(vars eq var)])
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        pos_var = (var eq 'thg_asf_ewo')? prefix+'fmlt_t89': prefix+'fmlat_t89'
        get_data, pos_var, xxs, yys
        oplot, xxs, yys, color=sgcolor('grey')
    endforeach

    if keyword_set(test) then stop
    sgclose
end
