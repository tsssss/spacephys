
    event_info = _2015_0218_02_load_data()
    test = 0


    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']

    ps_time_range = time_double(['2015-02-17/23:55','2015-02-18/06:25'])

    long_time_range = time_double(['2015-02-17/22:00','2015-02-18/08:00'])
    rbsp_read_en_spec, long_time_range, probe=probe
    rbsp_read_pa_spec, long_time_range, probe=probe
    rbsp_read_orbit, long_time_range, probe=probe
    get_data, prefix+'r_gse', times, r_gse
    dis_var = prefix+'dis'
    store_data, dis_var, times, snorm(r_gse)
    add_setting, dis_var, /smart, dictionary($
        'short_name', 'RBSP-'+strupcase(probe), $
        'yrange', [1,6], $
        'ytickv', [2,4,6], $
        'yticks', 2, $
        'yminor', 2, $
        'unit', 'Re', $
        'display_type', 'scalar' )

    rgb = constant('rgb')
    xyz = constant('xyz')
    fac_labels = ['para','west','out']
    label_size = 0.8


    tvar = prefix+'e_en_spec'
    yrange = [15,4e4]
    yminor = 10
    log_ytickv = make_bins(alog10(yrange), 1, /inner)
    ytickv = 10^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(alog10(ytickv),format='(I0)')
    options, tvar, 'ytitle', 'Energy!C(eV)'
    options, tvar, 'yrange', yrange
    options, tvar, 'yticks', yticks
    options, tvar, 'ytickv', ytickv
    options, tvar, 'yminor', yminor
    options, tvar, 'ytickname', ytickn
    options, tvar, 'zcharsize', label_size
    zrange = [1e4,1e10]
    log_zrange = alog10(zrange)
    ztickv = 10^make_bins(log_zrange, 1)
    zticks = n_elements(ztickv)
    ztickn = '10!U'+string(alog10(ztickv),format='(I0)')
    ztickn[0:*:2] = ' '
    options, tvar, 'zrange', zrange
    options, tvar, 'zticks', zticks
    options, tvar, 'ztickv', ztickv
    options, tvar, 'ztickname', ztickn



    tvar = prefix+'e_pa_spec'
    yrange = [0,180]
    ystep = 90
    yminor = 4
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    options, tvar, 'ytitle', 'Pitch!C(deg)'
    options, tvar, 'yrange', yrange
    options, tvar, 'yticks', yticks
    options, tvar, 'ytickv', ytickv
    options, tvar, 'yminor', yminor
    options, tvar, 'zcharsize', label_size
    zrange = [1e5,1e10]
    log_zrange = alog10(zrange)
    ztickv = 10^make_bins(log_zrange, 1)
    zticks = n_elements(ztickv)
    ztickn = '10!U'+string(alog10(ztickv),format='(I0)')
    ztickn[0:*:2] = ' '
    options, tvar, 'zrange', zrange
    options, tvar, 'zticks', zticks
    options, tvar, 'ztickv', ztickv
    options, tvar, 'ztickname', ztickn


    plot_file = join_path([homedir(),'2015_0218_fig_overview1.pdf'])
    if test eq 1 then plot_file = 0
    sgopen, plot_file, xsize=4, ysize=3, xchsz=xchsz, ychsz=ychsz

    vars = prefix+['e_en_spec','e_pa_spec','dis']
    nvar = n_elements(vars)
    fig_labels = letters(nvar)+'. '+[$
        'e EN','e PA','|R|']
    margins = [9,4,7,1]
    ypans = [1,1,0.6]
    poss = sgcalcpos(nvar, margins=margins, ypans=ypans)

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


    tplot, vars, trange=long_time_range, position=poss, /novtitle
    for pan_id=0,nvar-1 do begin
        tpos = poss[*,pan_id]
        tx = tpos[0]+xchsz*0.5
        tx = xchsz*1
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,/normal, fig_labels[pan_id]
    endfor

    var = prefix+'dis'
    tpos = poss[*,where(vars eq var)]
    get_data, var, xxs, yys, limits=lim
    xrange = long_time_range
    yrange = lim.yrange
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, noerase=1, nodata=1
    bar_thick = keyword_set(test)? 4: 8
    plots, time_range, yrange[0]+[0,0], thick=bar_thick
    tmp = convert_coord(time_range[0], yrange[0], /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    xyouts, tx+xchsz*0.5,ty+ychsz*0.2,/normal, 'This event', charsize=label_size


    tys = 3+[0,0]
    txs = ps_time_range
    plots, txs, tys, /data
    for ii=0,1 do begin
        tmp = convert_coord(txs[ii],tys[ii], /data, /to_normal)
        tx = tmp[0]
        ty = tmp[1]
        plots, tx+[0,0], ty+[-1,1]*ychsz*0.25, /normal
    endfor
    tmp = convert_coord(mean(txs),tys[0], /data, /to_normal)
    ty = tmp[1]+ychsz*0.2
    tx = tmp[0]
    xyouts, tx,ty,/normal, alignment=0.5, 'SC in plasma sheet', charsize=label_size

    if keyword_set(test) then stop
    sgclose
end
