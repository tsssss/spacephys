function fig_kaw_v01_ebratio_panel, time_range, position=tpos, $
    e_mor_var=e_mor_var, b_mor_var=b_mor_var, event_info=event_info, show_ytick=show_ytick


    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]
    label_size = 0.8

    dp_tr = time_range
    prefix = get_prefix(e_mor_var)

    pow = get_var_data(e_mor_var, fs, times=times, in=dp_tr, limits=lim)
    ntime = n_elements(times)
    e_gws = total(pow,1)/ntime
    pow = get_var_data(b_mor_var, fs, times=times, in=dp_tr, limits=lim)
    ntime = n_elements(times)
    b_gws = total(pow,1)/ntime
    ebratio = sqrt(e_gws/b_gws)*1e3

    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    window0 = ((event_info['rbsp'])['rbspa'])['b0_window']
    xrange = [20,2e5]
    yrange = get_setting(b_mor_var, 'yrange')
    yrange = minmax(1d3/[1d/8,window0])
    log_yrange = alog10(yrange)
    log_ytickv = make_bins(log_yrange,1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    yminor = 9
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    index = where(log_ytickv eq 0, count)
    if count ne 0 then ytickn[index] = '1'
    index = where(log_ytickv eq 1, count)
    if count ne 0 then ytickn[index] = '10'
    if keyword_set(show_ytick) then begin
        ytitle = 'Freq!C(mHz)'
        ytickformat = ''
    endif else begin
        ytitle = ''
        ytickformat = '(A1)'
    endelse

    
    xtitle = 'E/B (km/s)'

    plot, xrange, yrange, $
        xstyle=1, ystyle=1, $
        xrange=xrange, yrange=yrange, $
        xlog=1, ylog=1, $
        position=tpos, nodata=1, noerase=1, $
        xticklen=xticklen, yticklen=yticklen, $
        ytickname=ytickn, ytickv=ytickv, yminor=yminor, yticks=yticks, $
        xtitle=xtitle, ytickformat=ytickformat, ytitle=ytitle

    yys = fs
    xxs = ebratio
    oplot, xxs, yys, color=sgcolor('silver'), linestyle=1

    f_spin = 1d3/11
    oplot, xrange, f_spin+[0,0], linestyle=1
    if keyword_set(show_ytick) then begin
        tmp = convert_coord(xrange[0],f_spin, data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]
        msg = 'f!Dspin!N'
        xyouts, tx,ty,msg, normal=1;, charsize=label_size, alignment=0
    endif
    index = where_pro(fs,'[]',f_spin*[1d/3,6], count=count)
    if count ne 0 then begin
        ebratio[index] = !values.f_nan
    endif

    yys = fs
    xxs = ebratio
    oplot, xxs, yys;, color=sgcolor('light_salmon')
    
    

    va = mean(get_var_data(prefix+'va', times=times, in=dp_tr))
    oplot, va+[0,0], yrange, linestyle=1, color=sgcolor('red')
    if keyword_set(show_ytick) then begin
        tmp = convert_coord(va,yrange[1], data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]+ychsz*0.5
        msg = 'v!DA!N'
        xyouts, tx,ty,msg, normal=1, alignment=0.5, color=sgcolor('red')
    endif
    

    f_sc = fs*1e-3     ; in Hz.
    f_g0 = 1.6e-19*1e-9/1.67e-27/2/!dpi   ; in Hz.
    bmag = median(snorm(get_var_data(prefix+'b_gsm', in=dp_tr)))
    avg_mass = median(get_var_data(prefix+'avg_mass', in=dp_tr))
    f_gi = f_g0*bmag/avg_mass
    tavg = median(get_var_data(prefix+'p_t', in=dp_tr))
    vi = sqrt(tavg*1.6e-19/(avg_mass*1.67e-27))*1e-3    ; in km/s
    ;vf = median(snorm(get_var_data(prefix+'u_gsm', in=dp_tr)))
    ravg = median(snorm(get_var_data(prefix+'r_gsm', in=dp_tr)))
    vf = 50d
    ebr_theory = va*sqrt(1+(f_sc/f_gi*(vi/vf))^2)
    oplot, ebr_theory, fs, color=sgcolor('red'), linestyle=0
    if keyword_set(show_ytick) then begin
        tmp = convert_coord(va,yrange[0], data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.5
        msg = 'KAW'
        xyouts, tx,ty,msg, normal=1, alignment=0, color=sgcolor('red')
    endif

    return, 1

end

function fig_kaw_v01, event_info=event_info, test=test

    version = 'v01'
    id = '2015_0317'
    if n_elements(event_info) eq 0 then event_info = low_lshell_outflow_load_data(id)
    label_size = 0.8
    tmp = get_abs_chsz()
    abs_xchsz = tmp[0]
    abs_ychsz = tmp[1]
    uniform_ticklen = -abs_ychsz*0.15

    plot_dir = event_info.plot_dir
    plot_file = join_path([plot_dir,'fig_kaw_'+id+'_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    if keyword_set(test) then thick = 4 else thick = 16

    time_range = time_double(['2015-03-17/06:00','2015-03-17/15:00'])
    probe = 'a'
    prefix = 'rbsp'+probe+'_'
    
    kaw_tr_list = list()
    kaw_tr_list.add, time_double(['2015-03-17/08:10','2015-03-17/08:25'])
    kaw_tr_list.add, time_double(['2015-03-17/08:30','2015-03-17/09:00'])
    kaw_tr_list.add, time_double(['2015-03-17/11:05','2015-03-17/11:20'])
    kaw_tr_list.add, time_double(['2015-03-17/11:30','2015-03-17/11:40'])
    nkaw_tr = n_elements(kaw_tr_list)

    e_var = prefix+'edot0_fac_this'
    copy_data, prefix+'edot0_fac', e_var
    data = get_var_data(e_var, times=times, in=time_range)
    store_data, e_var, times, data

    b_var = prefix+'b1_fac_this'
    copy_data, prefix+'b1_fac', b_var
    data = get_var_data(b_var, times=times, in=time_range)
    store_data, b_var, times, data


    e_vars = stplot_split(e_var)
    yrange = [-1,1]*150
    ystep = 100
    yminor = 4
    foreach var, e_vars do begin
        set_ytick, var, yrange=yrange, ystep=ystep, yminor=yminor
    endforeach
    options, e_vars, 'constant', [-2,-1,0,1,2]*50
    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    foreach var, e_vars, id do options, var, 'labels', 'E!D'+fac_labels[id]+'!N'

    the_e_var = e_vars[2]
    e_mor_var = the_e_var+'_mor'
    if check_if_update(e_mor_var) then begin
        e_mor_var = stplot_mor_new(the_e_var, scale_info=scale_info)
        var = e_mor_var
        options, var, 'ztitle', '(mV/m)!U2!N'
        options, var, 'zrange', [1,1e4]
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        options, var, 'yrange', lim.yrange*1e3
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', 10^[1d,2,3]
    endif


    b_vars = stplot_split(b_var)
    the_b_var = b_vars[1]
    b_mor_var = the_b_var+'_mor'
    if check_if_update(b_mor_var) then begin
        b_mor_var = stplot_mor_new(the_b_var, scale_info=scale_info)
        var = b_mor_var
        options, var, 'ztitle', '(nT)!U2!N'
        ;options, var, 'zrange', [1,1e4]
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        options, var, 'yrange', lim.yrange*1e3
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', 10^[1d,2,3]
    endif
    
    
    ; Load density.
    all_species = ['e','o','p','he']
    foreach species, all_species do tmp = rbsp_read_hope_moments(time_range, probe=probe, species=species)
    vars = prefix+all_species+'_n'
    options, vars, 'ylog', 1
    o_dens_var = prefix+'o_n'
    p_dens_var = prefix+'p_n'
    avg_mass_var = prefix+'avg_mass'
    o_dens = get_var_data(o_dens_var, times=times)
    p_dens = get_var_data(p_dens_var, times=times)
    ion_dens = (o_dens+p_dens)
    o_ratio = o_dens/ion_dens
    p_ratio = p_dens/ion_dens
    o_mass = 16d
    p_mass = 1d
    avg_mass = p_ratio*p_mass+o_ratio*o_mass
    store_data, avg_mass_var, times, avg_mass

    va_var = prefix+'va'
    avg_mass = get_var_data(avg_mass_var, times=times)
    density = get_var_data(prefix+'density', at=times)
    mhd_rho = density*avg_mass
    bmag = snorm(get_var_data(prefix+'b0_gsm', at=times))
    va0 = 22.0d     ; km/s, for B in nT, n in cc, m in atomic mass.
    va = va0*bmag/sqrt(mhd_rho)
    store_data, va_var, times, va

    fig_size = [6d,4]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    margins = [10d,4,6,1]
    ypans = [[1,1]*0.7,1]
    ypads = [0.4,4]
    nypan = 3
    fig_labels = letters(nypan)
    all_poss = sgcalcpos(nypan, ypans=ypans, ypad=ypads, margins=margins)

    poss = all_poss[*,0:1]
    plot_vars = e_vars[1:2]
    foreach pid, [0,1] do begin
        tpos = poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        var = plot_vars[pid]
        options, var, 'xticklen', xticklen
        options, var, 'yticklen', yticklen
        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*0.7
        msg = fig_labels[pid]+')'
        xyouts, tx,ty, msg, normal=1
    endforeach

    tplot, plot_vars, position=poss, noerase=1, trange=time_range
    tpos = poss[*,-1]
    yrange = [0,1]
    set_axis, plot_vars[-1], position=tpos, xrange=time_range, yrange=yrange
    for pid=0, nkaw_tr-1 do begin
        tr = time_double(kaw_tr_list[pid])
        plots, tr, yrange[0]+[0,0], data=1, color=sgcolor('red'), thick=thick
    endfor
    
    tpos = poss[*,0]
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = strupcase('rbsp-'+probe)
    xyouts, tx,ty,msg, normal=1


    low_margins = [8,0,3,0]
    low_pos = all_poss[*,nypan-1]
    low_pos[0] = 0 & low_pos[2] = 1
    low_poss = sgcalcpos(1,nkaw_tr, region=low_pos, xpad=1, margins=low_margins)
    for pid=0, nkaw_tr-1 do begin
        tpos = low_poss[*,pid]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = fig_labels[nypan-1]+'-'+string(pid+1,format='(I0)')+')'
        xyouts, tx,ty,msg, normal=1
        show_ytick = (pid eq 0)? 1: 0
        tmp = fig_kaw_v01_ebratio_panel(kaw_tr_list[pid], position=tpos, $
            e_mor_var=e_mor_var, b_mor_var=b_mor_var, event_info=event_info, show_ytick=show_ytick)
    endfor


    if keyword_set(test) then stop
    sgclose
    return, plot_file


end

print, fig_kaw_v01(event_info=event_info, test=0)
end
