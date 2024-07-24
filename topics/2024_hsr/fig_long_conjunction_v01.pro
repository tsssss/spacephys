
function fig_long_conjunction_v01, event_info=event_info, test=test

    id = '2015_0317'
    version = 'v02'
    if n_elements(event_info) eq 0 then event_info = low_lshell_outflow_load_data(id)

    

    time_range = time_double(['2015-03-17/03:10','2015-03-17/10:30'])
    time_range = time_double(['2015-03-17/06:00','2015-03-17/08:30'])
    kaw_times = time_double('2015-03-17/'+['06:25','06:35','06:51','07:48','08:53'])

    ps_times = time_double(['2015-03-17/06:25','2015-03-17/06:35'])
    
    probe = 'b'
    prefix = 'rbsp'+probe+'_'
    external_models = ['t89','t96','t01','t04s']
    model_colors = sgcolor(['black','blue','green','purple'])
    


;---KEO.
    keo_var = 'thg_asf_keo'
    ytitle = 'MLat!C(deg)'
    model_suffix = '_dipole_t01_north'
    mlt_image_var = 'thg_asf_mlt_image_rect'

    mlt_images = get_var_data(mlt_image_var, times=times, settings=settings)
    ntime = n_elements(times)

    mlt_bins = settings['mlt_bins']
    mlat_bins = settings['mlat_bins']
    
    fmlt_var = prefix+'fmlt'+model_suffix
    fmlt_var = prefix+'mlt'
    fmlt = get_var_data(fmlt_var, at=times)+1
    fmlt_del = [-1,1]*0.2
    fmlt_del = [-1,1]
    nbin = n_elements(mlat_bins)
    keo = fltarr(ntime,nbin)
    foreach time, times, time_id do begin
        index = where_pro(mlt_bins,'[]',fmlt[time_id]+fmlt_del,count=count)
        if count eq 0 then continue
        keo[time_id,*] = total(mlt_images[time_id,index,*],2)/count
    endforeach
    index = where(keo eq 0, count)
    if count ne 0 then keo[index] = 0.1
    

    store_data, keo_var, times, keo, mlat_bins, limits={ytitle:ytitle,spec:1,zlog:0,color_table:49}
    set_ytick, keo_var, yrange=[59,71], ystep=5
    options, keo_var, ztickv=[100,1e3], zticklen=-0.5, $
        zticks=1, zminor=9, ztitle='Count (#)', ztickname=['10!U2','10!U3'], zrange=[50,5e3], zlog=1


;---EN spec.
    e_en_var = rbsp_read_en_spec(time_range, probe=probe, species='e', update=1)
    o_en_var = rbsp_read_en_spec(time_range, probe=probe, species='o')
    o_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='o', errmsg=errmsg, energy_range=[50,5e4], update=1)

    
;---E field spec.
    id = 'e'
    out_var = prefix+'e_spec'
    spec_vars = list()
    spec_vars.add, rbsp_read_wave_spec_mhz(time_range, probe=probe)
    spec_vars.add, rbsp_read_wave_spec_khz(time_range, probe=probe, id=id)
    spec_vars.add, rbsp_read_wave_spec_hz(time_range, probe=probe, id=id, use_e34=1, update=0)
    freq_ranges = [[1e4,5e4],[10,1e4],[0.1,10]]

    ;---Combine them.
    yrange = [0.1,5e4]
    dfreq = 1.15
    freqs = smkgmtrc(yrange[0],yrange[1],dfreq, 'dx')
    nfreq = n_elements(freqs)
    xrange = time_range
    time_step = 1d
    common_times = make_bins(xrange,time_step)
    ntime = n_elements(common_times)
    specs = fltarr(ntime,nfreq)
    foreach var, spec_vars, var_id do begin
        if var eq '' then continue
        data = get_var_data(var, vals, at=common_times, limits=lim)
        freq_range = freq_ranges[*,var_id]
        index = where_pro(freqs, '[)', freq_range, count=count)
        if count eq 0 then message, 'Inconsistency ...'
        specs[*,index] = transpose(sinterpol(transpose(data),vals,freqs[index]))
    endforeach
    unit = lim.unit
    zrange = [1e-8,1e1]

    log_ytickv = make_bins(minmax(alog10(yrange)),1,inner=1)
    ytickv = 10d^log_ytickv
    ytickname = get_short_log_tickname(log_ytickv)
    yminor = 9

    store_data, out_var, common_times, specs, freqs
    add_setting, out_var, smart=1, dictionary($
        'requested_time_range', time_range, $
        'no_interp', 1, $
        'display_type', 'spec', $
        'unit', unit, $
        'ytitle', 'Freq (Hz)', $
        'yrange', yrange, $
        'ylog', 1, $
        'ytickv', ytickv, $
        'ytickname', ytickname, $
        'yminor', yminor, $
        'zlog', 1, $
        'zrange', zrange, $
        'short_name', 'E' )
    e_spec_var = out_var
    
    ; Cyclotron freqs.
    fc_vars = list()
    foreach species, ['e','o','he','p'] do fc_vars.add, rbsp_read_gyro_freq(time_range, probe=probe, species=species)
    var = prefix+'fce_half'
    fce = get_var_data(prefix+'fce', times=times)
;        store_data, var, times, fce*0.5
;        fc_vars.add, var
    options, prefix+'fce', labels='f!Dc,e'
    options, prefix+'fcp', labels='f!Dc,H'
    options, prefix+'fco', labels='f!Dc,O'
    options, prefix+'fche', labels='f!Dc,He'

    var = prefix+'flh'
    fcp = get_var_data(prefix+'fcp', times=times)
    store_data, var, times, fcp*43, limits={labels:'f!DLH!N'}
    fc_vars.add, var
    fc_vars = fc_vars.toarray()
    fc_colors = get_color(n_elements(fc_vars))
    foreach var, fc_vars, ii do options, var, 'colors', fc_colors[ii]

    e_spec_combo = e_spec_var+'_combo'
    store_data, e_spec_combo, data=[e_spec_var,fc_vars]
    options, e_spec_combo, 'yrange', get_setting(e_spec_var,'yrange')
    options, e_spec_combo, 'labels', ''


    b_spec_var = prefix+'b_spec_hz'
    b_spec_combo = b_spec_var+'_combo'
    store_data, b_spec_combo, data=[b_spec_var,fc_vars]
    options, b_spec_combo, 'yrange', get_setting(b_spec_var,'yrange')
    options, b_spec_combo, 'labels', ''
    

;---tilt angle.
    bmod_gsm_vars = prefix+'bmod_gsm_igrf_'+external_models
    b_gsm_var = prefix+'b_gsm'

    b_tilt_vars = list()
    foreach var, [b_gsm_var,bmod_gsm_vars] do begin
        var_out = streplace(var,'gsm','tilt')
        b_tilt_vars.add, lets_calc_vec_elev(var, coord='sm', var_info=var_out)
    endforeach
    b_tilt_vars = b_tilt_vars.toarray()

    labels = ['Obs',strupcase(external_models)]
    colors = [sgcolor('red'),model_colors]
    nmodel = n_elements(external_models)
    b_tilt_combo_var = prefix+'b_tilt_combo'
    b_tilt_combo_var = stplot_merge(b_tilt_vars, output=b_tilt_combo_var, labels=labels[0:nmodel], colors=colors[0:nmodel])
    options, b_tilt_combo_var, ytitle='(deg)'
    add_setting, b_tilt_combo_var, smart=1, dictionary($
        'labels', labels, $
        'colors', colors, $
        'constant', [20,30,40], $
        'yrange', [15,48], $
        'ytickv', [20,30,40], $
        'yticks', 2, $
        'yminor', 4 )

        
;    ; magnetic pressure.
;    pmag_var = rbsp_read_magnetic_pressure(b_var=b_gsm_var)
;    options, pmag_var, 'ylog', 1
;    ; thermal_pressure.
;    pthermal_var = prefix+'p_thermal'
;    foreach species, ['e','p','o'] do begin
;        vars = rbsp_read_hope_moments(time_range, probe=probe, species=species)
;    endforeach
;    
;    species = 'e'
;    p_thermal = calc_thermal_pressure($
;        get_var_data(prefix+species+'_n',times=times), $
;        get_var_data(prefix+species+'_t',at=times))
;    foreach species, ['p','o'] do begin
;        p_thermal += calc_thermal_pressure($
;            get_var_data(prefix+species+'_n',at=times), $
;            get_var_data(prefix+species+'_t',at=times))
;    endforeach
;    store_data, pthermal_var, times, p_thermal
;    add_setting, pthermal_var, smart=1, dictionary($
;        'ylog', 1, $
;        'display_type', 'scalar', $
;        'short_name', 'P!Dth', $
;        'unit', 'nPa' )
;    
;    
;    p_thermal = get_var_data(pthermal_var, times=times)
;    p_mag = get_var_data(pmag_var, at=times)
;    beta = p_thermal/p_mag
;    beta_var = prefix+'beta'
;    store_data, beta_var, times, beta
;    add_setting, beta_var, smart=1, dictionary($
;        'ylog', 1, $
;        'yrange', [1e-4,1e-1]*7, $
;        'constant', [1e-2,1e-1], $
;        'display_type', 'scalar', $
;        'short_name', tex2str('beta'), $
;        'unit', '#' )
    
    
    ; FMLat.
;    fmlat_north_var = prefix+'fmlat'+suffix+'_north'
;    fmlat_south_var = prefix+'fmlat'+suffix+'_south'
;    fmlat_var = prefix+'fmlat'
;    fmlat_var = stplot_merge([fmlat_north_var,fmlat_south_var], output=fmlat_var)
;    get_data, fmlat_var, times, data
;    store_data, fmlat_var, times, abs(data)
;    add_setting, fmlat_var, smart=1, dictionary($
;        'display_type', 'stack', $
;        'ytitle', '(deg)', $
;        'yrange', [40,70], $
;        'ytickv', [50,60], $
;        'yticks', 1, $
;        'yminor', 2, $
;        'constant', [50,60], $
;        'labels', ['North','South'], $
;        'colors', sgcolor(['red','blue']) )



    ; dis.
    r_gsm_var = prefix+'r_gsm'
    dis_var = prefix+'dis'
    get_data, r_gsm_var, times, data
    store_data, dis_var, times, snorm(data)
    add_setting, dis_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', '|R|', $
        'unit', 'Re' )
    var = dis_var
    options, var, 'yrange', [1,6]
    options, var, 'ytickv', [2,4,6]
    options, var, 'yticks', 2
    options, var, 'yminor', 2
    options, var, 'constant', [2,4]


    ; MLat.
    mlat_var = rbsp_read_mlat(time_range, probe=probe)
    options, mlat_var, 'yrange', [-1,1]*20
    options, mlat_var, 'constant', [-1,0,1]*10
    options, mlat_var, 'ytickv', [-1,0,1]*10
    options, mlat_var, 'yticks', 2
    options, mlat_var, 'yminor', 2


    ; MLT.
    mlt_var = prefix+'mlt'
    yrange = [-1,1]*12
    ytickv = [-1,0,1]*6
    yticks = n_elements(ytickv)-1
    yminor = 3
    constants = [-1,0,1]*6
    options, mlt_var, 'yrange', yrange
    options, mlt_var, 'constant', constants
    options, mlt_var, 'ytickv', ytickv
    options, mlt_var, 'yticks', yticks
    options, mlt_var, 'yminor', yminor

    ; Tilt.
    var = b_tilt_combo_var
    options, var, 'yrange', [10,90]
    options, var, 'ytickv', [30,60]
    options, var, 'constant', [30,60]
    options, var, 'yticks', 1
    options, var, 'yminor', 3
;    var = dtheta_var
;    yrange = [-45,10]
;    ystep = 20
;    ytickv = make_bins(yrange,ystep, inner=1)
;    yticks = n_elements(ytickv)-1
;    options, var, 'ytickv', ytickv
;    options, var, 'yticks', yticks
;    options, var, 'yrange', yrange
;    options, var, 'yminor', 2
;    options, var, 'constant', ytickv

    ; E
    e_var = prefix+'edot0_fac'
    var = e_var
    options, var, 'yrange', [-1,1]*20
    options, var, 'ytickv', [-1,0,1]*10
    options, var, 'yticks', 2
    options, var, 'yminor', 4
    options, var, 'constant', [-1,0,1]*10
    
    dens_var = prefix+'density'
    e_en_var = prefix+'e_en_spec'
    p_en_var = prefix+'p_en_spec'
    o_en_var = prefix+'o_en_spec'
    b_theta_var = prefix+'b_gsm_theta'
    b_var = prefix+'b_sm'
    pthermal_var = prefix+'p_thermal'
    
    plot_vars = [keo_var,b_tilt_combo_var,e_en_var,p_en_var,o_en_var,o_pa_var,$
        e_spec_var]
    nvar = n_elements(plot_vars)
    fig_letters = letters(nvar)
    fig_labels = fig_letters+') '+['Aurora','B elev','e- EN','H+ EN','O+ EN','O+ PA','E spec']
    ypans = fltarr(nvar)+1.
    index = where(plot_vars eq e_spec_var, count)
    if count ne 0 then ypans[index] = 1.6
    index = where(plot_vars eq b_tilt_combo_var, count)
    if count ne 0 then ypans[index] = 1.2
;    index = where(plot_vars eq beta_var, count)
;    if count ne 0 then ypans[index] = 0.8
    index = where(plot_vars eq mlat_var, count)
    if count ne 0 then ypans[index] = 0.6
    index = where(plot_vars eq mlt_var, count)
    if count ne 0 then ypans[index] = 0.6
    index = where(plot_vars eq dis_var, count)
    if count ne 0 then ypans[index] = 0.6
    index = where(plot_vars eq o_pa_var, count)
    if count ne 0 then ypans[index] = 0.8
    index = where(plot_vars eq keo_var, count)
    if count ne 0 then ypans[index] = 1.2
    
    the_vars = [e_en_var,p_en_var,o_en_var]
    options, the_vars, 'ytitle', 'Energy!C(eV)'
    the_vars = [p_en_var,o_en_var]
    options, the_vars, ytickname=['10!U2','10!U3','10!U4'], $
        ytickv=[100,1000,10000], yticks=2, yminor=9, $
        yrange=[20,3e4]

    the_vars = e_en_var 
    ztickn = '10!U'+['4','5','6','7','8','9','10']
    ztickn[1:*:2] = ' '       
    options, the_vars, ytickname=['10!U2','10!U3','10!U4'], $
        ytickv=[100,1000,10000], yticks=2, yminor=9, $
        ztickv=10d^[4,5,6,7,8,9,10], zticks=6, ztickname=ztickn
        
    options, b_spec_var, ztickv=10d^[-3,-2,-1,0,1], zticks=4, $
        ztickname=['10!U-3','10!U-2','0.1','1','10']

    zticklen = -0.5
    the_vars = [e_en_var,p_en_var,o_en_var,o_pa_var,b_spec_var,e_spec_var]
    options, the_vars, zticklen=zticklen, zminor=9
    
    
    the_vars = o_pa_var
    options, the_vars, 'ytitle', 'PA!C(deg)'
    
    the_vars = [b_spec_combo,e_spec_combo]
    options, the_vars, 'ytitle', 'Freq!C(Hz)'
    
    options, mlt_var, yrange=[-1,1]*7
    options, mlat_var, yrange=[-1,1]*15
    
    pansize = [6,0.9]
    margins = [10,6,7,1]
    
    if n_elements(plot_dir) eq 0 then plot_dir = srootdir()
    plot_file = join_path([plot_dir,'fig_long_conjunction_'+id+'_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    poss = panel_pos(plot_file, ypans=ypans, fig_size=fig_size, $
        nypan=nvar, pansize=pansize, panid=[0,1], margins=margins)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz, inch=1

    xticklen_chsz = -0.25   ; in ychsz.
    yticklen_chsz = -0.40   ; in xchsz.
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        options, plot_vars[ii], 'xticklen', xticklen
        options, plot_vars[ii], 'yticklen', yticklen
    endfor
    
    
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        is_spec = get_var_setting(plot_vars[ii],'spec', exist)
        if not exist then continue
        if not is_spec then continue
        cbpos = tpos
        cbpos[0] = cbpos[2]+xchsz*0.5
        cbpos[2] = cbpos[0]+xchsz*0.8
        zticklen = -0.5
        options, plot_vars[ii], zposition=cbpos, zticklen=zticklen, zcharsize=0.8
    endfor
    
    
    ; Add PS label
    pid = where(plot_vars eq b_tilt_combo_var)
    tpos = poss[*,pid]
    tpos[1] = poss[1,pid+1]
    ps_color = sgcolor('wheat')
    xrange = time_range
    yrange = [0,1]
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
        position=tpos, nodata=1, noerase=1
    polyfill, ps_times[[0,1,1,0,0]], yrange[[0,0,1,1,0]], data=1, color=ps_color
    
    
    
    tplot_options, 'tickinterval', 600
    tplot_options, 'version', 3
    var_labels = prefix+['mlt','mlat','dis']
    options, prefix+'mlt', 'ytitle', 'MLT (h)'
    options, prefix+'mlat', 'ytitle', 'MLat (deg)'
    options, prefix+'dis', 'ytitle', '|R| (Re)'
    tplot, plot_vars, trange=time_range, position=poss, $
        noerase=1, var_label=var_labels, vlab_margin=9
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*(margins[0]-1)
        ty = tpos[3]-ychsz*0.7
        msg = fig_labels[pid]
        xyouts, tx,ty,normal=1, msg
    endfor
    
    
    foreach spec_var, [b_spec_var,e_spec_var] do begin
        pid = where(plot_vars eq spec_var, count)
        if count eq 0 then continue
        tpos = poss[*,pid]
        spec_settings = get_var_setting(spec_var)
        xrange = time_range
        yrange = spec_settings['yrange']
        ylog = spec_settings['ylog']
        plot, xrange, yrange, $
            ylog=ylog, xrange=xrange, yrange=yrange, $
            xstyle=5, ystyle=5, $
            position=tpos, noerase=1, nodata=1
        
        nfreq_var = n_elements(fc_vars)
        colors = get_color(nfreq_var)
        foreach var, fc_vars, var_id do begin
            color = colors[var_id]
            yys = get_var_data(var, times=xxs, settings=settings)
            tx = time_double('2015-03-17/06:00')
            ty = interpol(yys, xxs, tx)
            if product(ty-yrange) ge 0 then continue

            msg = settings.labels
            color = settings.colors
            oplot, xxs, yys, color=color
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]-xchsz*1.5
            ty = tmp[1]-ychsz*0.35
            tx = tmp[0]+xchsz*1
            ty = tmp[1]+ychsz*0.35
            xyouts, tx,ty,normal=1, msg, charsize=label_size, color=color
        endforeach
    endforeach
    
;    pid = where(plot_vars eq b_spec_var)
;    tpos = poss[*,pid]
    kaw_color = sgcolor('red')
    timebar, kaw_times, linestyle=1, color=kaw_color
;    xrange = time_range
;    yrange = [0,1]
;    plot, xrange, yrange, $
;        xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
;        position=tpos, nodata=1, noerase=1
;    foreach kaw_time, kaw_times, id do begin
;        tmp = convert_coord(kaw_time,yrange[1], data=1, to_normal=1)
;        tx = tmp[0]
;        ty = tmp[1]-ychsz*1.0
;        msg = string(id+1,format='(I0)')
;        xyouts, tx,ty,normal=1, msg, alignment=0.5, color=kaw_color
;        
;        if id eq 0 then begin
;            msg = 'KAW times'
;            tx = tx-xchsz*10
;            xyouts, tx,ty,normal=1, msg, charsize=label_size, color=kaw_color
;        endif
;    endforeach
    
    
    ; Add PS label
    pid = where(plot_vars eq b_tilt_combo_var)
    tpos = poss[*,pid]
    ps_color = sgcolor('wheat')
    xrange = time_range
    yrange = [0,1]
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
        position=tpos, nodata=1, noerase=1
    tx = mean(ps_times)
    ty = 1
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]+xchsz*2.5
    ty = tmp[1]-ychsz*1
    msg = 'PS entry'
    xyouts, tx,ty, normal=1, msg, alignment=0, charsize=label_size, color=sgcolor('goldenrod')
    

;---Add FMLat.
    pid = where(plot_vars eq keo_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        yrange = get_setting(keo_var, 'yrange')
        set_axis, keo_var, position=tpos, xrange=time_range, yrange=yrange
        foreach external_model, external_models, id do begin
            if external_model ne 't89' then continue
            fmlat_var = prefix+'fmlat_dipole_'+external_model+'_north'
            fmlats = abs(get_var_data(fmlat_var, times=times, in=time_range))
;            oplot, times, fmlats, color=model_colors[id]
        endforeach
        oplot, time_range, [0,0]+63, linestyle=1
    endif

    
    if keyword_set(test) then stop
    sgclose

    
    return, plot_file


end

test = 0
print, fig_long_conjunction_v01(event_info=event_info, test=test)
end