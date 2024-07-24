;+
; Showing RBSP conjunction with ASI.
;-

function fig_2013_0501_0850_ps_arc_v02, plot_file, event_info=event_info

;---Load data and settings.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_0850_load_data()
test = 0

    probe = event_info['probe']
    prefix = event_info['prefix']

    time_range = time_double(['2013-05-01/07:00','2013-05-01/09:30'])
    bar_times = time_double('2013-05-01/'+['07:38','08:44'])
    model_setting = event_info['model_setting']
    internal_model = model_setting['internal_model']
    external_model = model_setting['external_model']
    all_models = model_setting['models']

    if n_elements(plot_file) eq 0 then begin
        plot_dir = event_info['plot_dir']
        plot_file = join_path([plot_dir,$
            'fig_2013_0501_0850_ps_arc_v02.pdf'])
    endif
    if keyword_set(test) then plot_file = 0

    
    
;---plot_vars options
    label_size = 0.8

    keo_var = 'thg_asf_keo'
    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
    fmlt_var = prefix+'fmlt_'+internal_model+'_north'
    dmlt = 0.1
    mlat_range = [59,67]

    get_data, mlt_image_rect_var, times, data, limits=lim
    ntime = n_elements(times)
    model_index = where(all_models eq external_model)
    fmlts = (get_var_data(fmlt_var, at=times))[*,model_index]

    mlt_bins = lim.mlt_bins
    mlat_bins = lim.mlat_bins
    mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nybin)
    if nybin eq 0 then begin
        errmsg = 'Invalid MLat range ...'
        return, retval
    endif
    ybins = mlat_bins[mlat_index]

    hehe = fltarr(ntime,nybin)
    foreach time, times, time_id do begin
        the_mlt_range = fmlts[time_id]+[-1,1]*dmlt
        mlt_index = where_pro(mlt_bins, '[]', the_mlt_range, count=nxbin)
        hehe[time_id,*] = total(data[time_id,mlt_index,mlat_index],2)/nxbin
    endforeach
    yrange = mlat_range
    ytitle = 'MLat (deg)'
    
    if n_elements(ystep) eq 0 then ystep = 2
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = ystep
    store_data, keo_var, times, hehe, ybins
    add_setting, keo_var, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'ASI Count', $
        'unit', '#', $
        'ytitle', 'MLat!C(deg)', $
        'zrange', [1e3,1e4], $
        'zlog', 1, $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'color_table', 49 )
    
    
    edot0_fac_var = prefix+'edot0_fac'
    tmp = stplot_split(edot0_fac_var)
    e_var = edot0_fac_var+'_comp3'
    var = e_var
    add_setting, var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'E', $
        'unit', 'mV/m' )
    options, var, 'yrange', [-200,200]
    options, var, 'yticks', 2
    options, var, 'yminor', 5
    options, var, 'labels', 'E!D'+tex2str('perp')
   
    var = prefix+'density_hope'
    options, var, 'yrange', [0.1,3]
    options, var, 'labels', 'e- Density!C  >200 eV'

;    ; B tilt var.
;    var = prefix+'b_gsm'
;    get_data, var, times, b_gsm
;    b_sm = cotran(b_gsm, times, 'gsm2sm')
;    b_sm_var = prefix+'b_sm'
;    store_data, b_sm_var, times, b_sm
;    add_setting, b_sm_var, smart=1, dictionary($
;        'display_type', 'vector', $
;        'short_name', 'B', $
;        'unit', 'nT', $
;        'coord', 'SM' )
;    b_tilt = asin(b_sm[*,2]/snorm(b_sm))*constant('deg')
;    b_tilt_var = prefix+'b_tilt'
;    store_data, b_tilt_var, times, b_tilt
;    add_setting, b_tilt_var, smart=1, dictionary($
;        'display_type', 'scalar', $
;        'short_name', 'arcsin(B!Dz!N/|B|)', $
;        'unit', 'deg', $
;        'yrange', [10,40], $
;        'ytickv', [10,30], $
;        'yticks', 1, $
;        'yminor', 4, $
;        'ystyle', 1 )
    
    r_var = prefix+'r_gsm'
    options, r_var, 'mission_probe', 'rbsp'+probe
    mlat_vars = lets_read_mlat_vars(orbit_var=r_var)
    options, prefix+'mlt', 'ytitle', 'MLT (h)'
    options, prefix+'mlat', 'ytitle', 'MLat (deg)'
    options, prefix+'dis', 'ytitle', '|R| (Re)'
    
    
    external_models = ['t89','t96','t01','t04s']
    model_colors = sgcolor(['black','blue','green','purple'])
    bmod_gsm_vars = prefix+'bmod_gsm_'+external_models+'_igrf'
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
        'constant', [20,40], $
        'yrange', [10,50], $
        'ytickv', [20,40], $
        'yticks', 1, $
        'yminor', 4 )
    
    
    

    ; e- spec.
    e_en_spec_var = prefix+'e_en_spec'
    e_en_spec_var = rbsp_read_en_spec(time_range, probe=probe, species='e', update=1)
    ct_electron = 65
    zrange = [1e5,1e10]
    log_zrange = alog10(zrange)
    log_ztickv = make_bins(log_zrange,2,inner=1)
    ztickv = 10.^log_ztickv
    zticks = n_elements(ztickv)-1

    var = e_en_spec_var
    options, var, 'zrange', zrange
    options, var, 'ztickv', ztickv
    options, var, 'zticks', zticks
    options, var, 'zcharsize', label_size
    options, var, 'color_table', ct_electron

    yrange = [15,5e4]
    log_yrange = alog10(yrange)
    log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
    log_ytickv = make_bins(log_yrange,1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    yminor = 10
    foreach tx, ytickv, ii do begin
        if tx eq 1 then begin
            ytickn[ii] = '1'
        endif else if tx eq 10 then begin
            ytickn[ii] = '10'
        endif
    endforeach
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'ytickname', ytickn
    options, var, 'constant', ytickv
    options, var, 'ytitle', 'Energy!C(eV)'

    get_data, var, times, data, val, limits=lim
    index = where(finite(data,nan=1) or data eq 0, count)
    if count ne 0 then begin
        data[index] = 0.001
        store_data, var, times, data, val, limits=lim
    endif


    ; H+ spec.
    p_spec_var = prefix+'p_en_spec'
    ct_proton = 63
    zrange_proton = [1e4,1e6]

    zrange = zrange_proton
    log_zrange = alog10(zrange)
    log_ztickv = make_bins(log_zrange,1,inner=1)
    ztickv = 10.^log_ztickv
    zticks = n_elements(ztickv)-1

    var = p_spec_var
    options, var, 'zrange', zrange
    options, var, 'ztickv', ztickv
    options, var, 'zticks', zticks
    options, var, 'zcharsize', label_size
    options, var, 'color_table', ct_proton

    yrange = [15,5e4]
    log_yrange = alog10(yrange)
    log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
    log_ytickv = make_bins(log_yrange,1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    yminor = 10
    foreach tx, ytickv, ii do begin
        if tx eq 1 then begin
            ytickn[ii] = '1'
        endif else if tx eq 10 then begin
            ytickn[ii] = '10'
        endif
    endforeach
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'ytickname', ytickn
    options, var, 'constant', ytickv
    options, var, 'ytitle', 'Energy!C(eV)'

    get_data, var, times, data, val, limits=lim
    index = where(finite(data,nan=1) or data eq 0, count)
    if count ne 0 then begin
        data[index] = 0.001
        store_data, var, times, data, val, limits=lim
    endif
    
    
    ; E field spec.
    e_spec_var = rbsp_read_efield_spec(time_range, probe=probe)
    fc_vars = list()
    ;foreach species, ['e','o','he','p'] do fc_vars.add, rbsp_read_gyro_freq(time_range, probe=probe, species=species)
    foreach species, ['e','p'] do fc_vars.add, rbsp_read_gyro_freq(time_range, probe=probe, species=species)
    var = prefix+'fce_half'
    fce = get_var_data(prefix+'fce', times=times)
    options, prefix+'fce', labels='f!Dc,e'
    options, prefix+'fcp', labels='f!Dc,H'
    ;options, prefix+'fco', labels='f!Dc,O'
    ;options, prefix+'fche', labels='f!Dc,He'

    ;var = prefix+'flh'
    ;fcp = get_var_data(prefix+'fcp', times=times)
    ;store_data, var, times, fcp*43, limits={labels:'f!DLH!N'}
    ;fc_vars.add, var
    fc_vars = fc_vars.toarray()
    fc_colors = get_color(n_elements(fc_vars))
    foreach var, fc_vars, ii do options, var, 'colors', fc_colors[ii]

    e_spec_combo_var = e_spec_var+'_combo'
    store_data, e_spec_combo_var, data=[e_spec_var,fc_vars]
    options, e_spec_combo_var, 'yrange', get_setting(e_spec_var,'yrange')
    options, e_spec_combo_var, 'labels', ''
    

   

;---Set plot_vars.
    rbsp_vars = [e_var,prefix+['b_tilt','density_hope','e_spec_combo']]
    rbsp_vars = [b_tilt_combo_var,prefix+['density_hope','e_en_spec'],e_var,prefix+'e_spec_combo']
    asi_vars = 'thg_asf_keo'
    plot_vars = [asi_vars,rbsp_vars]
    nvar = n_elements(plot_vars)
    fig_labels = letters(nvar)+') '+['Aurora','B elev','N!De!N','e- spec','E','E spec']
    ypans = [1.2,1.2,0.8,1.2,1,1.8]
    
    vars = prefix+['e_spec_combo','e_en_spec']
    options, vars, zticklen=-0.5, zminor=9, zlog=1

   
;---Plot.
    margins = [12,6,8,2]
    poss = panel_pos(plot_file, margins=margins, nypan=nvar, ypans=ypans, pansize=[5d,0.8], fig_size=fig_size)
    
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    var_labels = prefix+['mlt','mlat','dis']
    
    tplot_options, 'tickinterval', 1200
    tplot_options, 'version', 3
    tplot, plot_vars, position=poss, trange=time_range, var_label=var_labels, vlab_margin=10
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty, normal=1, fig_labels[pid]
    endfor


    ; Add labels.
    yrange = get_setting(asi_vars[0], 'yrange')
    growth_times = time_double('2013-05-01/'+['07:17','08:18'])
    pid = where(plot_vars eq asi_vars[0], count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        set_axis, asi_vars[0], yrange=yrange, position=tpos
        
        color = sgcolor('red')
        timebar, bar_times, color=color, linestyle=1
        msg = 'Expansion!CPhase'
        foreach time, bar_times, tid do begin
            tx = time
            ty = yrange[1]
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]
            ty = tmp[1]+ychsz*1.2
            xyouts, tx,ty,msg, normal=1, alignment=0.5, charsize=label_size, color=color
        endforeach
        
        color = sgcolor('black')
        timebar, growth_times, color=color, linestyle=1
        msg = 'Growth!CPhase'
        foreach time, growth_times, tid do begin
            tx = time
            ty = yrange[1]
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]
            ty = tmp[1]+ychsz*1.2
            xyouts, tx,ty,msg, normal=1, alignment=0.5, charsize=label_size, color=color
        endforeach
        
        tx = growth_times[0]
        ty = 65
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*2
        ty = tmp[1]+ychsz*0.5
        alignment = 0
        msg = 'PBIs &!Cstreamer'
        xyouts, tx,ty,msg, normal=1, charsize=label_size, alignment=alignment
        plots, tx-[0.5,2]*xchsz,ty-[0,0.5]*ychsz, normal=1
        
        tx = growth_times[0]
        ty = 61
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*3
        ty = tmp[1]-ychsz*1
        alignment = 1
        msg = 'Equatorward arc'
        xyouts, tx,ty,msg, normal=1, charsize=label_size, alignment=alignment
        plots, tx-[5,4]*xchsz,ty+[0.5,1.2]*ychsz,normal=1
        
        tx = bar_times[0]
        ty = 62
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*3
        ty = tmp[1]-ychsz*0
        alignment = 0
        msg = 'Discrete arcs'
        xyouts, tx,ty,msg, normal=1, charsize=label_size, alignment=alignment
        plots, tx-[1,3]*xchsz,ty+[0.5,1.5]*ychsz,normal=1
        
    endif
    
    
    var = b_tilt_combo_var
    yrange = [0,1]
    xrange = time_range
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        set_axis, var, yrange=yrange, xrange=xrange, position=tpos

        color = sgcolor('red')
        timebar, bar_times, color=color, linestyle=1
        msg = 'Dipolarization'
        foreach time, bar_times, tid do begin
            tx = time
            ty = yrange[0]
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]
            ty = tmp[1]+ychsz*0.1
            xyouts, tx,ty,msg, normal=1, alignment=0.5, charsize=label_size, color=color
        endforeach
    endif
    
    
    var = e_spec_combo_var
    yrange = [0,1]
    xrange = time_range
    pid = where(plot_vars eq var, count)

        
    if count ne 0 then begin
        tpos = poss[*,pid]
        set_axis, var, yrange=yrange, xrange=xrange, position=tpos, ylog=0

        kaw_times = orderedhash($
            '2013-05-01/07:45', time_double('2013-05-01/'+['07:38:10','07:43:40','07:49:40','07:57:20']), $
            ;'2013-05-01/08:20', time_double('2013-05-01/'+['08:19:40']), $
            '2013-05-01/08:45', time_double('2013-05-01/'+['08:38:00','08:44:50','08:53:20']))
        color = sgcolor('red')
        msg = 'KAWs & TDS'
        foreach time_str, kaw_times.keys() do begin
            tx0 = time_double(time_str)
            ty0 = yrange[1]
            tmp = convert_coord(tx0,ty0, data=1, to_normal=1)
            tx0 = tmp[0]
            ty0 = tpos[3]+ychsz*0.5
            xyouts, tx0,ty0+ychsz*0.2,msg, normal=1, alignment=0.5, charsize=label_size, color=color

            foreach t2, kaw_times[time_str] do begin
                tx1 = t2
                tmp = convert_coord(tx1,ty0, data=1, to_normal=1)
                tx1 = tmp[0]
                ty1 = ty0-ychsz*3.5
                plots, [tx0,tx1],[ty0,ty1], normal=1, color=color
            endforeach
        endforeach
        
        ech_times = orderedhash($
            '2013-05-01/07:22', time_double('2013-05-01/'+['07:20:40','07:23:10']), $
            '2013-05-01/09:22', time_double('2013-05-01/'+['09:19:20','09:25:00']))
        color = sgcolor('red')
        msg = 'ECH'
        foreach time_str, ech_times.keys() do begin
            tx0 = time_double(time_str)
            ty0 = yrange[1]
            tmp = convert_coord(tx0,ty0, data=1, to_normal=1)
            tx0 = tmp[0]
            ty0 = tpos[3]+ychsz*0.5
            xyouts, tx0,ty0+ychsz*0.2,msg, normal=1, alignment=0.5, charsize=label_size, color=color

            foreach t2, ech_times[time_str] do begin
                tx1 = t2
                tmp = convert_coord(tx1,ty0, data=1, to_normal=1)
                tx1 = tmp[0]
                ty1 = ty0-ychsz*1.5
                plots, [tx0,tx1],[ty0,ty1], normal=1, color=color
            endforeach
        endforeach
    endif
    
    
    

    ; Add FMLat.
    rbsp_color = sgcolor('purple')
    dmlat = [-0.5,0.5]+0.2
    pid = where(plot_vars eq keo_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        fmlat_var = prefix+'fmlat_'+internal_model+'_north'
        get_data, fmlat_var, times, fmlats
        model_index = (where(all_models eq external_model))[0]
        
        get_data, keo_var, limits=lim
        xrange = time_range
        yrange = lim.yrange
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        oplot, times, fmlats[*,model_index]+dmlat[0], color=rbsp_color, linestyle=2
        oplot, times, fmlats[*,model_index]+dmlat[1], color=rbsp_color, linestyle=2
        ;oplot, time_range, [0,0]+64, color=rbsp_color, linestyle=2
        ;oplot, time_range, [0,0]+63, color=rbsp_color, linestyle=2
        tx = times[0]
        ty = fmlats[0,model_index]+dmlat[0]
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.12
        xyouts, tx,ty,normal=1, 'RBSP-'+strupcase(probe), charsize=label_size, alignment=0, color=rbsp_color
    endif

    ; Add electron temperature.
    var = e_en_spec_var
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        
        get_data, prefix+'e_temp', times, data
        oplot, times, data
    endif

    ; Add proton temperature.
    var = p_spec_var
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        
        get_data, prefix+'p_temp', times, data
        oplot, times, data
    endif
    
    
    ; Add labels for fcs.
    pid = where(plot_vars eq e_spec_combo_var, count)
    spec_var = e_spec_combo_var
    if count ne 0 then begin
        tpos = poss[*,pid]
        spec_settings = get_var_setting(spec_var)
        xrange = time_range
        yrange = spec_settings['yrange']
        ylog = 1
        plot, xrange, yrange, $
            ylog=ylog, xrange=xrange, yrange=yrange, $
            xstyle=5, ystyle=5, $
            position=tpos, noerase=1, nodata=1

        nfreq_var = n_elements(fc_vars)
        colors = get_color(nfreq_var)
        foreach var, fc_vars, var_id do begin
            color = colors[var_id]
            yys = get_var_data(var, times=xxs, settings=settings)
            tx = time_double('2013-05-01/07:01')
            ty = interpol(yys, xxs, tx)
            ;if product(ty-yrange) ge 0 then continue

            msg = settings.labels
            color = settings.colors

            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tpos[0]-xchsz*1.5
            ty = tmp[1]-ychsz*0.25
            xyouts, tx,ty,normal=1, msg, charsize=label_size, color=color
        endforeach
    endif
    
    

    if keyword_set(test) then stop
    sgclose
    
    return, plot_file

end



print, fig_2013_0501_0850_ps_arc_v02(event_info=event_info)
end