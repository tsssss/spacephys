;+
; Showing RBSP conjunction with ASI.
;-

function fig_2015_0416_0800_ps_arc_v01, plot_file, event_info=event_info

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
test = 1


;    time_range = time_double(['2015-04-16/07:30','2015-04-16/08:20'])
    time_range = time_double(['2015-04-16/08:00','2015-04-16/08:25'])
    bar_times = time_double('2015-04-16/'+['08:10','08:11'])
    bar_times = !null
    dp_times = time_double('2015-04-16/'+['08:07:30','08:09:30'])
    rbsp_info = event_info.rbsp.rbspa
    model_setting = rbsp_info['model_setting']
    all_models = model_setting['models']
    internal_model = event_info.internal_model
    external_model = event_info.external_model
    model_index = where(all_models eq external_model)

    if n_elements(plot_file) eq 0 then begin
        plot_dir = event_info['plot_dir']
        plot_file = join_path([plot_dir,$
            'fig_'+event_info.id+'_ps_arc_'+version+'.pdf'])
    endif
    if keyword_set(test) then plot_file = 0


    obs_info = dictionary()
    break_time = time_double('2015-04-16/08:11:46')
    break_energy = 3200d
    obs_info['o1'] = dictionary($
        'species', 'o', $
        'spec_var', 'rbspa_o_en_spec_para', $
        'nanti_appearance', 1, $
        'npara_appearance', 1, $
;        'trs', time_double('2015-04-16/'+['08:10:53','08:14:02']), $
;        'ens', [8266d,972] )
        'trs', [time_double('2015-04-16/08:10:38'),break_time], $
        'ens', [9000d,break_energy] )
    obs_info['o2'] = dictionary($
        'species', 'o', $
        'spec_var', 'rbspa_o_en_spec_para', $
        'nanti_appearance', 1, $
        'npara_appearance', 1, $
        ;        'trs', time_double('2015-04-16/'+['08:10:53','08:14:02']), $
        ;        'ens', [8266d,972] )
        'trs', [break_time,time_double('2015-04-16/08:13:41')], $
        'ens', [break_energy,1200d] )
    
;---plot_vars options
    label_size = 0.8
    tmp = smkarthm(0,2*!dpi,50,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1
    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    
    
    ; ewogram.
    ewo_var = 'thg_asf_ewo'
    mlat_range = [62,64]
    mlt_range = [-4,-2]
    ewo_zrange = [5e2,1e4]

    
    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
    get_data, mlt_image_rect_var, times, data, limits=lim
    ntime = n_elements(times)
    mlt_bins = lim.mlt_bins
    mlt_index = where_pro(mlt_bins, '[]', mlt_range, count=nmlt_bin)
    mlt_bins = mlt_bins[mlt_index]
    index = where(mlt_bins le 0, count)
    if count ne 0 then mlt_bins[index] += 24
    mlat_bins = lim.mlat_bins
    mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nmlat_bin)
    mlat_bins = mlat_bins[mlat_index]
    ewo = total(data[*,mlt_index,mlat_index],3)/nmlat_bin>0.01
    store_data, ewo_var, times, ewo, mlt_bins
    yrange = mlt_range+24
    ystep = 1
    ytickv = make_bins(yrange,ystep,inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 4
    add_setting, ewo_var, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'ASI Count', $
        'unit', '#', $
        'ytitle', 'MLT!C(h)', $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'zrange', ewo_zrange, $
        'zlog', 1, $
        'color_table', 49 )
    
    ; keogram.
    keo_var = 'thg_asf_keo'
    mlat_range = [60,67]
    mlt_range = [-3.5,-3]
    keo_zrange = ewo_zrange

    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
    get_data, mlt_image_rect_var, times, data, limits=lim
    ntime = n_elements(times)
    mlt_bins = lim.mlt_bins
    mlt_index = where_pro(mlt_bins, '[]', mlt_range, count=nmlt_bin)
    mlt_bins = mlt_bins[mlt_index]
    index = where(mlt_bins le 0, count)
    if count ne 0 then mlt_bins[index] += 24
    mlat_bins = lim.mlat_bins
    mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nmlat_bin)
    mlat_bins = mlat_bins[mlat_index]
    keo = total(data[*,mlt_index,mlat_index],2)/nmlt_bin
    store_data, keo_var, times, keo, mlat_bins
    add_setting, keo_var, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'ASI Count', $
        'unit', '#', $
        'ytitle', 'MLat!C(deg)', $
        'zrange', keo_zrange, $
        'zlog', 1, $
        'color_table', 49 )

    
    ; keogram around sc mlt.
    dmlt = 0.1
    mlat_range = [60,68]

    foreach info, event_info.rbsp do begin
        prefix = info['prefix']
        keo_var = prefix+'keo'
        
        mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
        get_data, mlt_image_rect_var, times, data, limits=lim
        ntime = n_elements(times)
        mlt_bins = lim.mlt_bins
        mlat_bins = lim.mlat_bins
        mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nmlat_bin)
        mlat_bins = mlat_bins[mlat_index]
        data = data[*,*,mlat_index]
        
        
        fmlt_var = prefix+'fmlt_'+internal_model+'_north'
        fmlts = (get_var_data(fmlt_var, at=times))[*,model_index]
        hehe = fltarr(ntime,nmlat_bin)
        foreach time, times, time_id do begin
            the_mlt_range = fmlts[time_id]+[-1,1]*dmlt
            mlt_index = where_pro(mlt_bins, '[]', the_mlt_range, count=nmlt_bin)
            hehe[time_id,*] = total(data[time_id,mlt_index,*],2)/nmlt_bin
        endforeach
        yrange = mlat_range
        ytitle = 'MLat (deg)'

        ystep = 5
        ytickv = make_bins(yrange, ystep, inner=1)
        yticks = n_elements(ytickv)-1
        yminor = ystep
        store_data, keo_var, times, hehe, mlat_bins
        add_setting, keo_var, smart=1, dictionary($
            'display_type', 'spec', $
            'short_name', 'ASI Count', $
            'unit', '#', $
            'ytitle', 'MLat!C(deg)', $
            'zrange', keo_zrange, $
            'zlog', 1, $
            'yrange', yrange, $
            'ytickv', ytickv, $
            'yticks', yticks, $
            'yminor', yminor, $
            'color_table', 49 )
    endforeach
    
    vars = ['thg_asf_'+['ewo','keo'],'rbsp'+['a','b']+'_keo']
    zrange = [2e2,1e4]
    log_zrange = alog10(zrange)
    log_ztickv = make_bins(log_zrange,1,inner=1)
    ztickv = 10d^log_ztickv
    zticks = n_elements(ztickv)-1
    zminor = 9
    ztickn = '10!U'+string(log_ztickv,format='(I0)')
    
    options, vars, 'zrange', zrange
    options, vars, 'ztickv', ztickv
    options, vars, 'zticks', zticks
    options, vars, 'zminor', zminor
    options, vars, 'ztickname', ztickn


    ; spec.
    foreach info, event_info.rbsp do begin
        prefix = info['prefix']
        probe = info['probe']
    
    
        ; e- spec.
        e_spec_var = prefix+'e_en_spec'+['','_'+['anti','perp','para']]
        ct_electron = 65
        zrange = [1e5,1e10]
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,2,inner=1)
        ztickv = 10.^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        
        var = e_spec_var
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'zminor', zminor
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

        foreach tvar, var do begin
            get_data, tvar, times, data, val, limits=lim
            index = where(finite(data,nan=1) or data eq 0, count)
            if count ne 0 then begin
                data[index] = 0.001
                store_data, tvar, times, data, val, limits=lim
            endif
        endforeach


        ; H+ spec.
        p_spec_var = prefix+'p_en_spec'+['','_'+['anti','perp','para']]
        ct_proton = 63
        zrange_proton = [1e4,1e6]

        zrange = zrange_proton
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,1,inner=1)
        ztickv = 10.^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9

        var = p_spec_var
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'zminor', zminor
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

        foreach tvar, var do begin
            get_data, tvar, times, data, val, limits=lim
            index = where(finite(data,nan=1) or data eq 0, count)
            if count ne 0 then begin
                data[index] = 0.001
                store_data, tvar, times, data, val, limits=lim
            endif
        endforeach


        ; O+ spec.
        o_spec_var = prefix+'o_en_spec'+['','_'+['anti','perp','para']]
        ct_oxygen = 64
        zrange_proton = [1e4,1e6]

        zrange = zrange_proton
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,1,inner=1)
        ztickv = 10.^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        
        var = o_spec_var
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'zminor', zminor
        options, var, 'zcharsize', label_size
        options, var, 'color_table', ct_oxygen

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

        foreach tvar, var do begin
            get_data, tvar, times, data, val, limits=lim
            index = where(finite(data,nan=1) or data eq 0, count)
            if count ne 0 then begin
                data[index] = 0.001
                store_data, tvar, times, data, val, limits=lim
            endif
        endforeach
        
        ; density.
        var = prefix+'density_hope'
        options, var, 'yrange', [0.1,1]
        ;options, var, 'labels', 'e- Density!C  >200 eV'
        options, var, 'labels', ' '
        options, var, 'ytitle', '(cm!U-3!N)'
        options, var, 'ystyle', 1

        ; B var.
        var = prefix+'b_gsm'
        get_data, var, times, b_gsm
        b_sm = cotran(b_gsm, times, 'gsm2sm')
        b_sm_var = prefix+'b_sm'
        store_data, b_sm_var, times, b_sm
        add_setting, b_sm_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'B', $
            'unit', 'nT', $
            'coord', 'SM' )
        b_tilt = asin(b_sm[*,2]/snorm(b_sm))*constant('deg')
        b_tilt_var = prefix+'b_tilt'
        store_data, b_tilt_var, times, b_tilt
        add_setting, b_tilt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'arcsin(B!Dz!N/|B|)', $
            'unit', 'deg', $
            'yrange', [0,50], $
            'ytickv', [20,40], $
            'yticks', 1, $
            'yminor', 4, $
            'ystyle', 1 )
        
        ; B spec.
        b_fac_var = prefix+'b1_fac'
        b_vars = stplot_split(b_fac_var)
        the_b_var = b_vars[1]
        b_mor_var = stplot_mor_new(the_b_var, scale_info=scale_info)
        var = b_mor_var
        zrange = [1,1e6]
        log_ztickv = make_bins(alog10(zrange),1,inner=1)
        ztickv = 10d^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        ztickn = '10!U'+string(log_ztickv,format='(I0)')
        index = where(log_ztickv eq 0, count)
        if count ne 0 then ztickn[index] = '1'
        index = where(log_ztickv eq 1, count)
        if count ne 0 then ztickn[index] = '10'
        ztickn[0:*:2] = ' '
        
        log_ytickv = [0d,1,2,3]
        ytickv = 10d^log_ytickv
        yticks = n_elements(ytickv)-1
        yminor = 9
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        index = where(log_ytickv eq 0, count)
        if count ne 0 then ytickn[index] = '1'
        index = where(log_ytickv eq 1, count)
        if count ne 0 then ytickn[index] = '10'
        options, var, 'display_type', 'spec'
        options, var, 'ztitle', 'B!D'+fac_labels[1]+'!N (nT)!U2!N'
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'zminor', zminor
        options, var, 'ztickname', ztickn
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        options, var, 'yrange', lim.yrange*1e3
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'ytickname', ytickn
        options, var, 'yminor', yminor
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', ytickv
        options, var, 'zcharsize', label_size
    endforeach
    
    tilt_var = 'b_tilt_combo'
    rbsp_colors = list()
    rbsp_probes = list()
    foreach info, event_info.rbsp do begin
        rbsp_colors.add, info['sc_color']
        rbsp_probes.add, info['probe']
    endforeach
    rbsp_colors = rbsp_colors.toarray()
    rbsp_probes = rbsp_probes.toarray()
    index = sort(rbsp_probes)
    rbsp_probes = rbsp_probes[index]
    rbsp_colors = rbsp_colors[index]
    tmp = stplot_merge('rbsp'+rbsp_probes+'_b_tilt', output=tilt_var)
    add_setting, tilt_var, smart=1, dictionary($
        'display_type', 'stack', $
        'ytitle', '(deg)', $
        'yrange', [5,55], $
        'ytickv', [20,40], $
        'yticks', 1, $
        'yminor', 4, $
        'colors', rbsp_colors, $
        'labels', strupcase('rbsp-'+rbsp_probes) )
    ; calculate correlation.
;    ; time to be compared for RBSP-A.
;    the_tr = time_double('2015-04-16/'+['08:05','08:10'])
;    ; time shift to be applied to RBSP-B.
;    dt = 1
;    time_shifts = make_bins(-120+[-1,1]*30, dt)
;    ntime_shift = n_elements(time_shifts)
;    corr = fltarr(ntime_shift)
;    a_tilt = get_var_data('rbspa_b_tilt', in=the_tr, times=times)
;    foreach time_shift, time_shifts, ts_id do begin
;        b_tilt = get_var_data('rbspb_b_tilt', at=times)
;        corr[ts_id] = c_correlate(a_tilt, b_tilt, 0)
;    endforeach
;    max_corr = max(corr, index)
;    best_time_shift = time_shifts[index]
;    print, 'Best time shift (min): ', best_time_shift/60
;    mlt_a = get_var_data(rbsp_read_mlt(time_range, probe='a'), at=the_tr[0])
;    mlt_b = get_var_data(rbsp_read_mlt(time_range, probe='b'), at=the_tr[0]+abs(best_time_shift))
;    omega_b = abs(mlt_a-mlt_b)*15/(best_time_shift/60)
;    print, 'omega (deg/min): ', omega_b
;    stop
    
    
    edot0_fac_var = prefix+'edot0_fac_spinfit_phasef'
    e_vars = stplot_split(edot0_fac_var)
    e_var = e_vars[2]
    var = e_var
    add_setting, var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'E', $
        'unit', 'mV/m' )
    options, var, 'yrange', [-1,1]*120
    options, var, 'ytickv', [-1,1]*100
    options, var, 'yticks', 2
    options, var, 'yminor', 5
    options, var, 'labels', 'E!D'+tex2str('perp')+',out'
        
    
    
    ; Poynting flux.
    pf_var = 'rbspb_pfdot0_fac_spinfit_phasef_map'
    var = pf_var
    options, var, 'yrange', [-400,100]
    options, var, 'ytickv', [-400,-200,0]
    options, var, 'yticks', 2
    options, var, 'yminor', 4
    options, var, 'constant', [0]

    pf_var = 'rbspa_pfdot0_fac_spinfit_phasef_map'
    var = pf_var
    options, var, 'yrange', [-400,100]
    options, var, 'ytickv', [-400,-200,0]
    options, var, 'yticks', 2
    options, var, 'yminor', 4
    options, var, 'constant', [0]
    
    pf_var = 'rbspb_pfdot0_fac_survey_map'
    var = pf_var
    options, var, 'yrange', [-400,100]
    options, var, 'ytickv', [-400,-200,0]
    options, var, 'yticks', 2
    options, var, 'yminor', 4
    options, var, 'constant', [0]

;---Set plot_vars.
;    rbsp_vars = [e_var,pf_var,prefix+['density_hope','e_en_spec','o_en_spec_para']]
;    asi_vars = [ewo_var,tilt_var]
;    plot_vars = [asi_vars,rbsp_vars]
;    nvar = n_elements(plot_vars)
;    fig_labels = letters(nvar)+') '+['Aurora','B tilt','E','S','N!De!N & ','e-','O+']
;    ypans = [1.5,1,1,1,0.8,1,1]
    big_ypad = 1
    plot_info = orderedhash()
    plot_info[ewo_var] = dictionary($
        'fig_label', 'Aurora!CEwogram', $
        'ypan', 1.5, $
        'ypad', big_ypad )
    plot_info['rbspa_o_en_spec_para'] = dictionary($
        'fig_label', 'S-hem O+' )
    plot_info['rbspa_o_en_spec_anti'] = dictionary($
        'fig_label', 'N-hem O+' )    
    plot_info['rbspa_b1_fac_comp2_mor'] = dictionary($
        'fig_label', 'Wave!CRBSP-A', $
        'ypan', 1.2, $
        'ypad', big_ypad )
    plot_info[tilt_var] = dictionary($
        'fig_label', 'B tilt', $
        'ypad', big_ypad )
    plot_info['rbspb_b1_fac_comp2_mor'] = dictionary($
        'fig_label', 'Wave!CRBSP-B', $
        'ypan', 1.2 )
    plot_info['rbspb_pfdot0_fac_spinfit_phasef_map'] = dictionary($
        'fig_label', 'S @100 km' )
    plot_info['rbspb_pfdot0_fac_survey_map'] = dictionary($
        'fig_label', 'S @100 km' )
    plot_info['rbspa_pfdot0_fac_spinfit_phasef_map'] = dictionary($
        'fig_label', 'S @100 km' )    
;    plot_info['rbspb_keo'] = dictionary($
;        'fig_label', 'Aurora!CKeogram', $
;        'ypan', 1 )
;    plot_info['rbspa_p_en_spec_anti'] = dictionary($
;        'fig_label', 'H+ anti' )
;    plot_info['rbspa_p_en_spec_para'] = dictionary($
;        'fig_label', 'H+ para' )
    plot_vars = (plot_info.keys()).toarray()
    nvar = n_elements(plot_vars)
    fig_labels = strarr(nvar)
    ypans = fltarr(nvar)
    ypads = fltarr(nvar)
    foreach key, plot_info.keys(), var_id do begin
        info = plot_info[key]
        fig_labels[var_id] = info['fig_label']
        if ~info.haskey('ypan') then info['ypan'] = 1.
        if ~info.haskey('ypad') then info['ypad'] = 0.4
        ypans[var_id] = info['ypan']
        ypads[var_id] = info['ypad']
    endforeach
    ypads = ypads[0:nvar-2]
    fig_labels = letters(nvar)+') '+fig_labels
    
   
;---Plot.
    fig_size = [6,7]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    margins = [12,4,10,1]
    poss = sgcalcpos(nvar, margins=margins, ypans=ypans, ypad=ypads)
    uniform_ticklen = -ychsz*0.2*fig_size[0]
    
    
    ; ticklen.
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, plot_vars[pid], 'xticklen', xticklen
        options, plot_vars[pid], 'yticklen', yticklen
        var = plot_vars[pid]
        display_type = get_setting(var, 'display_type')
        if display_type eq 'spec' then begin
            zticklen = uniform_ticklen/xchsz*1/fig_size[0]
            options, var, 'zticklen', zticklen
        endif
    endfor
    symsize = 0.4
    
    
    tplot, plot_vars, position=poss, trange=time_range
    label_yshift = -ychsz*0.7
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]+label_yshift
        xyouts, tx,ty, normal=1, fig_labels[pid]
    endfor
    if n_elements(bar_times) ne 0 then timebar, bar_times, color=sgcolor('red'), linestyle=1
    
    ; Add bars.
    bar_info = dictionary()
    bar_info['rbspa'] = dictionary($
        'time_bar', dp_times[0], $
        'panel_range', [1,4] )
    bar_info['rbspb'] = dictionary($
        'time_bar', dp_times[1], $
        'panel_range', [4,6] )
    foreach key, bar_info.keys() do begin
        info = bar_info[key]
        sc_info = (event_info.rbsp)[key]
        color = sc_info['sc_color']
        panel_range = info['panel_range']
        tpos = poss[*,panel_range[0]]
        tpos[1] = poss[1,panel_range[1]]
        xrange = time_range
        yrange = [0,1]
        plot, xrange, yrange, $
            position=tpos, xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
            nodata=1, noerase=1
        foreach tx, info['time_bar'] do begin
            plots, tx+[0,0], yrange, color=color, linestyle=3
        endforeach
    endforeach
    

    ; Add FMLT.
    ewo_var = 'thg_asf_ewo'
    pid = where(plot_vars eq ewo_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        get_data, ewo_var, limits=lim
        xrange = time_range
        yrange = lim.yrange
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        
        dp_mlts = []
        foreach key, sort_uniq(event_info.rbsp.keys()), probe_id do begin
            info = (event_info.rbsp)[key]
            prefix = info['prefix']
            probe = info['probe']
            sc_color = info['sc_color']
            fmlt_var = prefix+'fmlt_'+internal_model+'_north'
            get_data, fmlt_var, times, fmlts
            fmlts = fmlts[*,model_index]+24
            oplot, times, fmlts, color=sc_color, linestyle=2
            ; Add sc.
            tx = time_range[0]
            ty = interpol(fmlts,times,tx)
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]+xchsz*0.5
            ty = tmp[1]+ychsz*0.3
            sc_name = info['sc_name']
            msg = strupcase(sc_name+'-'+probe)
            xyouts, tx,ty,normal=1, msg, charsize=label_size, alignment=0, color=sc_color
            ; Add dp time.
            tx = dp_times[probe_id]
            ty = interpol(fmlts,times,tx)
            plots, tx,ty, psym=8, symsize=0.5, color=sc_color
            dp_mlts = [dp_mlts,ty]
        endforeach
        
        plots, dp_times,dp_mlts, linestyle=1
        dp_omega = total(dp_mlts*[-1,1])/total(dp_times*[-1,1])*60*15
        tx = mean(dp_times)
        ty = mean(dp_mlts)
        tmp = convert_coord(tx, ty, data=1, to_normal=1)
        tx = tmp[0]-xchsz*1
        ty = tmp[1]-ychsz*0.35
        msg = tex2str('omega')+' = '+string(abs(dp_omega),format='(F3.1)')+' deg/min'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size
    endif
    

    ; Add FMLat.
    foreach key, sort_uniq(event_info.rbsp.keys()), probe_id do begin
        info = (event_info.rbsp)[key]
        prefix = info['prefix']
        probe = info['probe']
        sc_color = info['sc_color']
        
        keo_var = prefix+'keo'
        pid = where(plot_vars eq keo_var, count)
        if count eq 0 then continue
        
        tpos = poss[*,pid]
        get_data, keo_var, limits=lim
        xrange = time_range
        yrange = lim.yrange
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        
        
        fmlat_var = prefix+'fmlat_'+internal_model+'_north'
        get_data, fmlt_var, times, fmlats
        fmlats = fmlts[*,model_index]
        oplot, times, fmlats, color=sc_color, linestyle=2
        ; Add sc.
        tx = time_range[0]
        ty = interpol(fmlts,times,tx)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.3
        sc_name = info['sc_name']
        msg = strupcase(sc_name+'-'+probe)
        xyouts, tx,ty,normal=1, msg, charsize=label_size, alignment=0, color=sc_color
    endforeach
    
    
    rbsp_color = sgcolor('purple')
    pid = where(plot_vars eq keo_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        fmlat_var = prefix+'fmlat_'+internal_model+'_north'
        get_data, fmlat_var, times, fmlats
                
        get_data, keo_var, limits=lim
        xrange = time_range
        yrange = lim.yrange
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        oplot, times, fmlats[*,model_index], color=rbsp_color, linestyle=2
        tx = time_range[0]
        ty = interpol(fmlats[*,model_index],times,tx)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.3
        xyouts, tx,ty,normal=1, 'RBSP-'+strupcase(probe), charsize=label_size, alignment=0, color=rbsp_color
    endif
    
    

    ; Add electron temperature.
    var = e_spec_var
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
    
    ; Add oxygen labels.
    var = 'rbspa_o_en_spec_para'
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.3
        msg = 'Away from Earth in S-hem, PA [0,45] deg'
        xyouts, tx,ty,msg, normal=1, charsize=label_size
    endif
    
    var = 'rbspa_o_en_spec_anti'
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos

        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.3
        msg = 'Toward Earth in S-hem, PA [135,180] deg'
        xyouts, tx,ty,msg, normal=1, charsize=label_size
    endif
    
    
    
    
    ; Add labels for Pflux.
    var = 'rbspb_pfdot0_fac_spinfit_phasef_map'
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        
        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=0, $
            nodata=1, noerase=1, position=tpos
        
        filter = (rbsp_info['pflux_setting']).filter
        

        tx = tpos[2]-xchsz*0.5
        ty = tpos[1]+ychsz*0.3
        msg = 'Normalized to 100 km altitude'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size
        ty = tpos[1]+ychsz*1.3
        msg = 'Filtered in '+string(1d3/filter[1],format='(F3.1)')+'mHz-'+string(1d/filter[0],format='(I0)')+'Hz'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size
        
        tmp = convert_coord(xrange[0],0, data=1, to_normal=1)
        tx = tpos[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.3
        xyouts, tx,ty,normal=1, alignment=0, 'Away from Earth in S-hem', charsize=label_size;, color=sgcolor('red')
        ty = tmp[1]-ychsz*0.9
        xyouts, tx,ty,normal=1, alignment=0, 'Toward Earth in S-hem', charsize=label_size;, color=sgcolor('red')
    endif

    ; Add B tilt.
    var = prefix+'density_hope'
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        the_var = prefix+'b_tilt'
        the_color = sgcolor('gray')

        xrange = time_range
        yrange = get_setting(the_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=0, $
            nodata=1, noerase=1, position=tpos

        get_data, the_var, times, data
        oplot, times, data, color=the_color
        ytickv = get_setting(the_var, 'ytickv')
        yticks = get_setting(the_var, 'yticks')
        yminor = get_setting(the_var, 'yminor')
        yticklen = get_setting(var, 'yticklen')
        ytitle = get_setting(the_var, 'ytitle')
        axis, yaxis=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, color=the_color
        tx = tpos[0]-xchsz*5
        ty = tpos[3]+label_yshift
        msg = tex2str('theta')+'!DB!N'
        xyouts, tx,ty,normal=1, msg, color=the_color
        tx = tpos[2]-xchsz*0.5
        ty = tpos[1]+ychsz*0.3
        msg = tex2str('theta')+'!DB!N = arcsin(B!Dz!N/|B|)'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size, color=the_color
        ty = tpos[1]+ychsz*1.1
        msg = 'e- density >200 eV'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size;, color=the_color
    endif

    

    if keyword_set(test) then stop
    sgclose
    
    return, plot_file

end



print, fig_2015_0416_0800_ps_arc_v01(event_info=event_info)
end