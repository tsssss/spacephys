;+
; Plot quantities showing the observed waves are Alfven waves.
;-

function fig_2013_0501_0850_alfven_wave_v01, plot_file, event_info=event_info

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
            'fig_2013_0501_0850_alfven_wave_v01.pdf'])
    endif
    if keyword_set(test) then plot_file = 0


;---plot_vars options
    label_size = 0.8


    edot0_fac_var = prefix+'edot0_fac'
    stplot_split, edot0_fac_var
    e_vars = edot0_fac_var+'_comp'+['1','2','3']
    labels = get_setting(edot0_fac_var, 'coord_labels')
    foreach var, e_vars, var_id do begin
        add_setting, var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'E!D'+labels[var_id]+'!N', $
            'unit', 'mV/m' )
        options, var, 'yrange', [-200,200]
        options, var, 'yticks', 2
        options, var, 'yminor', 5
    endforeach

    pf_spec_var = prefix+'pf_fac_mor_spec_1'
    
    pf_var = prefix+'pf_fac_map'
    get_data, pf_var, times, pf_fac, limits=lim
    
    pf_para_var = prefix+'pf_fac_para'
    store_data, pf_para_var, times, [[snorm(pf_fac)],[pf_fac[*,0]]]
    add_setting, pf_para_var, smart=1, dictionary($
        'display_type', 'stack', $
        'labels', strarr(2)+' ', $
        'colors', sgcolor(['black','red']), $
        'the_labels', ['|S|','S!D||!N'], $
        'the_ytitle', '(mW/m!U2!N)', $
        'constant', 1, $
        'yrange', [-20,50], $
        'ytickv', [0,1]*50, $
        'yticks', 1, $
        'yminor', 5, $
        'ystyle', 1, $
        'ylog', 0, $
        'ytitle', ' ' )
    pf_para_var_above = pf_para_var+'_above'
    copy_data, pf_para_var, pf_para_var_above
    var = pf_para_var_above
    options, var, 'yrange', [50,400]
    options, var, 'ystyle', 1
    options, var, 'ytickv', [150,300]
    options, var, 'yticks', 1
    options, var, 'yminor', 3
    pf_para_vars = [pf_para_var_above, pf_para_var]

;---Alfven speed.
    va_var = prefix+'va'
    b_var = prefix+'b_gsm'
    n_var = prefix+'density_hope'
    o_var = prefix+'o_density'
    p_var = prefix+'p_density'
    bmag = snorm(get_var_data(b_var, times=times))
    num_dens = get_var_data(n_var, at=times)
    o_dens = get_var_data(o_var, at=times)
    p_dens = get_var_data(p_var, at=times)
    o_ratio = o_dens/(o_dens+p_dens)
    p_ratio = 1-o_ratio
    p_mass = 1d
    o_mass = 16d
    avg_mass = p_ratio*p_mass+o_ratio*o_mass

    va0 = 1e-9/sqrt(1e6*!dpi*4e-7*1.67e-27)*1e-3    ; km/s, B in nT, n in cc, m in m_p.
    va = va0*bmag/sqrt(num_dens*avg_mass)
    store_data, va_var, times, va
    add_setting, va_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'v!DA!N', $
        'unit', 'km/s' )


;---keo var.
    keo_var = 'thg_asf_keo'
    mlt_image_rect_var = themis_read_asf_mlt_image_rect(time_range, get_name=1)
    fmlt_var = prefix+'fmlt_'+internal_model+'_north'
    dmlt = 0.1
    mlat_range = [59,67]

    get_data, mlt_image_rect_var, times, data, limits=lim
    ntime = n_elements(times)
    model_index = where(all_models eq external_model)
    fmlts = (get_var_data(fmlt_var, at=times))[*,model_index]

    mlt_bins = lim.mlt_bins
    mlat_bins = lim.mlat_bins
    mlat_index = lazy_where(mlat_bins, '[]', mlat_range, count=nybin)
    if nybin eq 0 then begin
        errmsg = 'Invalid MLat range ...'
        return, retval
    endif
    ybins = mlat_bins[mlat_index]

    hehe = fltarr(ntime,nybin)
    foreach time, times, time_id do begin
        the_mlt_range = fmlts[time_id]+[-1,1]*dmlt
        mlt_index = lazy_where(mlt_bins, '[]', the_mlt_range, count=nxbin)
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


;---ASI around footpoint.
    fpt_asi_var = 'thg_asf_around_fpt'
    fmlat_var = prefix+'fmlat_'+internal_model+'_north'
    fmlats = (get_var_data(fmlat_var, at=times))[*,model_index]
    haha = fltarr(ntime)
    dmlat = 0.1d
    foreach time, times, time_id do begin
        the_mlat_range = fmlats[time_id]+[-1,1]*dmlat
        mlat_index = lazy_where(ybins, '[]', the_mlat_range, count=nybin)
        haha[time_id] = mean(hehe[time_id,mlat_index],nan=1)
    endforeach
    store_data, fpt_asi_var, times, haha
    add_setting, fpt_asi_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'ASI Count', $
        'unit', '#', $
        'yrange', [1e2,5e4], $
        'ylog', 1 )

    
;---E/B ratio.
    e_fac_var = edot0_fac_var
    e_index = 2
    e_var = prefix+'edot0_fac_out'
    e_var = stplot_index(e_fac_var, e_index, output=e_var)

    b_fac_var = prefix+'b1_fac'
    b_index = 1
    b_var = prefix+'b1_fac_west'
    b_var = stplot_index(b_fac_var, b_index, output=b_var)

    ebr_var = prefix+'ebr'
    pflux_setting = event_info['pflux_setting']
    scale_info = pflux_setting['scale_info']
    ebr_var = stplot_calc_ebratio(time_range, e_var=e_var, b_var=b_var, output=ebr_var, scale_info=scale_info)
    
    ebr_spec_var = prefix+'ebr_spec'
    e_spec = get_var_data(e_var+'_mor', fs, times=times, limits=lim)
    b_spec = get_var_data(b_var+'_mor')
    ebr_spec = abs(e_spec)/abs(b_spec)*1e3
    foreach tmp, fs, ii do ebr_spec[*,ii] /= va
    store_data, ebr_spec_var, times, ebr_spec, fs*1e3, limits=lim
    add_setting, ebr_spec_var, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', '[E/B]/v!DA!N', $
        'ylog', 1, $
        'zlog', 1, $
        'zrange', [0.1,10], $
        'ytitle', 'Freq!C(mHz)', $
        'color_table', 64, $
;        'color_table', 38, $
        'unit', '#' )
    store_data, e_var+'_psd', times, abs(e_spec), fs*1e3, limits=lim
    add_setting, e_var+'_psd', smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'PS', $
        'ylog', 1, $
        'zlog', 1, $
        'ytitle', 'Freq!C(mHz)', $
        'unit', 'mV/m' )
    store_data, b_var+'_psd', times, abs(b_spec), fs*1e3, limits=lim
    add_setting, b_var+'_psd', smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'PS', $
        'ylog', 1, $
        'zlog', 1, $
        'ytitle', 'Freq!C(mHz)', $
        'unit', 'nT' )
    
    eb_spec_vars = [e_var,b_var]+'_psd'
    vars = [pf_spec_var,ebr_spec_var,eb_spec_vars]
    options, vars, 'constant', [1e1,1e2,1e3]
    options, vars, 'yrange', [1,1e3]
    options, vars, 'ytickv', [1e0,1e1,1e2]
    options, vars, 'yticks', 2
    options, vars, 'yminor', 10
    
    b1_fac_var = prefix+'b1_fac'
    var = b1_fac_var
    options, var, 'yrange', [-1,1]*60
    options, var, 'ytickv', [-1,0,1]*40
    options, var, 'yticks', 2
    options, var, 'yminor', 4
    options, var, 'labels', 'dB!D'+['||',tex2str('perp')+','+['west','out']]
    
;    plot_vars = [keo_var, fpt_asi_var, e_vars[[2:1]],b1_fac_var,pf_para_vars,pf_spec_var,ebr_spec_var,[e_var,b_var]+'_psd']
    plot_vars = [e_vars[[2:1]],b1_fac_var,pf_para_vars,pf_spec_var,ebr_spec_var,[e_var,b_var]+'_psd']
    nvar = n_elements(plot_vars)
    fig_labels = letters(nvar-1)+') '+['','','dB','S','S!D||!N spec','E/B',$
        '|E|!D'+tex2str('perp')+',out','|B|!D'+tex2str('perp')+',west']
    margins = [10,4,10,1]
    ypads = fltarr(nvar-1)+0.4
    ypans = fltarr(nvar)+1
    pid = where(plot_vars eq pf_para_var, count)
    if count ne 0 then begin
        ypads[pid-1] = 0
        ypans[pid] = 0.7
        ypans[pid-1] = 0.5
    endif
    
    sgopen, plot_file, fig_size=[6,6]
    poss = sgcalcpos(nvar, margins=margins, ypad=ypads, ypans=ypans, xchsz=xchsz, ychsz=ychsz)
    
    tplot, plot_vars, trange=time_range, position=poss

    ; Add FMLat.
    rbsp_color = sgcolor('purple')
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
        oplot, times, fmlats[*,model_index], color=rbsp_color, linestyle=2
        tx = times[0]
        ty = fmlats[0,model_index]
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.35
        xyouts, tx,ty,normal=1, 'RBSP-'+strupcase(probe), charsize=label_size, alignment=0, color=rbsp_color
    endif

    
    ; Add label for pf spec.
    pid = where(plot_vars eq pf_spec_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = 'Positive is toward Earth'
        xyouts, tx,ty,msg, normal=1, charsize=label_size
    endif
    
    ; Add labels for pf para.
    pid = where(plot_vars eq pf_para_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        tpos[3] = poss[3,pid-1]
        
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = 'Normalized to 100 km altitude'
        xyouts, tx,ty,msg, normal=1, charsize=label_size
        
        get_data, pf_para_var, limits=lim
        ytitle = lim.the_ytitle
        tx = tpos[0]-xchsz*5
        ty = (tpos[1]+tpos[3])*0.5
        msg = ytitle
        xyouts, tx, ty, msg, normal=1, orientation=90, alignment=0.5
        
        labels = lim.the_labels
        colors = lim.colors
        dy = (tpos[3]-tpos[1])/(n_elements(colors)+1)
        foreach msg, labels, ii do begin
            tx = tpos[2]+xchsz*0.5
            ty = tpos[1]+dy*(ii+1)
            xyouts, tx,ty, msg, normal=1, color=colors[ii]
        endforeach
    endif
    for ii=0,nvar-2 do begin
        tpos = poss[*,ii]
        if ii ge pid then begin
            tpos = poss[*,ii+1]
        endif
        msg = fig_labels[ii]
        tx = xchsz*2
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty, msg, normal=1
    endfor
    
    if keyword_set(test) then stop
    sgclose
    return, plot_file


end

print, fig_2013_0501_0850_alfven_wave_v01(event_info=event_info)
end