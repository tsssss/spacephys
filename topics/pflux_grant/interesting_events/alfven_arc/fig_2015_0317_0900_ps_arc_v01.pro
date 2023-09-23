;+
; Showing RBSP conjunction with ASI.
;-

function fig_2015_0317_0900_ps_arc_v01, plot_file, event_info=event_info

;---Load data and settings.
    version = 'v01'
    id = '2015_0317_0900'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
test = 1

    rbsp_info = event_info.rbsp.rbspa
    probe = rbsp_info['probe']
    prefix = rbsp_info['prefix']

    time_range = time_double(['2015-03-17/08:20','2015-03-17/09:20'])
    bar_times = time_double('2015-03-17/'+['08:38:30','08:52:00','09:05:30'])
    model_setting = rbsp_info['model_setting']
    all_models = model_setting['models']
    internal_model = event_info.internal_model
    external_model = event_info.external_model
    external_model = 't04s'

    if n_elements(plot_file) eq 0 then begin
        plot_dir = event_info['plot_dir']
        plot_file = join_path([plot_dir,$
            'fig_'+event_info.id+'_ps_arc_'+version+'.pdf'])
    endif
    if keyword_set(test) then plot_file = 0

    
    
;---plot_vars options
    label_size = 0.8

    keo_var = 'thg_asf_keo'
    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
    fmlt_var = prefix+'fmlt_'+internal_model+'_north'
    dmlt = 0.1
    mlat_range = [55,68]

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
        the_mlt_range = fmlts[time_id]+[-1,1]*dmlt+0.5
        mlt_index = where_pro(mlt_bins, '[]', the_mlt_range, count=nxbin)
        hehe[time_id,*] = total(data[time_id,mlt_index,mlat_index],2)/nxbin
    endforeach
    yrange = mlat_range
    ytitle = 'MLat (deg)'
    
    if n_elements(ystep) eq 0 then ystep = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = ystep
    store_data, keo_var, times, hehe, ybins
    add_setting, keo_var, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'ASI Count', $
        'unit', '#', $
        'ytitle', 'MLat!C(deg)', $
        'zrange', [1e3,2e4], $
        'zlog', 1, $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'color_table', 49 )
    
    
    edot0_fac_var = prefix+'edot0_fac'
    stplot_split, edot0_fac_var
    e_var = edot0_fac_var+'_comp3'
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
   
    var = prefix+'e_n'
    options, var, 'yrange', [0.1,10]
    ;options, var, 'labels', 'e- Density!C  >200 eV'
    options, var, 'labels', ' '
    options, var, 'ytitle', '(cm!U-3!N)'
    options, var, 'ystyle', 9

    ; B tilt var.
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
        'yrange', [10,40], $
        'ytickv', [20,40], $
        'yticks', 1, $
        'yminor', 4, $
        'ystyle', 1 )
    

    ; e- spec.
    e_spec_var = prefix+'e_en_spec'
    ct_electron = 65
    zrange = [1e5,1e10]
    log_zrange = alog10(zrange)
    log_ztickv = make_bins(log_zrange,2,inner=1)
    ztickv = 10.^log_ztickv
    zticks = n_elements(ztickv)-1

    var = e_spec_var
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
    
    
    ; O+ spec.
    o_spec_var = prefix+'o_en_spec_para'
    ct_oxygen = 64
    zrange_proton = [1e4,1e6]

    zrange = zrange_proton
    log_zrange = alog10(zrange)
    log_ztickv = make_bins(log_zrange,1,inner=1)
    ztickv = 10.^log_ztickv
    zticks = n_elements(ztickv)-1

    var = o_spec_var
    options, var, 'zrange', zrange
    options, var, 'ztickv', ztickv
    options, var, 'zticks', zticks
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
    
    get_data, var, times, data, val, limits=lim
    index = where(finite(data,nan=1) or data eq 0, count)
    if count ne 0 then begin
        data[index] = 0.001
        store_data, var, times, data, val, limits=lim
    endif
    
    
    ; Poynting flux.
    pf_var = prefix+'pfdot0_fac_map'
    var = pf_var
    options, var, 'yrange', [-300,100]
    options, var, 'ytickv', [-200,-100,0,100]
    options, var, 'yticks', 3
    options, var, 'yminor', 2
    options, var, 'constant', [0]

   

;---Set plot_vars.
    rbsp_vars = [e_var,pf_var,prefix+['e_n','e_en_spec','o_en_spec_para']]
    asi_vars = 'thg_asf_keo'
    plot_vars = [asi_vars,rbsp_vars]
    nvar = n_elements(plot_vars)
    fig_labels = letters(nvar)+') '+['Aurora','E','S','N!De!N & '+tex2str('theta')+'!DB!N','e-','O+']
    fig_labels = letters(nvar)+') '+['Aurora','E','S','N!De!N & ','e-','O+']
    ypans = [1.2,1,1,0.8,1,1]

   
;---Plot.
    sgopen, plot_file, size=[6,5], xchsz=xchsz, ychsz=ychsz
    margins = [12,4,10,1]
    poss = sgcalcpos(nvar, margins=margins, ypans=ypans, ypad=0.2)
    
    ; ticklen.
    xticklen_chsz = -0.2
    yticklen_chsz = -0.5
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        options, plot_vars[pid], 'xticklen', xticklen
        options, plot_vars[pid], 'yticklen', yticklen
    endfor
    
    
    tplot, plot_vars, position=poss, trange=time_range
    label_yshift = -ychsz*0.5
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]+label_yshift
        xyouts, tx,ty, normal=1, fig_labels[pid]
    endfor
    timebar, bar_times, color=sgcolor('red'), linestyle=1


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
    var = o_spec_var
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
    
    ; Add labels for Pflux.
    var = pf_var
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
        ty = tpos[3]-ychsz*0.9
        msg = 'Filtered in '+string(1d3/filter[1],format='(F3.1)')+'mHz-'+string(1d/filter[0],format='(I0)')+'Hz'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size
        
        tmp = convert_coord(xrange[0],0, data=1, to_normal=1)
        tx = tpos[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.3
        xyouts, tx,ty,normal=1, alignment=0, 'Away from Earth', charsize=label_size, color=sgcolor('red')
        ty = tmp[1]-ychsz*0.9
        xyouts, tx,ty,normal=1, alignment=0, 'Toward Earth', charsize=label_size, color=sgcolor('red')
    endif

    ; Add B tilt.
    var = prefix+'e_n'
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



print, fig_2015_0317_0900_ps_arc_v01(event_info=event_info)
end