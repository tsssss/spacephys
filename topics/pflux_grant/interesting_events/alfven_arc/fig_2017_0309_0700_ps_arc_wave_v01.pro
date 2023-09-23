;+
; Show conjunction between Alfven wave and arc.
;-

function fig_2017_0309_0700_ps_arc_v01, plot_file, event_info=event_info

;---Load data and settings.
    version = 'v01'
    id = '2017_0309_0700'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
test = 1

    time_range = time_double(['2017-03-09/07:00','2017-03-09/09:00'])
    bar_times = time_double(['2017-03-09/'+['07:25','08:01','08:20']])
    plot_dir = event_info['plot_dir']




    foreach sc_info, event_info.themis do begin
    ;---Test B0, S, and FAC.
        prefix = sc_info['prefix']
        probe = sc_info['probe']
        sc_name = sc_info['sc_name']
        sc_color = sc_info['sc_color']


        if n_elements(plot_file) eq 0 then begin
            plot_file = join_path([plot_dir,$
                'fig_'+event_info.id+'_ps_arc_'+sc_name+probe+'_'+version+'.pdf'])
        endif
        if keyword_set(test) then plot_file = 0


        b_gsm_var = prefix+'b_gsm'
        e_dsl_var = prefix+'e_themis_dsl'
        r_gsm_var = prefix+'r_gsm'
        model_settings = sc_info['model_setting']
        external_model = 't89'
        internal_model = 'dipole'
        models = model_settings['models']
        bmod_gsm_var = prefix+'bmod_gsm_'+external_model+'_'+internal_model
        
        
        ; Calculate B0.
        b0_gsm_var = prefix+'b0_gsm'
        window = 120.
        b_gsm = get_var_data(b_gsm_var, times=times)
        bmod_gsm = get_var_data(bmod_gsm_var, at=times)
        b1_gsm = b_gsm-bmod_gsm
        time_step = sdatarate(times)
        width = window/time_step
        ndim = 3
        for ii=0, ndim-1 do begin
            b1_gsm[*,ii] -= smooth(b1_gsm[*,ii], width, nan=1, edge_mirror=1)
        endfor
        b0_gsm = b_gsm-b1_gsm
        store_data, b0_gsm_var, times, b0_gsm
        add_setting, b0_gsm_var, id='bfield', dictionary('coord', 'GSM')
        
        bmag_var = prefix+'bmag'
        bmag = snorm(b_gsm)
        store_data, bmag_var, times, bmag
        add_setting, bmag_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', '|B|', $
            'unit', 'nT' )

        
        ; B1.
        b1_gsm_var = prefix+'b1_gsm'
        b1_gsm = b_gsm-b0_gsm
        store_data, b1_gsm_var, times, b1_gsm
        add_setting, b1_gsm_var, id='bfield', dictionary('coord', 'GSM')
        
        ; Calculate B tilt.
        b_tilt_var = prefix+'b_tilt'
        b_sm = cotran_pro(b0_gsm, times, 'gsm2sm')
        b_tilt = asin(b_sm[*,2]/snorm(b_sm))*constant('deg')
        store_data, b_tilt_var, times, b_tilt
        add_setting, b_tilt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'B tilt', $
            'unit', 'deg' )
        
        ; Calculate Edot0.
        edot0_dsl_var = prefix+'edot0_themis_dsl'
        edot0_dsl= get_var_data(e_dsl_var, times=times)
        b0_gsm = get_var_data(b0_gsm_var, at=times)
        b0_dsl = cotran_pro(b0_gsm, times, 'gsm2themis_dsl', probe=probe)
        edot0_dsl[*,2] = -(edot0_dsl[*,0]*b0_dsl[*,0]+edot0_dsl[*,1]*b0_dsl[*,1])/b0_dsl[*,2]
        store_data, edot0_dsl_var, times, edot0_dsl
        add_setting, edot0_dsl_var, id='efield', dictionary('coord', 'THEMIS_DSL')
        
        edot0_gsm_var = prefix+'edot0_gsm'
        edot0_gsm = cotran_pro(edot0_dsl, times, 'themis_dsl2gsm', probe=probe)
        store_data, edot0_gsm_var, times, edot0_gsm
        add_setting, edot0_gsm_var, id='efield', dictionary('coord', 'GSM')
        
        edot0_angle_var = prefix+'edot0_angle'
        edot0_angle = asin(b0_dsl[*,2]/snorm(b0_dsl))*constant('deg')
        store_data, edot0_angle_var, times, edot0_angle
        add_setting, edot0_angle_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'edot0 angle', $
            'unit', 'deg', $
            'constant', [20,25] )
        
            
        ; Velocity.
        u_gsm_var = prefix+'u_gsm'
        uts = make_bins(minmax(times),3)
        u_gsm = get_var_data(prefix+'p_vbulk_gsm', at=uts)
        store_data, u_gsm_var, uts, u_gsm
        add_setting, u_gsm_var, id='velocity', dictionary('coord', 'GSM')
        
        
        ; Convert vars to FAC.
        define_fac, b0_gsm_var, r_gsm_var
        vars = prefix+['b','b1','edot0','u']
        foreach var, vars do begin
            gsm_var = var+'_gsm'
            fac_var = var+'_fac'
            to_fac, gsm_var, to=fac_var
            add_setting, fac_var, smart=1, dictionary('coord', 'FAC', $
                'coord_labels', fac_labels)
        endforeach
        
        
        ; Calculate pflux dot0.
        fac_labels = ['||',tex2str('perp')+','+['west','out']]
        pf_fac_var = prefix+'pfdot0_fac'
        edot0_fac_var = prefix+'edot0_fac'
        b_fac_var = prefix+'b1_fac'
        scale_info = dictionary($
            's0', time_step, $
            's1', 1200, $
            'dj', 1d/8, $
            'ns', 0d )
        filter = [time_step,300]
        stplot_calc_pflux_mor, edot0_fac_var, b_fac_var, pf_fac_var, $
            scaleinfo=scale_info, filter=filter
        add_setting, pf_fac_var, id='pflux', dictionary('coord', 'FAC', $
            'coord_labels', fac_labels)

        ; Map to 100 km.
        pf_fac_map_var = prefix+'pfdot0_fac_map'
        bf_var = prefix+'bf_gsm_'+external_model+'_'+internal_model+'_south'
        b0_var = prefix+'b0_gsm'
        b0_gsm = get_var_data(b0_var, times=times)
        bf_gsm = get_var_data(bf_var, at=times)
        cmap = snorm(bf_gsm)/snorm(b0_gsm)
        pf_fac = get_var_data(pf_fac_var)
        for ii=0,ndim-1 do pf_fac[*,ii] *= cmap
        store_data, pf_fac_map_var, times, pf_fac
        add_setting, pf_fac_map_var, id='pflux', dictionary('coord', 'FAC', 'coord_labels', fac_labels)
        
        
        ; thermal pressure.
        pressure_var = prefix+'p_thermal'
        e_n = get_var_data(prefix+'e_n', times=times)
        e_t = get_var_data(prefix+'e_tavg', at=times)
        p_t = get_var_data(prefix+'p_tavg', at=times)
        p_thermal = e_n*(e_t+p_t)*constant('kb')*1e6*11604.525*1e9
        store_data, pressure_var, times, p_thermal
        add_setting, pressure_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'P!Dthermal!N', $
            'unit', 'nPa', $
            'ylog', 1 )
        
        ; magnetic pressure.
        bp_var = prefix+'p_b'
        b_gsm = get_var_data(prefix+'b_gsm', at=times)
        bmag = snorm(b_gsm) 
        p_magnetic = (bmag*1e-9)^2/(2*constant('mu0'))*1e9
        store_data, bp_var, times, p_magnetic
        add_setting, bp_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'P!DB!N', $
            'unit', 'nPa', $
            'ylog', 1 )
            
        ; total pressure.
        p_var = prefix+'p_total'
        p_total = p_magnetic+p_thermal
        store_data, p_var, times, [[p_total],[p_magnetic],[p_thermal]]
        add_setting, p_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'P', $
            'unit', 'nPa', $
            'ylog', 1, $
            'coord', '', $
            'coord_labels', ['total','B','th'] )
        
        
        ; plasma beta.
        beta = p_thermal/p_magnetic
        beta_var = prefix+'beta'
        store_data, beta_var, times, beta
        add_setting, beta_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'beta', $
            'ylog', 1, $
            'unit', '#' )


        ; e and b spec.
        pflux_setting = sc_info['pflux_setting']
        scale_info = pflux_setting['scale_info']
        e_var = edot0_fac_var
        e_vars = stplot_split(e_var)
        e_mor_vars = e_vars
        foreach var, e_vars, var_id do e_mor_vars[var_id] = stplot_mor_new(var, scale_info=scale_info)
        options, e_mor_vars, 'zrange', [1,1e4]
        options, e_mor_vars, 'ztitle', '(mV/m)!U2!N'
        options, e_mor_vars, 'constant', 1d/3*[1,2,3,4]

        b_var = b_fac_var
        b_vars = stplot_split(b_var)
        b_mor_vars = b_vars
        foreach var, b_vars, var_id do b_mor_vars[var_id] = stplot_mor_new(var, scale_info=scale_info)
        options, b_mor_vars, 'zrange', [1,1e4]
        options, b_mor_vars, 'ztitle', '(nT)!U2!N'
        
        vars = prefix+[['edot0','b1','u']+'_fac','pfdot0_fac_map',$
            'edot0_angle','b_tilt','b_gsm','b0_gsm','p_total','beta']
        tplot, vars, trange=time_range
        stop
        
        vars = [e_mor_vars,b_mor_vars, prefix+'pfdot0_fac_map']
        tplot, vars, trange=time_range
        stop


    ;---plot_vars options
        label_size = 0.8
    
        keo_var = 'thg_asf_keo'
        mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
        fmlt_var = prefix+'fmlt_'+internal_model+'_north'
        dmlt = 0.1
        mlat_range = [65,71]
    
        get_data, mlt_image_rect_var, times, data, limits=lim
        ntime = n_elements(times)
        model_index = where(models eq external_model)
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
            'zrange', [5e2,2e4], $
            'ztickv', [1e3,1e4], $
            'zticks', 1, $
            'ztickname', '10!U'+['3','4'], $
            'zminor', 9, $
            'zlog', 1, $
            'zticklen', zticklen, $
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
    
    
        density_var = prefix+'e_n'
        var = density_var
        options, var, 'yrange', [0.05,0.5]
        ;options, var, 'labels', 'e- Density!C  >200 eV'
        options, var, 'labels', ' '
        options, var, 'ytitle', '(cm!U-3!N)'
        options, var, 'ystyle', 9
    
        ; B tilt var.
        b_tilt_var = prefix+'b_tilt'
        add_setting, b_tilt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'arcsin(B!Dz!N/|B|)', $
            'unit', 'deg', $
            'yrange', [0,50], $
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
    
    
        ; Poynting flux.
        pf_var = prefix+'pfdot0_fac_map'
        var = pf_var
        options, var, 'yrange', [-250,100]
        options, var, 'ytickv', [-200,-100,0,100]
        options, var, 'yticks', 3
        options, var, 'yminor', 5
        options, var, 'constant', [0]
        
        ; Poynting flux spec.
        pf_spec_var = prefix+'pfdot0_fac_spec'
        copy_data, prefix+'pfdot0_fac_mor_spec_1', pf_spec_var
        var = pf_spec_var
        get_data, var, times, spec, ps
        fs = 1d3/ps
        spec *= 1e3
        store_data, var, times, spec, fs
        add_setting, var, smart=1, dictionary($
            'display_type', 'spec', $
            'short_name', 'S', $
            'unit', tex2str('mu')+'W/m!U2!N', $
            'zrange', [-1,1]*2, $
            'ztickv', [-2,0,2], $
            'zticks', 2, $
            'zminor', 2, $
            'zcharsize', label_size, $
            'color_table', 70, $
            'ytitle', '(mHz)', $
            'ylog', 1, $
            'yrange', [1,1e3], $
            'ytickv', [1e1,1e2], $
            'yticks', 1, $
            'yminor', 9, $
            'ytickname', ['10','10!U2'] )
    
    ;---Set plot_vars.
        sc_vars = [b_mor_vars[1],b_tilt_var,prefix+'e_en_spec',beta_var]
        asi_vars = 'thg_asf_keo'
        plot_vars = [asi_vars,sc_vars]
        nvar = n_elements(plot_vars)
        fig_labels = letters(nvar)+') '+['Aurora','B spec',tex2str('theta')+'!DB!N','e-','beta']
        ypans = [1.2,1,0.8,1,1]
    
       
    ;---Plot.
        fig_size = [6,5]
        sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
        margins = [12,4,10,1]
        poss = sgcalcpos(nvar, margins=margins, ypans=ypans, ypad=0.2)
        
        ; ticklen.
        uniform_ticklen = -ychsz*0.15*fig_size[0]
        for pid=0,nvar-1 do begin
            tpos = poss[*,pid]
            xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            options, plot_vars[pid], 'xticklen', xticklen
            options, plot_vars[pid], 'yticklen', yticklen
        endfor
        
        zticklen = uniform_ticklen/xchsz*1/fig_size[0]
        options, keo_var, 'zticklen', zticklen
        options, keo_var, 'zcharsize', label_size
        foreach var, prefix+['e_en_spec','i_en_spec'] do begin
            options, var, 'zticklen', zticklen
            options, var, 'zminor', 9
        endforeach
        
        options, pf_spec_var, 'zticklen', zticklen
        
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
        pid = where(plot_vars eq keo_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
    
            fmlat_var = prefix+'fmlat_'+internal_model+'_north'
            get_data, fmlat_var, times, fmlats
            model_index = (where(models eq external_model))[0]
            
            get_data, keo_var, limits=lim
            xrange = time_range
            yrange = lim.yrange
            plot, xrange, yrange, $
                xstyle=5, xrange=xrange, $
                ystyle=5, yrange=yrange, $
                position=tpos, nodata=1, noerase=1
            oplot, times, fmlats[*,model_index], color=sc_color, linestyle=2
            tx = time_range[0]
            ty = interpol(fmlats[*,model_index],times,tx)
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]+xchsz*0.5
            ty = tmp[1]+ychsz*0.3
            xyouts, tx,ty,normal=1, strupcase(sc_name+'-'+probe), charsize=label_size, alignment=0, color=sc_color
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
            
            filter = (sc_info['pflux_setting']).filter
            
    
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
        var = density_var
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
        
    endforeach
    
    
    return, plot_file    


end


print, fig_2017_0309_0700_ps_arc_v01(event_info=event_info)
end