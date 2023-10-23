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

        prefix = sc_info['prefix']
        probe = sc_info['probe']
        sc_name = sc_info['sc_name']
        sc_color = sc_info['sc_color']


        plot_file = join_path([plot_dir,$
            strlowcase('fig_'+event_info.id+'_ps_arc_'+sc_name+probe+'_'+version+'.pdf')])
        if keyword_set(test) then plot_file = 0


        ; Settings.
        fac_labels = ['||',tex2str('perp')+','+['west','out']]
        model_settings = sc_info['model_setting']
        external_model = 't89'
        internal_model = 'dipole'
        models = model_settings['models']
        
        
        
;        
;    ;---Test mor.
;        ; Calculate Edot0 angle.
;        b0_dsl_var = prefix+'b0_themis_dsl'
;        b0_dsl = get_var_data(b0_dsl_var, times=times)
;        edot0_angle_var = prefix+'edot0_angle'
;        edot0_angle = asin(b0_dsl[*,2]/snorm(b0_dsl))*constant('deg')
;        store_data, edot0_angle_var, times, edot0_angle
;        add_setting, edot0_angle_var, smart=1, dictionary($
;            'display_type', 'scalar', $
;            'short_name', 'edot0 angle', $
;            'unit', 'deg', $
;            'constant', [20,25] )
;        edot0_angle = get_var_data(prefix+'edot0_angle', times=times, in=time_range)
;        index = where(edot0_angle ge 25, count)
;        trs = times[time_to_range(index,time_step=1)]
;        durs = trs[*,1]-trs[*,0]
;        index = where(durs ge 180, ntr)
;        snapshot_trs = trs[index,*]
;        nsnapshot = n_elements(snapshot_trs[*,0])
;    
;        for pan_id=0,nsnapshot-1 do begin
;            the_tr = reform(snapshot_trs[pan_id,*])
;            time_step = 3d
;            time_step = 1d/8
;            common_times = make_bins(time_range, time_step)
;            b_var = themis_read_bfield(time_range, probe=probe, update=1, coord='gsm', id='fgl')
;            e_var = themis_read_efield(time_range, probe=probe, update=1, coord='gsm', id='survey', edot0_e56=1)
;            e_vars = stplot_split(prefix+'edot0_fac')
;            b_vars = stplot_split(prefix+'b1_fac')
;            e_var = e_vars[2]
;            b_var = b_vars[1]
;            
;            
;            pflux_setting = sc_info['pflux_setting']
;            scale_info = pflux_setting['scale_info']
;            
;            ; settings for wavelet transform.
;            s0 = time_step*4
;            dj = 1d/8
;            s1 = 1200d
;            j1 = floor(alog(s1/s0)/alog(2)/dj)
;            s1 = s0*2d^(dj*j1)
;            ns = j1+1
;            w0 = 6d
;            cdelta = 0.776d
;            psi0 = !dpi^(-0.25)
;        
;            foreach tvar, [e_var,b_var] do begin
;                interp_time, tvar, common_times
;                dat = get_var_data(tvar)
;                ;dat = snorm(dat)
;                index = where(finite(dat,nan=1), count)
;                if count ne 0 then dat[index] = 0
;                mor = wavelet(dat, time_step, pad=1, $
;                    s0=s0, dj=dj, j=j1, mother='Morlet', param=w0, period=ps, scale=ss)
;                psd = abs(mor)^2
;                idx = where(common_times ge the_tr[0] and common_times le the_tr[1], tnrec)
;                ;idx = where(uts ge time_range[0] and uts le time_range[1], tnrec)
;                psd = psd[idx,*]
;                gws = total(psd,1)/tnrec^2
;                ngws = (gws/ss)*(time_step*dj/cdelta)*tnrec
;                store_data, tvar+'_tmp', ps, [[gws],[ngws]]
;            endforeach
;        
;            get_data, e_var+'_tmp', ps, edat
;            get_data, b_var+'_tmp', ps, bdat
;            ebratio = sqrt(edat[*,0]/bdat[*,0])*1e3
;    
;    
;            avg_mass = 1
;            va0 = 22.0d     ; km/s, for B in nT, n in cc, m in atomic mass.
;            bmag = median(snorm(get_var_data(prefix+'b_gsm', in=the_tr)))
;            density = median(get_var_data(prefix+'e_n', in=the_tr))
;           ; density = 0.2
;            mhd_rho = density*avg_mass
;            va = va0*bmag/sqrt(mhd_rho)
;            
;            xxs = 1d3/ps
;            
;            f_sc = xxs*1e-3     ; in Hz.
;            f_g0 = 1.6e-19*1e-9/1.67e-27/2/!dpi   ; in Hz.
;            f_gi = f_g0*bmag/avg_mass
;            tavg = median(get_var_data(prefix+'p_tavg', in=the_tr))
;            vi = sqrt(tavg*1.6e-19/(avg_mass*1.67e-27))*1e-3    ; in km/s
;            vf = median(snorm(get_var_data(prefix+'p_vbulk_gsm', in=the_tr)))
;            ebr_theory = va*sqrt(1+(f_sc/f_gi*(vi/vf))^2)
;            stop
;            
;            plot, xxs, ebratio, ylog=1, xlog=1
;            oplot, xxs, ebr_theory, color=sgcolor('salmon'), linestyle=0
;
;        endfor
        

        
        
        ; Calculate B tilt.
        b_tilt_var = prefix+'b_tilt'
        b0_gsm_var = prefix+'b0_gsm'
        b0_gsm = get_var_data(b0_gsm_var, times=times)
        b_sm = cotran_pro(b0_gsm, times, 'gsm2sm')
        b_tilt = asin(b_sm[*,2]/snorm(b_sm))*constant('deg')
        store_data, b_tilt_var, times, b_tilt
        add_setting, b_tilt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'B tilt', $
            'unit', 'deg' )
        

        ; Calculate Edot0 angle.
        b0_dsl_var = prefix+'b0_themis_dsl'
        b0_dsl = get_var_data(b0_dsl_var, times=times)
        edot0_angle_var = prefix+'edot0_angle'
        edot0_angle = asin(b0_dsl[*,2]/snorm(b0_dsl))*constant('deg')
        store_data, edot0_angle_var, times, edot0_angle
        add_setting, edot0_angle_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'edot0 angle', $
            'unit', 'deg', $
            'constant', [20,25] )
        edot0_angle = get_var_data(prefix+'edot0_angle', times=times, in=time_range)
        index = where(edot0_angle ge 25, count)
        trs = times[time_to_range(index,time_step=1)]
        durs = trs[*,1]-trs[*,0]
        index = where(durs ge 180, ntr)
        snapshot_trs = trs[index,*]
        nsnapshot = n_elements(snapshot_trs[*,0])
            
        ; Velocity.
        u_gsm_var = prefix+'u_gsm'
        uts = make_bins(minmax(times),3)
        u_gsm = get_var_data(prefix+'p_vbulk_gsm', at=uts)
        store_data, u_gsm_var, uts, u_gsm
        add_setting, u_gsm_var, id='velocity', dictionary('coord', 'GSM')
        

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
        
        ; Bmag and B pressure.
        b_gsm_var = prefix+'b_gsm'
        b_gsm = get_var_data(b_gsm_var, at=times)
        bmag_var = prefix+'bmag'
        bmag = snorm(b_gsm)
        store_data, bmag_var, times, bmag
        add_setting, bmag_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', '|B|', $
            'unit', 'nT' )

        bp_var = prefix+'p_b'
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
            'yrange', [0.005,20], $
            'ytickv', [0.01,0.1,1,10], $
            'yticks', 3, $
            'ytickname', ['0.01','0.1','1','10'], $
            'yminor', 9, $
            'constant', [0.1,1], $
            'display_type', 'scalar', $
            'short_name', 'beta', $
            'ylog', 1, $
            'unit', '#' )

        ; alfven speed.
        va_var = prefix+'va'
        va0 = 1e-9/sqrt(1e6*!dpi*4e-7*1.67e-27)*1e-3    ; km/s, B in nT, n in cc, m in m_p.
        num_dens = get_var_data(prefix+'e_n', at=times)
        avg_mass = 1d
        va = va0*bmag/sqrt(num_dens*avg_mass)
        store_data, va_var, times, va
        add_setting, va_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'v!DA!N', $
            'unit', 'km/s' )


        ; pflux, e and b spec.
        pflux_setting = sc_info['pflux_setting']
        scale_info = pflux_setting['scale_info']

        b_fac_var = prefix+'b1_fac'
        b_vars = stplot_split(b_fac_var)
        the_b_var = b_vars[1]
        b_mor_var = stplot_mor_new(the_b_var, scale_info=scale_info)
        var = b_mor_var
        options, var, 'ztitle', '(nT)!U2!N'
        options, var, 'zrange', [1,1e4]
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        options, var, 'yrange', lim.yrange*1e3
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', 10^[1d,2,3]

        e_fac_var = prefix+'edot0_fac'
        e_vars = stplot_split(e_fac_var)
        the_e_var = e_vars[2]
        e_mor_var = stplot_mor_new(the_e_var, scale_info=scale_info)
        var = e_mor_var
        options, var, 'ztitle', '(mV/m)!U2!N'
        options, var, 'zrange', [1,1e4]
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        options, var, 'yrange', lim.yrange*1e3
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', 10^[1d,2,3]

        pf_fac_var = prefix+'pfdot0_fac'
        pf_vars = stplot_split(pf_fac_var)
        the_pf_var = pf_vars[0]
        pf_mor_var = stplot_mor_new(the_pf_var, scale_info=scale_info)
        options, var, 'ztitle', '('+tex2str('mu')+'W/m!U2!N)!U2!N'
        options, var, 'zrange', [1,1e4]
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        options, var, 'yrange', lim.yrange*1e3
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', 10^[1d,2,3]



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
    
    
            ; Add FMLat.
            fmlat_var = prefix+'fmlat_'+internal_model+'_north'
            the_fmlat_var = prefix+'fmlat_north_smooth'
            dfmlat = [0d,1]
            model_index = (where(models eq external_model))[0]
            get_data, fmlat_var, uts, fmlats
            fmlats = fmlats[*,model_index]
            store_data, the_fmlat_var, uts, fmlats
;            window = 600.
;            width = window/sdatarate(uts)
;            fmlats = smooth(fmlats, width, edge_mirror=1, nan=1)
                        
;            get_data, keo_var, times, keo, mlats, limits=lim
;            fmlats = get_var_data(the_fmlat_var, at=times)
;            nmlat = n_elements(mlats)
;            nmlat1 = nmlat*10
;            mlats = rebin(mlats,nmlat1)
;            keo = rebin(keo,ntime,nmlat1)
;            
;            transparency = 80
;            c_trans = transparency*0.01
;            bg_color = sgcolor('white')
;            ct = lim.color_table
;            foreach time, times, time_id do begin
;                index = where_pro(mlats, '[]', fmlats[time_id]+dfmlat, count=count)
;                if count eq 0 then continue
;                keo[time_id,index] *= c_trans
;            endforeach
;            store_data, keo_var, times, keo, mlats



;        edot0_fac_var = prefix+'edot0_fac'
;        vars = stplot_split(edot0_fac_var)
;        e_var = edot0_fac_var+'_comp3'
;        var = e_var
;        add_setting, var, smart=1, dictionary($
;            'display_type', 'scalar', $
;            'short_name', 'E', $
;            'unit', 'mV/m' )
;        options, var, 'yrange', [-1,1]*120
;        options, var, 'ytickv', [-1,1]*100
;        options, var, 'yticks', 2
;        options, var, 'yminor', 5
;        options, var, 'labels', 'E!D'+tex2str('perp')+',out'
    
    
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
            'yrange', [5,55], $
            'ytickv', [10,30,50], $
            'yticks', 2, $
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
        
        
        ; Morlet spec.
        var = b_mor_var
        ct_mor = 64

        ztitle = 'dB!D'+fac_labels[1]+'!N (nT!U2!N)'
        zrange = [1,1e4]
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,1,inner=1)
        ztickv = 10.^log_ztickv
        zticks = n_elements(ztickv)-1
        ztickn = '10!U'+string(log_ztickv,format='(I0)')
        ztickn[0:*:2] = ' '
        yrange = get_setting(var, 'yrange')
        yrange = [1,3e3]
        log_yrange = alog10(yrange)
        log_ytickv = make_bins(log_yrange,1,inner=1)
        ytickv = 10.^log_ytickv
        yticks = n_elements(ytickv)-1
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        index = where(ytickn eq '10!U1', count)
        if count ne 0 then ytickn[index] = '10'
        index = where(ytickn eq '10!U0', count)
        if count ne 0 then ytickn[index] = '1'
        ;ytickn[1:*:2] = ' '
        ;ytickn[0] = ' '

        options, var, 'ztitle', ztitle
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'ztickname', ztickn
        options, var, 'minor', 9
        options, var, 'zcharsize', label_size
        options, var, 'color_table', ct_mor
    
        options, var, 'yrange', yrange
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'ytickname', ytickn
        options, var, 'yminor', 9
    
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
        sc_vars = [b_mor_var,b_tilt_var,beta_var]
        asi_vars = 'thg_asf_keo'
        plot_vars = [asi_vars,sc_vars]
        nvar = n_elements(plot_vars)
        fig_labels = letters(nvar)+') '+['Aurora','Wave','B tilt','Plasma '+tex2str('beta')]
        ypans = [1,1,0.7,0.7]
    
       
    ;---Plot.
        fig_size = [6,5.5]
        sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
        uniform_ticklen = -ychsz*0.2*fig_size[0]
        margins = [12,5,10,1]
        all_poss = sgcalcpos(nvar+1, margins=margins, ypans=[ypans,0.8], ypad=[fltarr(nvar-1)+0.4,4])
        poss = all_poss[*,0:nvar-1]
        
        tpos = poss[*,-1]
        plot, time_range, [0,1], $
            xstyle=5, ystyle=5, nodata=1, noerase=1, position=tpos
        thick = (keyword_set(test))? 1: 4
        for ii=0,ntr-1 do begin
            txs = reform(snapshot_trs[ii,*])
            tys = fltarr(2)
            for jj=0,1 do begin
                tmp = convert_coord(txs[jj],0, data=1, to_normal=1)
                txs[jj] = tmp[0]
                tys[jj] = tmp[1]
            endfor
            plots, txs,tys-ychsz*0.2,normal=1, color=sgcolor('red'), thick=thick
            tx = mean(txs)
            ty = tys[0]-ychsz*1.1
            msg = 'T'+string(ii+1,format='(I0)')
            xyouts, tx,ty, msg, normal=1, alignment=0.5
        endfor
        label_yshift = -ychsz*0.5

        
        low_poss = sgcalcpos(1,nsnapshot*2, position=all_poss[*,nvar], xpad=[1,3,1])
        f_window = minmax([1d3/sc_info['b0_window'],1d3/3])
        for pan_id=0,nsnapshot-1 do begin
            tpos = low_poss[*,pan_id*2]
            xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            msg = (letters(nvar+pan_id))[-1]+'-'+string(1,format='(I0)')+')'
            xyouts, tx,ty,normal=1, msg
            
            xtitle = 'E!D'+fac_labels[2]+'!N/B!D'+fac_labels[1]+'!N!C(km/s)'
            ytitle = 'Freq!C(mHz)'
            the_tr = reform(snapshot_trs[pan_id,*])
            ;the_tr = mean(the_tr)+[-1,1]*60
            

            pow = get_var_data(e_mor_var, fs, times=times, in=the_tr, limits=lim)
            ntime = n_elements(times)
            e_gws = total(pow,1)/ntime
            pow = get_var_data(b_mor_var, fs, times=times, in=the_tr, limits=lim)
            ntime = n_elements(times)
            b_gws = total(pow,1)/ntime
            ebratio = sqrt(e_gws/b_gws)*1e3

            ; Calc V_alfven.
            avg_mass = 1
            va0 = 22.0d     ; km/s, for B in nT, n in cc, m in atomic mass.
            bmag = median(snorm(get_var_data(prefix+'b_gsm', in=the_tr)))
            density = median(get_var_data(prefix+'e_n', in=the_tr))
            mhd_rho = density*avg_mass
            va = va0*bmag/sqrt(mhd_rho)

            yys = fs
            xxs = ebratio
            xrange = [200,2e5]
            yrange = get_setting(b_mor_var, 'yrange')
            ytickv = get_setting(b_mor_var, 'ytickv')
            yticks = get_setting(b_mor_var, 'yticks')
            ytickn = get_setting(b_mor_var, 'ytickname')
            yminor = get_setting(b_mor_var, 'yminor')
            ytickformat = ''
            
            if pan_id ne 0 then begin
                ytitle = ' '
                ytickformat = '(A1)'
            endif
            
            plot, xxs, yys, $
                xstyle=1, xlog=1, xtitle=xtitle, xrange=xrange, $
                ystyle=1, ylog=1, ytitle=ytitle, yrange=yrange, $
                ytickv=ytickv, yticks=yticks, ytickname=ytickn, yminor=yminor, ytickformat=ytickformat, $
                position=tpos, noerase=1, nodata=1, $
                xticklen=xticklen, yticklen=yticklen
            index = where_pro(fs, '[', f_window[0])
            oplot, xxs[index], yys[index]
            ;oplot, xxs[0:index[0]], yys[0:index[0]], color=sgcolor('silver')
            oplot, [0,0]+va, yrange, linestyle=1, color=sgcolor('salmon')
            tmp = convert_coord(va, yrange[0], data=1, to_normal=1)
            tx = tmp[0]+xchsz*0.5
            ty = tmp[1]+ychsz*0.5
            msg = 'v!DA!N'
            xyouts, tx,ty, msg, normal=1, color=sgcolor('salmon')
            
            
            f_sc = fs*1e-3     ; in Hz.
            f_g0 = 1.6e-19*1e-9/1.67e-27/2/!dpi   ; in Hz.
            f_gi = f_g0*bmag/avg_mass
            tavg = median(get_var_data(prefix+'p_tavg', in=the_tr))
            vi = sqrt(tavg*1.6e-19/(avg_mass*1.67e-27))*1e-3    ; in km/s
            vf = median(snorm(get_var_data(prefix+'u_gsm', in=the_tr)))
            ebr_theory = va*sqrt(1+(f_sc/f_gi*(vi/vf))^2)
            oplot, ebr_theory, fs, color=sgcolor('salmon'), linestyle=0
            
            tx = (tpos[0]+low_poss[2,pan_id*2+1])*0.5
            ty = tpos[3]+ychsz*0.4
            msg = 'T'+string(pan_id+1,format='(I0)')+': '+strjoin(time_string(the_tr,tformat='hh:mm'),'-')+' UT'
            xyouts, tx,ty,normal=1, msg, charsize=label_size, alignment=0.5
            
            
            ; pflux power.
            tpos = low_poss[*,pan_id*2+1]
            xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

            tx = tpos[2]-xchsz*0.5
            ty = tpos[3]-ychsz*1
            msg = (letters(nvar+pan_id))[-1]+'-'+string(2,format='(I0)')+')'
            xyouts, tx,ty,normal=1, msg, alignment=1
            
            ytickformat = '(A1)'
            ytitle = ' '
            xtitle = 'Power S!D'+fac_labels[0]+'!C(mW/m!U2!N)'
            
            pow = get_var_data(pf_mor_var, fs, times=times, in=the_tr, limits=lim)
            ntime = n_elements(times)
            pf_gws = total(pow,1)/ntime
            pf_mor_info_var = pf_mor_var+'_fft_info'
            get_data, pf_mor_info_var, tmp, info
            psd = 2*info.c_tau*info.dt/info.cdelta*pf_gws*1e3
            fs = info.fs*1e3*0.5

            index = where_pro(fs, '[', f_window[0])
            xxs = psd[index]
            yys = fs[index]
            xrange = minmax(xxs)
            log_xrange = minmax(make_bins(alog10(xrange),1))
            log_xrange = [-6d,2]
            xrange = 10^log_xrange
            log_xtickv = make_bins(log_xrange,1)
            xtickv = 10^log_xtickv
            xticks = n_elements(xtickv)-1
            xminor = 9
            xtickn = '10!U'+string(log_xtickv,format='(I0)')
            xtickn[0:*:2] = ' '
            xtickformat = ''
            
            plot, xrange, yrange, $
                xstyle=1, xlog=1, xtitle=xtitle, xrange=xrange, $
                ystyle=1, ylog=1, ytitle=ytitle, yrange=yrange, $
                xtickv=xtickv, xticks=xticks, xtickname=xtickn, xminor=xminor, xtickformat=xtickformat, $
                ytickv=ytickv, yticks=yticks, ytickname=ytickn, yminor=yminor, ytickformat=ytickformat, $
                position=tpos, noerase=1, nodata=1, $
                xticklen=xticklen, yticklen=yticklen
            
            oplot, xxs, yys
        endfor

        
        
        
        ; ticklen.
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
        foreach var, prefix+['e_en_spec','p_en_spec','b1_fac_comp2_mor'] do begin
            options, var, 'zticklen', zticklen
            options, var, 'zminor', 9
        endforeach
        
        options, pf_spec_var, 'zticklen', zticklen
        
        tplot, plot_vars, position=poss, trange=time_range, noerase=1
        for pid=0,nvar-1 do begin
            tpos = poss[*,pid]
            tx = tpos[0]-xchsz*10
            ty = tpos[3]+label_yshift
            xyouts, tx,ty, normal=1, fig_labels[pid]
        endfor
        ;timebar, , color=sgcolor('red'), linestyle=1


        ; Add FMLat.
        pid = where(plot_vars eq keo_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
    
            fmlats = get_var_data(the_fmlat_var, times=times)
            get_data, keo_var, limits=lim
            xrange = time_range
            yrange = lim.yrange
            plot, xrange, yrange, $
                xstyle=5, xrange=xrange, $
                ystyle=5, yrange=yrange, $
                position=tpos, nodata=1, noerase=1
;            alpha = 0.5
;            white = sgcolor('white')
;            the_color = alpha*sc_color+(1-alpha)*white
            the_color = sc_color
            oplot, times, fmlats+dfmlat[0], color=the_color, linestyle=2
            oplot, times, fmlats+dfmlat[1], color=the_color, linestyle=2
            tx = time_range[0]
            ty = interpol(fmlats,times,tx)+mean(dfmlat)
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]+xchsz*0.5
            ty = tmp[1]-ychsz*0.3
            xyouts, tx,ty,normal=1, strupcase(sc_name+'-'+probe), charsize=label_size, alignment=0, color=sc_color
        endif

        
    
        if keyword_set(test) then stop
        sgclose
        
    endforeach
    
    
    return, plot_file    


end


print, fig_2017_0309_0700_ps_arc_v01(event_info=event_info)
end