;+
; Part of fig_2015_0416_0800_rbsp_v01
;;+
; Plot ewo, pflux, dp.
;-

function to_power_scale, data, power
    if n_elements(power) eq 0 then power = 0.5
    return, abs(data)^power*sign(data)
end
function from_power_scale, data, power
    if n_elements(power) eq 0 then power = 0.5
    return, abs(data)^(1d/power)*sign(data)
end

function fig_2015_0416_0800_extended_data_pflux_v01_rbsp_panel, rbsp_poss, event_info=event_info

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)


    time_range = time_double(['2015-04-16/08:00','2015-04-16/08:25'])
    ;time_range = event_info['time_range']

    rbsp_info = event_info.rbsp.rbspa
    model_setting = rbsp_info['model_setting']
    all_models = model_setting['models']
    internal_model = event_info.internal_model
    external_model = event_info.external_model
    model_index = where(all_models eq external_model)

    dp_times = event_info['dp_times']
    dp_mlts = fltarr(2)
    foreach key, sort_uniq(event_info.rbsp.keys()), probe_id do begin
        info = (event_info.rbsp)[key]
        prefix = info['prefix']
        probe = info['probe']
        sc_color = info['sc_color']
        fmlt_var = prefix+'fmlt_'+internal_model+'_north'
        tx = dp_times[probe_id]
        dp_mlts[probe_id] = (get_var_data(fmlt_var, at=tx))[model_index]
    endforeach
    dp_omega = total(dp_mlts*[-1,1])/total(dp_times*[-1,1])*60*15
    event_info['dp_omega'] = dp_omega

    power = 1d/4
    label_size = 0.8
    tmp = smkarthm(0,2*!dpi,20,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1

;---B spec.
    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    foreach info, event_info.rbsp do begin
        prefix = info['prefix']
        probe = info['probe']

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
        
        time_step = info['field_time_step']
        b0_window = info['b0_window']
        yrange = minmax(1d3/[b0_window*0.5,time_step*4])
        log_ytickv = [1,2,3]
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
        options, var, 'zcharsize', label_size
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        options, var, 'yrange', yrange
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'ytickname', ytickn
        options, var, 'yminor', yminor
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', ytickv
        options, var, 'zcharsize', label_size
        
        
        ; Pflux.
        pf_spinfit_var = prefix+'pfdot0_fac_spinfit_phasef_map'
        var = pf_spinfit_var
        options, var, 'yrange', [-60,30]
        options, var, 'ytickv', [-60,-20,20]
        options, var, 'yticks', 2
        options, var, 'yminor', 2
        options, var, 'constant', [0]
        pf_survey_var = prefix+'pfdot0_fac_survey_map'
        if probe eq 'a' then pf_survey_var = prefix+'pfdot0_fac_v24_map'
        
        var = prefix+'pfdot0_fac_spinfit_map'
        copy_data, pf_spinfit_var, var
        var = prefix+'pfdot0_fac_survey_map'
        copy_data, pf_survey_var, var
        options, var, 'yrange', [-360,-60]
        options, var, 'ytickv', [-360,-180]
        options, var, 'yticks', 1
        options, var, 'yminor', 3
        
        vars = prefix+'pfdot0_fac_'+['spinfit','survey']+'_map'
        options, vars, 'ytitle', ' '
        options, vars, 'labels', [' ',' ',' ']
        
        var = prefix+'pfdot0_fac_survey_map'
        pfdot0 = get_var_data(var, times=times, limits=lim)
        pfpara = pfdot0[*,0]
        if probe eq 'a' then pfpara[*] = !values.f_nan
        var = prefix+'pfdot0_fac_survey_map_para'
        store_data, var, times, pfpara, limits=lim
        options, var, 'labels', ' '
        options, var, 'colors', sgcolor('light_salmon')
    endforeach


;---Plot.
    rbsp_letters = ['a','b','c','d']
    foreach key, sort_uniq(event_info.rbsp.keys()), rbsp_id do begin
        sc_info = (event_info.rbsp)[key]
        rbsp_letter = rbsp_letters[rbsp_id*2]
        rbsp_letter2 = rbsp_letters[rbsp_id*2+1]
        prefix = sc_info['prefix']
        probe = sc_info['probe']
        sc_color = sc_info['sc_color']
        
        plot_info = orderedhash()
        plot_info[prefix+'b1_fac_comp2_mor'] = dictionary($
            'fig_label', 'Wave', $
            'ypan', 1.0 )
        plot_info[prefix+'pfdot0_fac_spinfit_map'] = dictionary($
            'fig_label', 'S FAC', $
            'ypan', 0.6, $
            'ypad', 0. )
        plot_info[prefix+'pfdot0_fac_survey_map_para'] = dictionary($
            'fig_label', ' ', $
            'ypan', 0.4 )
        
        my_pos = rbsp_poss[*,rbsp_id]
        poss = my_pos

        fig_size = double([!d.x_size,!d.y_size])
        tmp = get_charsize()
        xchsz = tmp[0]
        ychsz = tmp[1]
        uniform_ticklen = -ychsz*0.15*fig_size[1]
    
        
        tpos = poss[*,-1]
        msg = strupcase('RBSP-'+probe)
        tx = total(tpos[[0,2]])*0.5
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,normal=1, msg, color=sc_color, alignment=0.5
        

        ; Add pflux spec and ebratio.
        low_pos = poss[*,-1]
        ;low_pos[[0,2]] += [0,2]*xchsz
        low_poss = sgcalcpos(1,2, position=low_pos, xpad=2)

        ; ebratio.
        tpos = low_poss[*,0]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

        e_fac_var = prefix+'edot0_fac_spinfit_phasef'
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

        pow = get_var_data(e_mor_var, fs, times=times, in=dp_tr, limits=lim)
        ntime = n_elements(times)
        e_gws = total(pow,1)/ntime
        pow = get_var_data(b_mor_var, fs, times=times, in=dp_tr, limits=lim)
        ntime = n_elements(times)
        b_gws = total(pow,1)/ntime
        ebratio = sqrt(e_gws/b_gws)*1e3

        xtitle = 'E!D'+fac_labels[2]+'!N/B!D'+fac_labels[1]+'!N (km/s)'
        ytitle = 'Freq!C(mHz)'

        ; Calc V_alfven.
        np = get_var_data(prefix+'p_n', in=dp_tr)
        no = get_var_data(prefix+'o_n', in=dp_tr)
        o_ratio = mean(no/(no+np),nan=1)
        p_ratio = 1-o_ratio
        p_mass = 1d
        o_mass = 16d
        avg_mass = p_ratio*p_mass+o_ratio*o_mass
        print, avg_mass
        
        va0 = 22.0d     ; km/s, for B in nT, n in cc, m in atomic mass.
        bmag = median(snorm(get_var_data(prefix+'b_gsm', in=dp_tr)))
        ;var = rbsp_read_density(time_range, probe=probe, id='emfisis')
        tmp = get_var_data(prefix+'density_combo', in=dp_tr, limits=lim)
        tmp = tmp[*,where(lim.labels eq 'N EMFISIS', count)]
        if count eq 0 then message, 'Inconsistency ...'
        density = median(tmp)
        mhd_rho = density*avg_mass
        va = va0*bmag/sqrt(mhd_rho)

        yys = fs
        xxs = ebratio
        xrange = [20,2e5]
        yrange = get_setting(b_mor_var, 'yrange')
        ytickv = get_setting(b_mor_var, 'ytickv')
        yticks = get_setting(b_mor_var, 'yticks')
        ytickn = get_setting(b_mor_var, 'ytickname')
        yminor = get_setting(b_mor_var, 'yminor')
        ytickformat = ''
        
        
        plot, xxs, yys, $
            xstyle=1, xlog=1, xtitle=xtitle, xrange=xrange, $
            ystyle=1, ylog=1, ytitle=ytitle, yrange=yrange, $
            ytickv=ytickv, yticks=yticks, ytickname=ytickn, yminor=yminor, ytickformat=ytickformat, $
            position=tpos, noerase=1, nodata=1, $
            xticklen=xticklen, yticklen=yticklen
        f_spin = 1d3/11
        index = where_pro(fs, ']', f_spin)
        oplot, xxs[index], yys[index]
        oplot, xrange, f_spin+[0,0], linestyle=1

        
        if probe eq 'b' then begin
            e_fac_var = prefix+'edot0_fac_survey'
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

            pow = get_var_data(e_mor_var, fs, times=times, in=dp_tr, limits=lim)
            ntime = n_elements(times)
            e_gws = total(pow,1)/ntime
            pow = get_var_data(b_mor_var, fs, times=times, in=dp_tr, limits=lim)
            ntime = n_elements(times)
            b_gws = total(pow,1)/ntime
            ebratio = sqrt(e_gws/b_gws)*1e3
            
            yys = fs
            xxs = ebratio
            oplot, xxs, yys, color=sgcolor('light_salmon')

            f_sc = fs*1e-3     ; in Hz.
            f_g0 = 1.6e-19*1e-9/1.67e-27/2/!dpi   ; in Hz.
            bmag = median(snorm(get_var_data(prefix+'b_gsm', in=dp_tr)))
            f_gi = f_g0*bmag/avg_mass
            tavg = median(get_var_data(prefix+'p_t', in=dp_tr))
            vi = sqrt(tavg*1.6e-19/(avg_mass*1.67e-27))*1e-3    ; in km/s
            ;vf = median(snorm(get_var_data(prefix+'u_gsm', in=dp_tr)))
            ravg = median(snorm(get_var_data(prefix+'r_gsm', in=dp_tr)))
            dp_omega = event_info['dp_omega']
            vf = dp_omega*constant('rad')/60*ravg*constant('re')
            ebr_theory = va*sqrt(1+(f_sc/f_gi*(vi/vf))^2)
            oplot, ebr_theory, fs, color=sgcolor('light_salmon'), linestyle=0
        endif

        ;oplot, xxs[0:index[0]], yys[0:index[0]], color=sgcolor('silver')
        oplot, [0,0]+va, yrange, linestyle=1, color=sgcolor('red')
        tmp = convert_coord(va, yrange[1], data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = 'v!DA!N'
        xyouts, tx,ty, msg, normal=1, color=sgcolor('red')
        tx = tpos[0]-xchsz*7
        label_yshift = -0.6*ychsz
        ty = tpos[3]+label_yshift
        msg = rbsp_letter+') E/B'
        xyouts, tx,ty, msg, normal=1
        

        ; pflux power.
        tpos = low_poss[*,1]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        ytitle = ' '
        xtitle = 'Power S!D'+fac_labels[0]+'!N (mW/m!U2!N)!U2'
        ytickformat = '(A1)'
        
        pf_var1 = prefix+'pfdot0_fac_spinfit_phasef'
        pf_vars = stplot_split(pf_var1)
        the_pf_var = pf_vars[0]
        pf_mor_var = stplot_mor_new(the_pf_var, scale_info=scale_info)
        options, var, 'ztitle', '('+tex2str('mu')+'W/m!U2!N)!U2!N'
        options, var, 'zrange', [1,1e4]
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        ;options, var, 'yrange', lim.yrange*1e3
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', 10^[1d,2,3]
        
        pow = get_var_data(pf_mor_var, fs, times=times, in=dp_tr, limits=lim)
        ntime = n_elements(times)
        pf_gws = total(pow,1)/ntime
        info = lim.cwt_info
        psd = 2*info.c_tau*info.dt/info.cdelta*pf_gws*1e3
        fs = info.fs*1e3*0.5

        xxs = psd
        yys = fs
        xrange = minmax(xxs)
        log_xrange = minmax(make_bins(alog10(xrange),1))
        log_xrange = [-3d,4]
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
        
        f_spin = 1d3/11
        index = where_pro(fs, ']', f_spin)
        oplot, xxs[index], yys[index]
        oplot, xrange, f_spin+[0,0], linestyle=1
        tmp = convert_coord(xrange[1],f_spin, data=1, to_normal=1)
        msg = 'f!Dspin!N'
        tx = tmp[0]+xchsz*1
        ty = tmp[1]-ychsz*0.3
        if probe eq 'b' then xyouts, tx,ty, msg, normal=1, charsize=label_size, alignment=0
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = rbsp_letter2+') Power'
        xyouts, tx,ty, msg, normal=1, alignment=1
        
        
        if probe eq 'b' then begin
            pf_var2 = prefix+'pfdot0_fac_survey'
            pf_vars = stplot_split(pf_var2)
            the_pf_var = pf_vars[0]
            pf_mor_var = stplot_mor_new(the_pf_var, scale_info=scale_info)
            options, var, 'ztitle', '('+tex2str('mu')+'W/m!U2!N)!U2!N'
            options, var, 'zrange', [1,1e4]
            get_data, var, times, mor, fs, limits=lim
            store_data, var, times, mor, fs*1e3
            ;options, var, 'yrange', lim.yrange*1e3
            options, var, 'ytitle', 'Freq!C(mHz)'
            options, var, 'constant', 10^[1d,2,3]

            pow = get_var_data(pf_mor_var, fs, times=times, in=dp_tr, limits=lim)
            ntime = n_elements(times)
            pf_gws = total(pow,1)/ntime
            info = lim.cwt_info
            psd2 = 2*info.c_tau*info.dt/info.cdelta*pf_gws*1e3
            fs2 = info.fs*1e3*0.5

            xxs = psd2
            yys = fs2
            
            oplot, xxs, yys, color=sgcolor('light_salmon')
        endif

    endforeach

end






function fig_2015_0416_0800_extended_data_pflux_v01, event_info=event_info, test=test

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
    event_info['dp_times'] = time_double('2015-04-16/'+['08:07:42','08:10:00'])

;---Plot file.
    test_plot_panels = 0
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,$
        strlowcase('fig_'+event_info.id+'_extended_data_pflux_'+version+'.pdf')])
    if keyword_set(test) then plot_file = 0


;---Figure out panel size.
    tmp = get_abs_chsz()
    abs_xchsz = tmp[0]
    abs_ychsz = tmp[1]
    uniform_ticklen = -abs_ychsz*0.15

    fig_size = [6,1.8]
    sgopen, plot_file, size=fig_size

    bottom_margins = [9,4,3,2]
    bottom_poss = sgcalcpos(1,2, xpad=10, margins=bottom_margins)
    left_pos = bottom_poss[*,0]
    right_pos = bottom_poss[*,1]
    rbsp_poss = [[left_pos],[right_pos]]
    

;--Panels.

    tpos = fig_2015_0416_0800_extended_data_pflux_v01_rbsp_panel(rbsp_poss, event_info=event_info)

    if keyword_set(test) then stop
    sgclose
    
    return, plot_file

end

test = 1
print, fig_2015_0416_0800_extended_data_pflux_v01(event_info=event_info, test=test)
end