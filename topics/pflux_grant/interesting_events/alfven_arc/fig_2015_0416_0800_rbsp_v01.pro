;+
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

function fig_2015_0416_0800_rbsp_v01_rbsp_panel, rbsp_poss, event_info=event_info

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)


    time_range = time_double(['2015-04-16/08:00','2015-04-16/08:25'])
    ;time_range = event_info['time_range']
    dp_times = event_info['dp_times']
    
    power = 1d/4
    label_size = 0.8
    tmp = smkarthm(0,2*!dpi,20,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1


;---keogram.
    dmlt = 0.1
    mlat_range = [60,68]
    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    keo_zlog = 1
    if keo_zlog eq 0 then begin
        zrange = [-1,1]*1e3
        zrange = to_power_scale([-1,1]*1e4, power)
        ct = 70
        ztickv = [-1,0,1]*abs(max(zrange))
        zticks = n_elements(ztickv)-1
        zminor = 10
        ztickn = string(ztickv,format='(I0)')
    endif else begin
        zrange = [1e1,1e4]
        ct = 49
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,1,inner=1)
        ztickv = 10^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        ztickn = '10!U'+string(log_ztickv,format='(I0)')
        index = where(ztickn eq '10!U1', count)
        if count ne 0 then ztickn[index] = '10'
        index = where(ztickn eq '10!U0', count)
        if count ne 0 then ztickn[index] = '1'
    endelse

    foreach info, event_info.rbsp do begin
        prefix = info['prefix']
        keo_var = prefix+'keo'

        model_setting = info['model_setting']
        all_models = model_setting['models']
        internal_model = event_info.internal_model
        external_model = event_info.external_model
        model_index = where(all_models eq external_model)
        
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
        if keo_zlog eq 1 then begin
            hehe = hehe>0.01
        endif else begin
            ;hehe = to_power_scale(hehe, power)
            hehe = hehe
        endelse
        store_data, keo_var, times, hehe, mlat_bins
        add_setting, keo_var, smart=1, dictionary($
            'display_type', 'spec', $
            'short_name', 'ASI Count', $
            'unit', '#', $
            'ytitle', 'MLat!C(deg)', $
            'yrange', yrange, $
            'ytickv', ytickv, $
            'yticks', yticks, $
            'yminor', yminor, $
            'zlog', keo_zlog, $
            'zrange', zrange, $
            'ztickv', ztickv, $
            'zticks', zticks, $
            'zminor', zminor, $
            'ztickname', ztickn, $
            'color_table', ct )
    endforeach

;---B spec.
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
    rbsp_letters = ['c','d']
    foreach key, sort_uniq(event_info.rbsp.keys()), rbsp_id do begin
        sc_info = (event_info.rbsp)[key]
        rbsp_letter = rbsp_letters[rbsp_id]
        prefix = sc_info['prefix']
        probe = sc_info['probe']
        sc_color = sc_info['sc_color']
        
        plot_info = orderedhash()
;        plot_info[prefix+'keo'] = dictionary($
;            'fig_label', 'Keo', $
;            'ypan', 1.0 )
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
        
        plot_vars = (plot_info.keys()).toarray()
        nvar = n_elements(plot_vars)
        fig_labels = strarr(nvar)
        ypans = fltarr(nvar)
        ypads = fltarr(nvar)
        foreach key, plot_info.keys(), var_id do begin
            info = plot_info[key]
            fig_labels[var_id] = info['fig_label']
            if probe eq 'b' then fig_labels[var_id] = ' '
            if ~info.haskey('ypan') then info['ypan'] = 1.
            if ~info.haskey('ypad') then info['ypad'] = 0.4
            ypans[var_id] = info['ypan']
            ypads[var_id] = info['ypad']
        endforeach
        ypads = ypads[0:nvar-2]
        fig_labels = rbsp_letter+'-'+string(findgen(nvar)+1,format='(I0)')+') '+fig_labels
        fig_labels[-1] = ' '

        my_pos = rbsp_poss[*,rbsp_id]
        poss = sgcalcpos(nvar+1, margins=margins, ypans=[ypans,1], ypad=[ypads,3], position=my_pos)

        fig_size = double([!d.x_size,!d.y_size])
        tmp = get_charsize()
        xchsz = tmp[0]
        ychsz = tmp[1]
        uniform_ticklen = -ychsz*0.15*fig_size[1]
    
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
        
        ; Add labels for pflux.
        pf_var1 = prefix+'pfdot0_fac_spinfit_map'
        pf_var2 = prefix+'pfdot0_fac_survey_map'
        pid = where(plot_vars eq pf_var1, count)
        if count ne 0 then begin
            pos1 = poss[*,pid]
            pos2 = poss[*,pid+1]
            tpos = pos1
            tpos[1] = pos2[1]


            labels = 'S!D'+fac_labels+'!N'
            colors = constant('rgb')
            nlabel = n_elements(labels)
            dy = (total(tpos[[1,3]]*[-1,1])-ychsz*1)/(nlabel)
            label_poss = tpos[3]-findgen(nlabel+1)*dy
            for ii=0,nlabel-1 do begin
                tx = tpos[2]+xchsz*0.5
                ty = label_poss[ii]-ychsz*0.8
                msg = labels[ii]
                xyouts, tx,ty,msg, normal=1, color=colors[ii]
            endfor

            if probe eq 'b' then begin
                get_data, 'rbspb_pfdot0_fac_survey_map', times, pf_fac
                pf_para = pf_fac[*,0]
                the_color = colors[0]
                the_color = sgcolor('light_salmon')
                the_linestyle = 0

                get_data, pf_var1, limits=lim
                xrange = time_range
                yrange = lim.yrange
                plot, xrange, yrange, $
                    xstyle=5, xrange=xrange, $
                    ystyle=5, yrange=yrange, $
                    position=pos1, nodata=1, noerase=1
                oplot, times, pf_para, color=the_color, linestyle=the_linestyle

                get_data, pf_var2, limits=lim
                xrange = time_range
                yrange = lim.yrange
                plot, xrange, yrange, $
                    xstyle=5, xrange=xrange, $
                    ystyle=5, yrange=yrange, $
                    position=pos2, nodata=1, noerase=1
                oplot, times, pf_para, color=the_color, linestyle=the_linestyle
                msg = 'S!D||!N hires'
                tx = pos2[2]+xchsz*0.5
                ty = label_poss[-1]-ychsz*0.8
                xyouts, tx,ty,msg, normal=1, color=the_color
            endif
        endif

        tplot_options, 'version', 2
        tplot, plot_vars, position=poss, trange=time_range, noerase=1, single_line_uttick=1
        dp_tr = dp_times[rbsp_id]+[0,1.5]*60
        timebar, dp_tr, linestyle=1, color=sc_color
        label_yshift = -ychsz*0.7
        for pid=0,nvar-1 do begin
            tpos = poss[*,pid]
            tx = tpos[0]-xchsz*10
            if probe eq 'b' then tx = tpos[0]-xchsz*6
            ty = tpos[3]+label_yshift
            xyouts, tx,ty, normal=1, fig_labels[pid]
        endfor
        if n_elements(bar_times) ne 0 then timebar, bar_times, color=sgcolor('red'), linestyle=1
        
        tpos = poss[*,0]
        msg = strupcase('RBSP-'+probe)
        tx = total(tpos[[0,2]])*0.5
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,normal=1, msg, color=sc_color, alignment=0.5
        

        ; Add FMLat.
        keo_var = prefix+'keo'
        pid = where(plot_vars eq keo_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            get_data, keo_var, limits=lim
            xrange = time_range
            yrange = lim.yrange
            plot, xrange, yrange, $
                xstyle=5, xrange=xrange, $
                ystyle=5, yrange=yrange, $
                position=tpos, nodata=1, noerase=1
            
            fmlat_var = prefix+'fmlat_'+internal_model+'_north'
            get_data, fmlat_var, times, fmlats
            fmlats = fmlats[*,model_index]
            oplot, times, fmlats, color=sc_color, linestyle=2
            ; Add sc.
            tx = time_range[0]
            ty = interpol(fmlats,times,tx)
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]+xchsz*0.5
            ty = tmp[1]-ychsz*1
            sc_name = sc_info['sc_name']
            msg = strupcase(sc_name+'-'+probe)
            xyouts, tx,ty,normal=1, msg, charsize=label_size, alignment=0, color=sc_color
        endif

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
        tmp = tmp[*,where(lim.labels eq 'N EMFISIS')]
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
        tx = tpos[0]-xchsz*8
        ty = tpos[3]+label_yshift
        msg = rbsp_letter+'-3) E/B'
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
        pf_mor_info_var = pf_mor_var+'_fft_info'
        get_data, pf_mor_info_var, tmp, info
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
        xyouts, tx,ty, msg, normal=1, charsize=label_size, alignment=0
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = rbsp_letter+'-4) Power'
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
            pf_mor_info_var = pf_mor_var+'_fft_info'
            get_data, pf_mor_info_var, tmp, info
            psd2 = 2*info.c_tau*info.dt/info.cdelta*pf_gws*1e3
            fs2 = info.fs*1e3*0.5

            xxs = psd2
            yys = fs2
            
            oplot, xxs, yys, color=sgcolor('light_salmon')
        endif

    endforeach

end




function fig_2015_0416_0800_rbsp_v01_mid_panel, my_pos, event_info=event_info

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)


    time_range = time_double(['2015-04-16/08:00','2015-04-16/08:25'])
    dp_times = event_info['dp_times']
    
    label_size = 0.8
    tmp = smkarthm(0,2*!dpi,20,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1


;---Plot ewo and dp.
    ; ewogram.
    ewo_var = 'thg_asf_ewo'
    mlat_range = [62,64]
    mlt_range = [-4,-2]
    ewo_zlog = 0
    if ewo_zlog eq 0 then begin
        zrange = [-1,1]*5e3
        ct = 70
        ztickv = [-1,0,1]*abs(max(zrange))
        zticks = n_elements(ztickv)-1
        zminor = 10
        ztickn = string(ztickv,format='(I0)')
    endif else begin
        zrange = [1e1,1e4]
        ct = 49
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,1,inner=1)
        ztickv = 10^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        ztickn = '10!U'+string(log_ztickv,format='(I0)')
        index = where(ztickn eq '10!U1', count)
        if count ne 0 then ztickn[index] = '10'
        index = where(ztickn eq '10!U0', count)
        if count ne 0 then ztickn[index] = '1'
    endelse
    
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
    ewo = total(data[*,mlt_index,mlat_index],3)/nmlat_bin
    if ewo_zlog eq 1 then ewo = ewo>0.01
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
        'zrange', zrange, $
        'zlog', ewo_zlog, $
        'ztickv', ztickv, $
        'zticks', zticks, $
        'zminor', zminor, $
        'ztickname', ztickn, $
        'color_table', ct )


    foreach info, event_info.rbsp do begin
        prefix = info['prefix']
        probe = info['probe']

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
    

;---Plot settings.
    plot_info = orderedhash()
    plot_info[ewo_var] = dictionary($
        'fig_label', 'Ewo', $
        'ypan', 1.5 )
    plot_info[tilt_var] = dictionary($
        'fig_label', 'B tilt' )

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
    fig_labels = 'b-'+string(findgen(nvar)+1,format='(I0)')+') '+fig_labels


    poss = sgcalcpos(nvar, ypans=ypans, ypad=ypads, position=my_pos)
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]
    
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
    
    tplot_options, 'version', 2
    tplot, plot_vars, position=poss, trange=time_range, noerase=1, single_line_uttick=1
    label_yshift = -ychsz*0.7
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*8
        ty = tpos[3]+label_yshift
        xyouts, tx,ty, normal=1, fig_labels[pid]
    endfor


    ; Add FMLT.
    ewo_var = 'thg_asf_ewo'
    pid = where(plot_vars eq ewo_var, count)
    
    rbsp_info = event_info.rbsp.rbspa
    model_setting = rbsp_info['model_setting']
    all_models = model_setting['models']
    internal_model = event_info.internal_model
    external_model = event_info.external_model
    model_index = where(all_models eq external_model)
    
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
            
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            txs = tmp[0]+[0,0]
            tys = [tmp[1],poss[1,-1]]
            plots, txs,tys, normal=1, linestyle=1, color=sc_color
        endforeach
        
        ; overall slope.
        plots, dp_times,dp_mlts, linestyle=1
        dp_omega = total(dp_mlts*[-1,1])/total(dp_times*[-1,1])*60*15
        event_info['dp_omega'] = dp_omega
        tx = mean(dp_times)
        ty = mean(dp_mlts)
        tmp = convert_coord(tx, ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*2
        ty = tmp[1]+ychsz*0.2
        msg = tex2str('omega')+' = '+string(abs(dp_omega),format='(F3.1)')+' deg/min'
        xyouts, tx,ty,normal=1, alignment=0, msg, charsize=label_size
    endif


    return, poss
end

function fig_2015_0416_0800_rbsp_v01_asi_panel, asi_poss, event_info=event_info

    label_size = 0.8
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]

    asi_times = event_info['dp_times']
    nasi_panel = n_elements(asi_times)
    asi_letters = 'a-'+string(findgen(nasi_panel)+1,format='(I0)')
    asi_setting = (event_info.ground)['asi_setting']
    asi_sites = sort_uniq(asi_setting.sites)


    ; Settings.
    mlt_range = [-1d,0]*6
    min_mlat = 50d
    mlat_range = [min_mlat,90]
    ct_asi = 70
    asi_zlog = 0
    color_top = 254
    nangle = 50
    dangle = !dpi/nangle
    angles = make_bins([-!dpi,0], dangle)


    ; xlabels.
    xtickn_pos = (90.-min_mlat)/(90.-min_mlat+6)
    xrange = mlt_range
    xstep = 6d
    xtickv = smkarthm(xrange[0]+1,xrange[1]-1, 2, 'dx')
    xtick_minor = make_bins(xrange, 1, inner=1)
    xtickn = xtickv
    index = where(xtickv lt 0, count)
    if count ne 0 then xtickn[index] += 24
    xtickn = string(xtickn,format='(I02)')
    
    ; ylabels.
    ytick_pos = (-1+1d/12)*!dpi
    ytick_pos = (-1)*!dpi
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange, 10, inner=1)
    ytick_minor = make_bins(yrange, ystep, inner=1)
    ytickn = string(ytickv,format='(I0)')
    
    
    ; ASI.
    asi_zrange_list = list()
    asi_zrange_list.add, [-1,1]*3e3
    asi_zrange_list.add, [-1,1]*5e3
    ;asi_zrange_list.add, [-1,1]*1e4
    ;asi_zrange_list.add, [-1,1]*1e4
    foreach asi_time, asi_times, asi_id do begin
        asi_mlt_image_var = 'thg_asf_mlt_image'
        mlt_image = get_var_data(asi_mlt_image_var, at=asi_time)
        npx = n_elements(mlt_image[0,*])
        mlt_image = mlt_image[0:npx*0.5-1,0:npx*0.5-1]
        
        asi_zrange = asi_zrange_list[asi_id]
        if asi_zlog eq 1 then begin
            asi_log_zrange = alog10(asi_zrange)
            asi_zzs = bytscl(alog10(mlt_image), min=asi_log_zrange[0],max=asi_log_zrange[1], top=color_top)
            asi_log_ztickv = make_bins(asi_log_zrange,1, inner=1)
            asi_ztickv = 10^asi_log_ztickv
            asi_zticks = n_elements(asi_ztickv)-1
            asi_zminor = 9
            asi_ztickn = '10!U'+string(asi_log_ztickv,format='(I0)')
        endif else begin
            ;zrange = to_power_scale(asi_zrange)
            ;asi_zzs = to_power_scale(mlt_image)
            zrange = asi_zrange
            asi_zzs = mlt_image
            asi_zzs = bytscl(asi_zzs, min=zrange[0],max=zrange[1], top=color_top)
            asi_ztickv = sort_uniq([asi_zrange,0])
            asi_zticks = n_elements(asi_ztickv)-1
            asi_zminor = 10
            asi_ztickn = string(asi_ztickv,format='(I0)')
        endelse

        the_info = dictionary($
        'msg', [asi_letters[asi_id]+') North | white light',strjoin(strupcase(asi_sites),' ')+' | '+$
        time_string(asi_time,tformat='hh:mm:ss')+' UT'], $
        'hemisphere', 'north', $
        'position', asi_poss[*,asi_id], $
        'zzs', asi_zzs, $
        'ct', ct_asi )

        tpos = the_info.position

        ; Draw data.
        zzs = the_info.zzs
        ct = the_info.ct
        sgtv, zzs, ct=ct, position=tpos

        ; Add labels, etc.
        ; Draw axes.
        plot, [-1,0], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos
        plots, [0,-1], [0,0], color=sgcolor('silver')
        plots, [0,0], [0,-1], color=sgcolor('silver')
        ;plots, [0,1], [0,0], color=sgcolor('silver')
        

        ; circles for ytickv.
        foreach yminor, ytick_minor, val_id do begin
            rr = (yminor-min_mlat)/(90-min_mlat)
            txs = rr*cos(angles)
            tys = rr*sin(angles)
            linestyle = 1
            index = where(ytickv eq yminor, count)
            if count ne 0 then linestyle = 0
            oplot, txs,tys, linestyle=linestyle, color=sgcolor('silver')
        endforeach


        ; lines for xickv.
        foreach xminor, xtick_minor, val_id do begin
            linestyle = 1
            index = where(xtickv eq xminor, count)
            if count ne 0 then linestyle = 0
            
            tt = (xminor*15-90)*constant('rad')
            txs = [0,1]*cos(tt)
            tys = [0,1]*sin(tt)

            plots, txs,tys, data=1, linestyle=linestyle, color=sgcolor('silver')
        endforeach
        
        ; add yticknames.
        foreach yminor, ytickv, val_id do begin
            rr = 1-(yminor-min_mlat)/(90-min_mlat)
            tt = ytick_pos
            tx = rr*cos(tt)
            ty = rr*sin(tt)
            msg = ytickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*label_size*0.35
            xyouts, tx,ty,normal=1, alignment=0.5, msg, charsize=label_size
        endforeach

        ; add xticknames.
        foreach xminor, xtickv, val_id do begin
            tmp = (xminor*15-90)*constant('rad')
            rr = xtickn_pos
            tx = rr*cos(tmp)
            ty = rr*sin(tmp)
            msg = xtickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            xyouts, tx,ty-ychsz*0.3,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
        endforeach
        
        ; Add panel label.
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.2
        msgs = the_info.msg
        xyouts, tx,ty,normal=1, msgs[1], charsize=label_size
        xyouts, tx,ty+ychsz*1,normal=1, msgs[0]

        ; Colorbar.
        asi_cbpos = asi_poss[*,asi_id]
        horizontal = 1
        if horizontal eq 1 then begin
            asi_cbpos[1] = asi_cbpos[3]+ychsz*0.4
            asi_cbpos[3] = asi_cbpos[1]+ychsz*0.4
            asi_cbpos[0] += xchsz*1
            asi_cbpos[2] -= xchsz*1
        endif else begin
            asi_cbpos[0] = asi_poss[2,1]+xchsz*1
            asi_cbpos[2] = asi_cbpos[0]+xchsz*1
        endelse
        asi_ztitle = 'N-hem ASI (#)'
        asi_linestyle = 1
        zticklen = uniform_ticklen/(asi_cbpos[3]-asi_cbpos[1])/fig_size[1]
        sgcolorbar, findgen(color_top), horizontal=horizontal, $
            ztitle=asi_ztitle, zrange=asi_zrange, ct=ct_asi, position=asi_cbpos, $
            ztickv=asi_ztickv, zticks=asi_zticks, ztickname=asi_ztickn, zminor=asi_zminor, zticklen=zticklen, log=asi_zlog
    endforeach

;---SC.
    rbsp_info = event_info.rbsp.rbspa
    model_setting = rbsp_info['model_setting']
    internal_model = event_info['internal_model']
    external_model = event_info['external_model']
    ;external_model = 't04s'
    models = model_setting.models
    model_index = where(models eq external_model)
    foreach asi_time, asi_times, asi_id do begin
        tpos = asi_poss[*,asi_id]
        plot, [-1,0], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos

        foreach sc_info, event_info.rbsp do begin
            prefix = sc_info['prefix']
            probe = sc_info['probe']
            sc_name = sc_info['sc_name']
            sc_color = sc_info['sc_color']

            suffix = '_'+internal_model+'_'+the_info.hemisphere
            fmlts = get_var_data(prefix+'fmlt'+suffix, at=asi_time)+24
            fmlats = get_var_data(prefix+'fmlat'+suffix, at=asi_time)
            fmlt = fmlts[model_index]
            fmlat = abs(fmlats[model_index])

            tr = (90-fmlat)/(90-min_mlat)
            tt = (fmlt*15-90)*!dtor
            tx = tr*cos(tt)
            ty = tr*sin(tt)
            tmp = convert_coord(tx, ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            plots, tx,ty,normal=1, psym=6, symsize=label_size, color=sc_color
            msg = strupcase(sc_name)+'-'+strupcase(probe)
            if the_info.hemisphere eq 'south' then msg = strupcase(probe)
            xyouts, tx-xchsz*0.5,ty-ychsz*0.5, alignment=1,normal=1, $
                msg, color=sc_color, charsize=sc_label_size
        endforeach
    endforeach

    

end




function fig_2015_0416_0800_rbsp_v01, event_info=event_info, test=test

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
    event_info['dp_times'] = time_double('2015-04-16/'+['08:07:42','08:10:00'])

;---Plot file.
    test_plot_panels = 0
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,$
        strlowcase('fig_'+event_info.id+'_rbsp_'+version+'.pdf')])
    if keyword_set(test) then plot_file = 0


;---Figure out panel size.
    tmp = get_abs_chsz()
    abs_xchsz = tmp[0]
    abs_ychsz = tmp[1]
    uniform_ticklen = -abs_ychsz*0.15

    asi_times = event_info['dp_times']
    nasi_panel = n_elements(asi_times)
    asi_xpan_size = 1.8
    asi_ypan_size = asi_xpan_size
    mid_xpan_size = 2.2
    top_xpads = [10d,8]
    top_margins = [1,2,1,3]
    top_ysize = asi_ypan_size+total(top_margins[[1,3]])*abs_ychsz
    top_xsize = asi_xpan_size*2+mid_xpan_size+(total(top_xpads)+total(top_margins[[0,2]]))*abs_xchsz

    margins = [0,0,0,0]
    ypad = 0
    bottom_ysize = 3.5

    all_poss = panel_pos(plot_file, ypans=[top_ysize,bottom_ysize], $
        margins=margins, panid=[0,0], pansize=[top_xsize,top_ysize], ypad=ypad, fig_size=fig_size)
    
    top_pos = all_poss[*,0]
    top_poss = sgcalcpos(1,nasi_panel+1, xpad=top_xpads, $
        xpans=[asi_xpan_size,mid_xpan_size,asi_xpan_size], margins=top_margins, region=top_pos, $
        xchsz=xchsz, ychsz=ychsz)
    asi_poss = top_poss[*,[0,2]]
    mid_pos = top_poss[*,1]+[0,1,0,1]*ychsz*(top_margins[3]-1)

    bottom_pos = all_poss[*,1]
    bottom_margins = [11,4,8,0]
    bottom_poss = sgcalcpos(1,2, xpad=14, margins=bottom_margins, region=bottom_pos)
    left_pos = bottom_poss[*,0]
    right_pos = bottom_poss[*,1]
    rbsp_poss = [[left_pos],[right_pos]]
 
    
;---Test panels.
    if keyword_set(test_plot_panels) then begin
        panel_list = list()
        panel_list.add, mid_pos
        panel_list.add, left_pos
        panel_list.add, right_pos
        for ii=0,nasi_panel-1 do panel_list.add, reform(asi_poss[*,ii])

        sgopen, 0, size=fig_size

        foreach tpos, panel_list do begin
            xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            plot, [0,1],[0,1], $
                xstyle=1, ystyle=1, $
                xtickformat='(A1)', ytickformat='(A1)', $
                xticklen=xticklen, yticklen=yticklen, $
                nodata=1, noerase=1, position=tpos
        endforeach
        stop
    endif
    
    

;--Panels.
    sgopen, plot_file, size=fig_size
    tpos = fig_2015_0416_0800_rbsp_v01_asi_panel(asi_poss, event_info=event_info)
    tpos = fig_2015_0416_0800_rbsp_v01_mid_panel(mid_pos, event_info=event_info)
    tpos = fig_2015_0416_0800_rbsp_v01_rbsp_panel(rbsp_poss, event_info=event_info)

    if keyword_set(test) then stop
    sgclose
    
    return, plot_file

end

test = 0
print, fig_2015_0416_0800_rbsp_v01(event_info=event_info, test=test)
end