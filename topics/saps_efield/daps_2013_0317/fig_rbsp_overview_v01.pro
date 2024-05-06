;+
;-

function fig_rbsp_overview_v01, test=test, event_info=event_info

    version = 'v01'
    id = '2013_0317'

    if n_elements(event_info) eq 0 then event_info = saps_efield_load_data(id)
    time_range = time_double(['2013-03-17','2013-03-18'])
    plot_time_range = time_double(['2013-03-17/08:30','2013-03-17/10:00'])
    plot_time_range = time_double(['2013-03-17/08:30','2013-03-17/10:00'])
    plot_time_range = time_double(['2013-03-17/05:00','2013-03-17/12:00'])
    ;plot_time_range = time_double(['2013-03-17/08:30','2013-03-17/12:00'])-3600
    probes = ['a']
    
    pp_time = time_double(['2013-03-17/09:31','2013-03-17/11:18'])
    pp_color = sgcolor('red')
    
    daps_time = time_double(['2013-03-17/08:49','2013-03-17/09:00'])
    daps_color = sgcolor('wheat')
    daps_color2 = sgcolor('black')

    saps_time = time_double(['2013-03-17/11:16','2013-03-17/11:30'])
    saps_color = sgcolor('wheat')
    saps_color2 = sgcolor('black')

    fc_text_time = plot_time_range[0]
    
    
    load_info = dictionary()
    update = 0

    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'

        o_en_var = rbsp_read_en_spec(time_range, probe=probe, species='o')
        o_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='o', errmsg=errmsg, energy_range=[50,5e4], update=update)

        
        ; E field spec.
        e_spec_var = prefix+'e_spec'
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
        

        ; Model field and tilt angle.
        external_model = 't89'
        internal_model = 'igrf'
        igrf = (internal_model eq 'igrf')? 1: 0
        suffix = '_'+internal_model+'_'+external_model
        bmod_gsm_var = prefix+'bmod_gsm'+suffix
        b_gsm_var = prefix+'b_gsm'
        foreach var, [b_gsm_var,bmod_gsm_var] do begin
            get_data, var, times, b_vec, limits=lim
            coord = strlowcase(lim.coord)
            if coord ne 'sm' then begin
                b_vec = cotran(b_vec, times, coord+'2sm')
            endif
            theta_var = var+'_theta'
            theta = atan(b_vec[*,2]/snorm(b_vec))*constant('deg')
            store_data, theta_var , times, theta
            add_setting, theta_var, smart=1, dictionary($
                'display_type', 'scalar', $
                'short_name', 'Tilt', $
                'unit', 'deg' )
            if var eq b_gsm_var then begin
                b_sm = b_vec
                b_sm_var = prefix+'b_sm'
                lim.coord = 'sm'
                store_data, b_sm_var, times, b_sm, limits=lim
            endif
        endforeach
        dtheta_var = prefix+'b_dtheta'
        b_theta = get_var_data(b_gsm_var+'_theta', times=times)
        bmod_theta = get_var_data(bmod_gsm_var+'_theta', at=times)
        dtheta = b_theta-bmod_theta
        store_data, dtheta_var, times, dtheta
        add_setting, dtheta_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'dTilt', $
            'unit', 'deg' )
        theta_combo_var = prefix+'b_theta_combo'
        
        store_data, theta_combo_var, times, [[b_theta],[bmod_theta]]
        add_setting, theta_combo_var, smart=1, dictionary($
            'display_type', 'stack', $
            'short_name', tex2str('theta'), $
            'unit', 'deg', $
            'labels', ['Obs','Model'], $
            'colors', ['red','blue'], $
            'constant', [30,40], $
            'yrange', [25,40], $
            'ytickv', [30,40], $
            'yticks', 1, $
            'yminor', 5 )
            
        ; magnetic pressure.
        pmag_var = rbsp_read_magnetic_pressure(b_var=b_gsm_var)
        options, pmag_var, 'ylog', 1
        ; thermal_pressure.
        pthermal_var = prefix+'p_thermal'
        foreach species, ['e','p','o'] do begin
            vars = rbsp_read_hope_moments(time_range, probe=probe, species=species)
        endforeach
        
        species = 'e'
        p_thermal = calc_thermal_pressure($
            get_var_data(prefix+species+'_n',times=times), $
            get_var_data(prefix+species+'_t',at=times))
        foreach species, ['p','o'] do begin
            p_thermal += calc_thermal_pressure($
                get_var_data(prefix+species+'_n',at=times), $
                get_var_data(prefix+species+'_t',at=times))
        endforeach
        store_data, pthermal_var, times, p_thermal
        add_setting, pthermal_var, smart=1, dictionary($
            'ylog', 1, $
            'display_type', 'scalar', $
            'short_name', 'P!Dth', $
            'unit', 'nPa' )
        
        
        p_thermal = get_var_data(pthermal_var, times=times)
        p_mag = get_var_data(pmag_var, at=times)
        beta = p_thermal/p_mag
        beta_var = prefix+'beta'
        store_data, beta_var, times, beta
        add_setting, beta_var, smart=1, dictionary($
            'ylog', 1, $
            'yrange', [1e-4,1e-1]*7, $
            'constant', [1e-2,1e-1], $
            'display_type', 'scalar', $
            'short_name', tex2str('beta'), $
            'unit', '#' )
        
        
        ; FMLat.
        fmlat_north_var = prefix+'fmlat'+suffix+'_north'
        fmlat_south_var = prefix+'fmlat'+suffix+'_south'
        fmlat_var = prefix+'fmlat'
        fmlat_var = stplot_merge([fmlat_north_var,fmlat_south_var], output=fmlat_var)
        get_data, fmlat_var, times, data
        store_data, fmlat_var, times, abs(data)
        add_setting, fmlat_var, smart=1, dictionary($
            'display_type', 'stack', $
            'ytitle', '(deg)', $
            'yrange', [40,70], $
            'ytickv', [50,60], $
            'yticks', 1, $
            'yminor', 2, $
            'constant', [50,60], $
            'labels', ['North','South'], $
            'colors', sgcolor(['red','blue']) )



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
        var = bmod_gsm_var+'_theta'
        options, var, 'yrange', [0,90]
        options, var, 'ytickv', [40,80]
        options, var, 'constant', [40,80]
        options, var, 'yticks', 1
        options, var, 'yminor', 4
        var = dtheta_var
        yrange = [-45,10]
        ystep = 20
        ytickv = make_bins(yrange,ystep, inner=1)
        yticks = n_elements(ytickv)-1
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'yrange', yrange
        options, var, 'yminor', 2
        options, var, 'constant', ytickv

        ; E
        e_var = prefix+'edot0_fac'
        copy_data, e_var, e_var+'_plot'
        e_var = e_var+'_plot'
        var = e_var
        options, var, 'yrange', [-1,1]*8
        options, var, 'ytickv', [-1,0,1]*5
        options, var, 'yticks', 2
        options, var, 'yminor', 5
        options, var, 'constant', [-1,0,1]*5
        fac_labels = ['||',tex2str('perp')+','+['west','out']]
        options, var, 'labels', 'E!D'+fac_labels
        get_data, e_var, times, data
        index = where(get_var_data(prefix+'lshell',at=times) le 2, count)
        if count ne 0 then begin
            data[index,*] = !values.f_nan
            store_data, e_var, times, data
        endif
        
        dens_var = prefix+'density'
        options, dens_var, 'constant', [10,100,1000]
        e_en_var = prefix+'e_en_spec'
        p_en_var = prefix+'p_en_spec'
        o_en_var = prefix+'o_en_spec'
        b_theta_var = prefix+'b_gsm_theta'
        b_var = prefix+'b_sm'
        pthermal_var = prefix+'p_thermal'
        
        
        models = ['t89','t96','t01','t04s']
        fmlat_var = prefix+'fmlat_dipole_south'
        fmlat_var = stplot_merge(prefix+'fmlat_dipole_'+models+'_south', output=fmlat_var, $
            labels=strupcase(['t89','t96','t01','t04s']), colors=sgcolor(['red','green','blue','black']))
        yrange = [-63,-49]
        ystep = 5
        yminor = ystep
        ytickv = make_bins(yrange,ystep, inner=1)
        set_ytick, fmlat_var, yrange=yrange, ytickv=ytickv, yminor=yminor
        options, fmlat_var, 'ytitle', 'FMLat!C(deg)'
        options, fmlat_var, 'constant', [-60,-55,-50]
        options, fmlat_var, 'labflag', -1
        
        
        
        ; KEO
        mlt_image_var = 'thg_asf_mlt_image_rect'
        mlt_images = get_var_data(mlt_image_var, times=times, settings=settings)
        ntime = n_elements(times)
        
        keo_var = 'thg_asf_keo'
        ytitle = 'MLat!C(deg)'

        mlat_range = [55,75]
        mlt_range = [1,5]
        mlt_bins = settings['mlt_bins']
        mlat_bins = settings['mlat_bins']
        index = where_pro(mlt_bins, '[]', mlt_range, count=nbin)
        ;bins = mlt_bins[index]
        ;keo = fltarr(ntime,nbin)
        keo = total(mlt_images[*,index,*],2)/nbin


        fmlt_var = prefix+'fmlt_igrf_t01_north'
        fmlt = get_var_data(fmlt_var, at=times)
        fmlt_del = [-2,0.5]
        fmlt_del = [-0.1,0.1]
        nbin = n_elements(mlat_bins)
        keo = fltarr(ntime,nbin)
        foreach time, times, time_id do begin
            index = where_pro(mlt_bins,'[]',fmlt[time_id]+fmlt_del,count=count)
            if count eq 0 then continue
            keo[time_id,*] = total(mlt_images[time_id,index,*],2)/count
        endforeach

        store_data, keo_var, times, keo, mlat_bins, limits={ytitle:ytitle,spec:1,zlog:0,color_table:49}
        options, keo_var, 'zrange', [0,5000]
        set_ytick, keo_var, yrange=[58,72], ytickv=[60,65,70], yminor=5
        
        
        vars = [theta_combo_var,e_en_var,p_en_var,o_en_var,dens_var,$
            e_var,e_spec_var,fmlat_var]
        vars = [e_en_var,p_en_var,o_en_var,dens_var,$
            e_var,keo_var,e_spec_combo]
        nvar = n_elements(vars)
        fig_letters = letters(nvar)
        fig_labels = fig_letters+') '+['e- spec','H+ spec','O+ spec','N','E','Keo','E spec']
        ypans = fltarr(nvar)+1.
        index = where(vars eq e_spec_combo, count)
        if count ne 0 then ypans[index] = 1.6
        index = where(vars eq theta_combo_var, count)
        if count ne 0 then ypans[index] = 0.8
        index = where(vars eq beta_var, count)
        if count ne 0 then ypans[index] = 0.8
        index = where(vars eq mlat_var, count)
        if count ne 0 then ypans[index] = 0.6
        index = where(vars eq mlt_var, count)
        if count ne 0 then ypans[index] = 0.6
        index = where(vars eq dis_var, count)
        if count ne 0 then ypans[index] = 0.6
        index = where(vars eq o_pa_var, count)
        if count ne 0 then ypans[index] = 0.8
        index = where(vars eq fmlat_var, count)
        if count ne 0 then ypans[index] = 0.9       
        
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
        
        pansize = [5,0.7]
        margins = [12,6,8,1]
        
        if n_elements(plot_dir) eq 0 then plot_dir = event_info.plot_dir
        plot_file = join_path([plot_dir,'fig_rbsp'+probe+'_overview_'+id+'_'+version+'.pdf'])
        print, plot_file
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
            options, vars[ii], 'xticklen', xticklen
            options, vars[ii], 'yticklen', yticklen
        endfor

        lshell_var = prefix+'lshell'
        var_labels = [mlt_var,lshell_var]
        options, mlt_var, 'ytitle', 'MLT (h)'
        options, lshell_var, 'ytitle', 'L (#)'
        
        pid = where(vars eq e_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            ;tpos[1] = poss[1,-1]
            yrange = [0,1]
            set_axis, position=tpos, xrange=plot_time_range, yrange=yrange
            polyfill, daps_time[[0,1,1,0,0]], yrange[[0,0,1,1,0]], color=daps_color
            polyfill, saps_time[[0,1,1,0,0]], yrange[[0,0,1,1,0]], color=saps_color
        endif
        
        
        ;margins = [12,6,10,1]
        ;poss = sgcalcpos(nvar, margins=margins, ypans=ypans)
        tplot_options, 'tickinterval', 3600
        tplot, vars, trange=plot_time_range, var_label=var_labels, position=poss, vlab_margin=10, noerase=1
        for pid=0,nvar-1 do begin
            tpos = poss[*,pid]
            tx = tpos[0]-xchsz*10
            ty = tpos[3]-ychsz*0.7
            msg = fig_labels[pid]
            xyouts, tx,ty,msg, normal=1
        endfor
        
        pid = 0
        tpos = poss[*,pid]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = strupcase('rbsp-'+probe)
        xyouts, tx,ty,msg, normal=1;, color=sgcolor('white')
        
        
        bar_times = make_bins(plot_time_range, 3600, inner=1)
        timebar, bar_times, linestyle=1, color=sgcolor('silver')
        
        
        timebar, pp_time, linestyle=2, color=pp_color
        pid = where(vars eq dens_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            yrange = [0,1]
            set_axis, vars[pid], position=tpos, xrange=plot_time_range, yrange=yrange
            tmp = convert_coord(mean(pp_time),yrange[0], data=1, to_normal=1)
            tx = tmp[0];+xchsz*0.5
            ty = tmp[1]+ychsz*1.5
            msg = 'Plasmapause'
            xyouts, tx,ty,msg, normal=1, color=pp_color, alignment=0.5
        endif
        
        
        pid = where(vars eq e_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            yrange = [0,1]
            set_axis, vars[pid], position=tpos, xrange=plot_time_range, yrange=yrange
            tmp = convert_coord(mean(daps_time),yrange[1], data=1, to_normal=1)
            tx = tmp[0]
            ty = tmp[1]-ychsz*1
            msg = 'DAPS E'
            xyouts, tx,ty,msg, normal=1, color=daps_color2, alignment=0.5
            
            tmp = convert_coord(mean(saps_time),yrange[0], data=1, to_normal=1)
            tx = tmp[0]
            ty = tmp[1]+ychsz*0.2
            msg = 'SAPS E'
            xyouts, tx,ty,msg, normal=1, color=saps_color2, alignment=0.5
        endif


        foreach spec_var, [e_spec_combo] do begin
            pid = where(vars eq spec_var, count)
            if count eq 0 then continue
            tpos = poss[*,pid]
            spec_settings = get_var_setting(spec_var)
            xrange = plot_time_range
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
                tx = fc_text_time
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
        
        
        if keyword_set(test) then stop
        sgclose
    endforeach
    
    return, plot_file

end

print, fig_rbsp_overview_v01(test=0, event_info=event_info)
end