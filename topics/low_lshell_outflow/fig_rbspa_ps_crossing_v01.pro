
function fig_rbspa_ps_crossing_v01, plot_dir, event_info=event_info, test=test

    version = 'v01'
    id = '2015_0317'
    if n_elements(event_info) eq 0 then event_info = low_lshell_outflow_load_data(id)

    time_range = time_double(['2015-03-17/03:10','2015-03-17/10:30'])
    time_range = time_double(['2015-03-17/07:10','2015-03-17/14:30'])
    kaw_times = time_double('2015-03-17/'+['08:18','08:39','08:58','09:24',$
        '11:13','11:33','12:11','12:36','13:00','13:23'])
    ps_times = !null
    fc_text_time = time_double('2015-03-17/10:30')
    
    sw_b_var = omni_read_sw_b(time_range, coord='gse')

    ;foreach probe, rbsp_probes do begin
    foreach probe, ['a'] do begin
        prefix = 'rbsp'+probe+'_'

        o_en_var = rbsp_read_en_spec(time_range, probe=probe, species='o')
        o_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='o', errmsg=errmsg, energy_range=[50,5e4], update=1)

        
        ; E field spec.
        e_spec_var = rbsp_read_efield_spec(time_range, probe=probe, update=update)
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
        internal_model = 'dipole'
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
            'constant', [20,30,40], $
            'yrange', [15,48], $
            'ytickv', [20,30,40], $
            'yticks', 2, $
            'yminor', 4 )
            
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
        
        ;vars = [e_var,dens_var,e_en_var,p_en_var,o_en_var,o_pa_var,$
        ;    e_spec_combo,b_spec_combo,b_theta_var,b_var]
        vars = [theta_combo_var,beta_var,e_en_var,p_en_var,o_en_var,o_pa_var,$
            b_spec_var,e_spec_var,mlt_var,mlat_var,dis_var]
        nvar = n_elements(vars)
        fig_letters = letters(nvar)
        fig_labels = fig_letters+') '+['B tilt angle','Beta','e- spec','H+ spec','O+ spec','O+ PA','B spec','E spec','MLT','MLat','|R|']
        ypans = fltarr(nvar)+1.
        index = where(vars eq e_spec_var, count)
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
        
        pansize = [6,0.7]
        margins = [10,3,7,1]
        
        if n_elements(plot_dir) eq 0 then plot_dir = event_info.plot_dir
        plot_file = join_path([plot_dir,'fig_rbsp'+probe+'_ps_crossing_'+id+'_'+version+'.pdf'])
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
        
        ; Add PS label
        if n_elements(ps_times) ne 0 then begin
            pid = where(vars eq theta_combo_var)
            tpos = poss[*,pid]
            tpos[1] = poss[1,pid+1]
            ps_color = sgcolor('wheat')
            xrange = time_range
            yrange = [0,1]
            plot, xrange, yrange, $
                xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                position=tpos, nodata=1, noerase=1
            polyfill, ps_times[[0,1,1,0,0]], yrange[[0,0,1,1,0]], data=1, color=ps_color
        endif
        
        
        
        tplot_options, 'tickinterval', 1800
        tplot_options, 'version', 3
        var_labels = ''
        tplot, vars, trange=time_range, position=poss, noerase=1
        for pid=0,nvar-1 do begin
            tpos = poss[*,pid]
            tx = tpos[0]-xchsz*(margins[0]-0.5)
            ty = tpos[3]-ychsz*0.7
            msg = fig_labels[pid]
            xyouts, tx,ty,normal=1, msg
        endfor
        
        
        foreach spec_var, [b_spec_var,e_spec_var] do begin
            pid = where(vars eq spec_var, count)
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
        
        pid = where(vars eq b_spec_var)
        tpos = poss[*,pid]
        kaw_color = sgcolor('red')
        timebar, kaw_times, linestyle=1, color=kaw_color
        xrange = time_range
        yrange = [0,1]
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        foreach kaw_time, kaw_times, id do begin
            tmp = convert_coord(kaw_time,yrange[1], data=1, to_normal=1)
            tx = tmp[0]
            ty = tmp[1]-ychsz*1.0
            msg = string(id+1,format='(I0)')
            xyouts, tx,ty,normal=1, msg, alignment=0.5, color=kaw_color
            
            if id eq 0 then begin
                msg = 'KAW times'
                tx = tx-xchsz*10
                xyouts, tx,ty,normal=1, msg, charsize=label_size, color=kaw_color
            endif
        endforeach
        
        
        ; Add PS label
        if n_elements(ps_times) ne 0 then begin
            pid = where(vars eq theta_combo_var)
            tpos = poss[*,pid]
            ps_color = sgcolor('wheat')
            xrange = time_range
            yrange = [0,1]
            plot, xrange, yrange, $
                xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                position=tpos, nodata=1, noerase=1
            tx = mean(ps_times)
            ty = 0.5
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]-xchsz*7
            ty = tmp[1]
            msg = 'PS entry due!Cto stretching'
            xyouts, tx,ty, normal=1, msg, alignment=0.5, charsize=label_size, color=sgcolor('goldenrod')
        endif

        
        
        if keyword_set(test) then stop
        sgclose
    endforeach
    
    return, plot_file
    
end

print, fig_rbspa_ps_crossing_v01(event_info=event_info, test=0)
end