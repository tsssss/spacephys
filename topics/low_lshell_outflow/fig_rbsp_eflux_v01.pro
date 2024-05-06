
function fig_rbsp_eflux_v01, plot_dir, event_info=event_info, test=test

    version = 'v01'
    id = '2015_0317'
    if n_elements(event_info) eq 0 then event_info = low_lshell_outflow_load_data(id)
    label_size = 0.8
    
    time_range = time_double(['2015-03-17','2015-03-18'])
    time_range = time_double(['2015-03-17/06:00','2015-03-17/15:00'])
    rbsp_probes = ['a','b']
    rbsp_probes = ['a']
    plasma_cloak_times = dictionary($
        'a', time_double([$
            '2015-03-17/13:46','2015-03-17/14:06',$
            '2015-03-17/15:48','2015-03-17/16:34',$
            '2015-03-17/22:44','2015-03-17/23:18']), $
        'b', time_double([$
            '2015-03-17/09:44','2015-03-17/10:24',$
            '2015-03-17/12:34','2015-03-17/12:58',$
            '2015-03-17/19:28','2015-03-17/19:45',$
            '2015-03-17/21:16','2015-03-17/21:46']) )

    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'

        p_en_var = prefix+'p_en_spec'
        o_en_var = rbsp_read_en_spec(time_range, probe=probe, species='o')
        o_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='o', errmsg=errmsg, energy_range=[50,5e4], update=1)

        the_vars = [p_en_var,o_en_var]
        options, the_vars, 'ytitle', 'Energy!C(eV)'
        the_vars = [p_en_var,o_en_var]
        options, the_vars, ytickname=['10!U2','10!U3','10!U4'], $
            ytickv=[100,1000,10000], yticks=2, yminor=9, $
            yrange=[20,3e4]

        the_vars = o_pa_var
        options, the_vars, 'ytitle', 'PA!C(deg)'

        zticklen = -0.5
        the_vars = [p_en_var,o_en_var,o_pa_var]
        options, the_vars, zticklen=zticklen, zminor=9
        
        
        ; test particles.
        if probe eq 'a' then begin
            target_time_range = time_double(['2015-03-17/07:00','2015-03-17/15:00'])
        endif else begin
            target_time_range = time_double(['2015-03-17/12:00','2015-03-17/20:00'])
            ;target_time_range = time_double(['2015-03-17/03:00','2015-03-17/20:00'])
            target_time_range = time_double(['2015-03-17/03:00','2015-03-17/11:00'])
        endelse
        target_energy_range = [20,3e4]
        spec_var = o_en_var
        get_data, spec_var, times, fluxs, energys
        time_index = where_pro(times, '[]', target_time_range)
        target_times = times[time_index]
        energy_index = where_pro(energys[0,*], '[]', target_energy_range)
        the_energys = energys[0,energy_index]

        ntarget = n_elements(target_times)
        target_energys = fltarr(ntarget)
        target_fluxs = fltarr(ntarget)
        foreach tid, time_index, target_id do begin
            the_fluxs = reform(fluxs[tid,energy_index])
            the_fluxs = smooth(the_fluxs,3, edge_zero=1)
            target_fluxs[target_id] = max(the_fluxs, the_energy_index)
            target_energys[target_id] = the_energys[the_energy_index]
        endforeach
        var = spec_var+'_test_particle'
        
        lshells = get_var_data(prefix+'lshell', at=target_times)
        index = where(target_fluxs le 1e6 or lshells le 3, count)
        if count ne 0 then begin
            target_times[index] = !values.f_nan
            target_energys[index] = !values.f_nan
        endif
        store_data, var, target_times, target_energys, target_fluxs
        
        
        ; nflux.
        pa_range = [0d,90]
        en_range = target_energy_range
        en_range = [20,1e4]
        
        ; Read L3 spec, with PA and EN info.
        species = 'o'
        files = rbsp_load_hope(time_range, id='l3%pa', probe=probe, errmsg=errmsg)

        var_list = list()
        species_suffix = (species eq 'e')? '_Ele': '_Ion'
        time_var = 'Epoch'+species_suffix
        energy_var = 'HOPE_ENERGY'+species_suffix
        flux_var = strupcase('f'+species+'du')
        var_list.add, dictionary($
            'in_vars', [energy_var,flux_var], $
            'time_var_name', time_var, $
            'time_var_type', 'Epoch' )
        read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
        pitch_angles = cdf_read_var('PITCH_ANGLE', filename=files[0])
        npitch_angle = n_elements(pitch_angles)
        dpas = pitch_angles-shift(pitch_angles,1)
        dpas = [pitch_angles[0],dpas[2:npitch_angle-2],pitch_angles[0]]
        
        
        ;nflux = flux * (2*!dpi*sin(pa)*dpa) * de * cos(pa), flux in #/cm^2-s-sr-keV, nflux in #/cm^2-s; de in keV, pa in rad.
        fluxs = get_var_data(flux_var, in=time_range, times=times, limits=lim)
        index = where(fluxs le 1e-30, count)
        if count ne 0 then fluxs[index] = 0d
        energys = get_var_data(energy_var)
        nenergy = n_elements(energys[0,*])
        ntime = n_elements(times)
        para_fluxs = fltarr(ntime)
        foreach time, times, time_id do begin
            the_energys = reform(energys[time_id,*])
            des = the_energys-shift(the_energys,1)
            de_es = des/the_energys
            de_e = mean(de_es[1:-2])
            
            ; in [energy,pitch_angle], #/cm^2-s-sr-keV
            the_fluxs = reform(fluxs[time_id,*,*])
            for ii=0,npitch_angle-1 do begin
                ; convert to keV/cm^2-s-sr-keV
                the_fluxs[*,ii] *= the_energys*1e-3
            endfor
            ; convert to #/cm^2-s-sr
            the_fluxs *= de_e

            energy_index = where_pro(the_energys, '[]', en_range, count=count)
            if count ne 0 then begin
                the_energys = the_energys[energy_index]
                the_fluxs = the_fluxs[energy_index,*]
            endif
            pa_index = where_pro(pitch_angles, '[]', pa_range, count=count)
            if count ne 0 then begin
                the_pas = pitch_angles[pa_index]*constant('rad')
                the_dpas = dpas[pa_index]*constant('rad')
                the_fluxs = the_fluxs[*,pa_index]
            endif

            ; integrate in energy and flux, in #/cm^2-s
            unit = '#/cm!U2!N-s'
            foreach the_pa, the_pas, pa_id do begin
                para_fluxs[time_id] += total(the_fluxs[*,pa_id]*2*!dpi*sin(the_pa)*dpas[pa_id]*cos(the_pa),nan=1)
            endforeach
        endforeach
        
        flux_var = prefix+'o_nflux'
        cmap_var = prefix+'cmap'
        cmap = get_var_data(cmap_var, at=times)
        norm_fluxs = para_fluxs*cmap
        store_data, flux_var, times, norm_fluxs
        
        yrange = [1e5,1e12]
        ;yrange = minmax(norm_fluxs)
        log_yrange = alog10(yrange)
        log_ytickv = make_bins(log_yrange,1,inner=1)
        ytickv = 10d^log_ytickv
        yticks = n_elements(ytickv)-1
        yminor = 9
        ylog = 1
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        ytickn[0:*:2] = ' '
        ytickn[-1] = ' '
        
        
        add_setting, flux_var, smart=1, dictionary($
            'ylog', 1, $
            'yrange', yrange, $
            'ytickv', ytickv, $
            'yticks', yticks, $
            'yminor', yminor, $
            'ytickname', ytickn, $
            'constant', ytickv, $
            'display_type', 'scalar', $
            'short_name', 'O+ Nflux', $
            'unit', unit )

        
        
        ; eflux.
        theta_pa = 30d
        sr = !dpi*sin(theta_pa*constant('rad'))^2
        nenergy_bin = 3     ; # of energy_bins around
        eflux = target_fluxs*target_energys^2*1e-3*sr*nenergy_bin  ; in eV/cm^2-s
        eflux *= 1.6e-19*1e4*1e3  ; convert from eV/cm^2-s to mW/m^2.
        cmap_var = prefix+'cmap'
        cmap = get_var_data(cmap_var, at=target_times)
        var = prefix+'o_eflux_map'
        store_data, var, target_times, eflux*cmap
        yrange = [0.05,50]
        add_setting, var, smart=1, dictionary($
            'yrange', yrange, $
            'ylog', 1, $
            'display_type', 'scalar', $
            'short_name', 'O+ KEflux', $
            'unit', 'mW/m!U2!N' )
        
        
        yrange = [0.05,50]
        pflux = get_var_data(prefix+'pf_fac', times=times)
        cmap = get_var_data(cmap_var, at=times)
        mlats = get_var_data(prefix+'mlat', at=times)
        sign = (mlats ge 0)*2-1
        var = prefix+'pf_fac_map'
        store_data, var, times, -pflux[*,0]*cmap*sign
        add_setting, var, smart=1, dictionary($
            'yrange', yrange, $
            'ylog', 1, $
            'display_type', 'scalar', $
            'short_name', 'S!DEarthward!N', $
            'unit', 'mW/m!U2!N' )
        
        vars = prefix+['pf_fac_map','o_eflux_map']
        options, vars, 'constant', [0.1,1,10]
        
        
        var = prefix+'eflux_combo'
        yrange = [2e-3,50]
        log_yrange = alog10(yrange)
        log_ytickv = fix(make_bins(log_yrange, 1, inner=1))
        ytickv = 10d^log_ytickv
        yticks = n_elements(ytickv)-1
        yminor = 9
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        index = where(ytickv eq 10, count)
        if count ne 0 then ytickn[index] = '10'
        index = where(ytickv eq 1, count)
        if count ne 0 then ytickn[index] = '1'        
        index = where(log_ytickv eq -1, count)
        if count ne 0 then ytickn[index] = '0.1'
        ytickn = snum2str(ytickv,shortfloat=1)
        ;ytickn[0] = '0.001'
        store_data, var, data=prefix+['o_eflux_map','pf_fac_map'], limits={$
            colors: sgcolor(['black','red']), labels: ['O+,Outward','S,Earthwd'], labflag:-1, $
            yrange:yrange, ylog:1, ytitle: '(mW/m!U2!N)', constant:[0.01,0.1,1,10], $
            yticks:yticks, ytickv:ytickv, yminor:yminor, ytickname:ytickn }

        ; dis.
        lshell_var = prefix+'lshell'
        var = lshell_var
        options, var, 'yrange', [1,7]
        options, var, 'ytickv', [2,4,6]
        options, var, 'yticks', 2
        options, var, 'yminor', 2
        options, var, 'constant', [2,4,6]
        options, var, 'ytitle', '(#)'
        options, var, 'labels', 'L-shell'
    endforeach

    plot_vars = []
    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'
        plot_vars = [plot_vars, prefix+['o_en_spec','o_pa_spec','eflux_combo','o_nflux','lshell']]
    endforeach
    nvar = n_elements(plot_vars)
    fig_letters = letters(nvar)
    fig_labels = []
    ypans = []
    foreach probe, rbsp_probes do begin
        the_labels = ['O+ EN','O+ PA','Eflux','Nflux','L-shell']
        fig_labels = [fig_labels,fig_letters+') '+the_labels]
        tmp = fltarr(n_elements(the_labels))+1
        index = where(plot_vars eq prefix+'lshell', count)
        if count ne 0 then tmp[index] = 0.6
        index = where(plot_vars eq prefix+'mlat', count)
        if count ne 0 then tmp[index] = 0.6
        ypans = [ypans,tmp]
    endforeach


    pansize = [6,0.7]
    margins = [10,3,8,1]
    
    if n_elements(plot_dir) eq 0 then plot_dir = event_info.plot_dir
    plot_file = join_path([plot_dir,'fig_rbsp_eflux_'+id+'_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    poss = panel_pos(plot_file, ypans=ypans, fig_size=fig_size, $
        nypan=nvar, pansize=pansize, panid=[0,1.2], margins=margins)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz, inch=1

    xticklen_chsz = -0.25   ; in ychsz.
    yticklen_chsz = -0.40   ; in xchsz.
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        options, plot_vars[ii], 'xticklen', xticklen
        options, plot_vars[ii], 'yticklen', yticklen
    endfor

    tplot_options, 'tickinterval', 3600*1
    tplot_options, 'version', 3
    var_labels = ''
    tplot, plot_vars, trange=time_range, position=poss
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*(margins[0]-0.5)
        ty = tpos[3]-ychsz*0.7
        msg = fig_labels[pid]
        xyouts, tx,ty,normal=1, msg
    endfor
    
    
    ; Add RBSP.
    pid = 0
    tpos = poss[*,pid]
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = strupcase('rbsp-'+probe)
    xyouts, tx,ty,normal=1, msg, color=sgcolor('white')


    ; Add PP times.
    pp_times = ((event_info['rbsp'])['rbspa'])['pp_times']
    pp_times = pp_times[where_pro(pp_times,'[]',time_range)]
    timebar, pp_times, linestyle=2, color=sgcolor('red')

    plasma_cloak_color = sgcolor('red')
    void_color = sgcolor('light_blue')
    thick = (keyword_set(test))? 2: 8
    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'
        foreach pa_var, prefix+['o_pa_spec','o_en_spec'] do begin
            pid = where(plot_vars eq pa_var, count)
            if count eq 0 then continue
            tpos = poss[*,pid]
            xrange = time_range
            yrange = [0,1]
            plot, xrange, yrange, $
                xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                position=tpos, nodata=1, noerase=1
            the_times = plasma_cloak_times[probe]
            nsector = n_elements(the_times)*0.5
            nsector = 1

            tys = mean(yrange)+[0,0]
            for ii=0,nsector-1 do begin
                txs = the_times[ii*2:ii*2+1]
                plots, txs, tys, data=1, color=plasma_cloak_color
                
                if pa_var eq prefix+'o_pa_spec' then begin
                    tx = mean(txs)
                    ty = tys[0]
                    tmp = convert_coord(tx,ty, data=1, to_normal=1)
                    tx = tmp[0]
                    ty = tmp[1]+ychsz*0.2
                    msg = 'PC'+string(ii+1,format='(I0)')
                    xyouts,tx,ty,normal=1, msg, alignment=0.5, color=plasma_cloak_color, charsize=label_size
                endif
                
                foreach tx, txs do begin
                    tmp = convert_coord(tx,tys[0],data=1,to_normal=1)
                    tx = tmp[0]
                    ty = tmp[1]
                    plots, tx+[0,0], ty+[-1,1]*ychsz*0.15, normal=1, color=plasma_cloak_color
                endforeach
            endfor
        endforeach


        spec_var = prefix+'o_en_spec'
        pid = where(plot_vars eq spec_var, count)
        test_particle_color = sgcolor('orange')
        if count ne 0 then begin
            tpos = poss[*,pid]
            xrange = time_range
            yrange = get_var_setting(spec_var, 'yrange')
            ylog = 1
            plot, xrange, yrange, $
                xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                position=tpos, nodata=1, noerase=1, ylog=1
            get_data, spec_var+'_test_particle', target_times, target_energys, target_fluxs
            plots, target_times, target_energys, psym=1, symsize=0.05, color=test_particle_color
        endif

        
        lshell_var = prefix+'lshell'
        pid = where(plot_vars eq lshell_var, count)
        if count eq 0 then continue
        tpos = poss[*,pid]
        xrange = time_range
        yrange = get_var_setting(lshell_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        the_times = plasma_cloak_times[probe]
        nsector = n_elements(the_times)*0.5
        nsector = 1
        for ii=0,nsector-1 do begin
            tr = the_times[ii*2:ii*2+1]
            tys = get_var_data(lshell_var, in=tr, times=txs)
            plots, txs,tys, color=plasma_cloak_color, thick=thick
        endfor
        
        
        ; Get the times for L in given range.
        lshell_range = [4,5]
        lshells = get_var_data(lshell_var, times=times)
        index = where_pro(lshells, '[]', lshell_range)
        trs = times[time_to_range(index,time_step=1)]
        for ii=1,1 do begin
            tr = reform(trs[ii,*])
            
            pid = where(plot_vars eq lshell_var, count)
            if count eq 0 then continue
            tpos = poss[*,pid]
            xrange = time_range
            yrange = get_var_setting(lshell_var, 'yrange')
            plot, xrange, yrange, $
                xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                position=tpos, nodata=1, noerase=1
            tys = get_var_data(lshell_var, in=tr, times=txs)
            oplot, txs,tys, color=void_color, thick=thick
            
            foreach pa_var, prefix+['o_pa_spec','o_en_spec'] do begin
                pid = where(plot_vars eq pa_var, count)
                if count eq 0 then continue
                tpos = poss[*,pid]
                xrange = time_range
                yrange = [0,1]
                plot, xrange, yrange, $
                    xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                    position=tpos, nodata=1, noerase=1
                

                tys = mean(yrange)+[0,0]
                txs = tr
                oplot, txs, tys, color=void_color
                
                if pa_var eq prefix+'o_pa_spec' then begin
                    tx = mean(txs)
                    ty = tys[0]
                    tmp = convert_coord(tx,ty, data=1, to_normal=1)
                    tx = tmp[0]
                    ty = tmp[1]+ychsz*0.2
                    msg = 'NoPC';+string(ii,format='(I0)')
                    xyouts,tx,ty,normal=1, msg, alignment=0.5, color=void_color, charsize=label_size
                endif
                
                
                foreach tx, txs do begin
                    tmp = convert_coord(tx,tys[0],data=1,to_normal=1)
                    tx = tmp[0]
                    ty = tmp[1]
                    plots, tx+[0,0], ty+[-1,1]*ychsz*0.15, normal=1, color=void_color
                endforeach
            endforeach
        endfor
        
        var = prefix+'eflux_combo'
        pid = where(plot_vars eq var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            xrange = time_range
            yrange = get_var_setting(var, 'yrange')
            ylog = 1
            plot, xrange, yrange, $
                xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                position=tpos, nodata=1, noerase=1, ylog=1

            msg = 'Flux norm. to 100 km by |B!Diono!N|/|B|'
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*0.8
            xyouts, tx,ty,normal=1, msg, charsize=label_size
        endif
        
        
        var = prefix+'o_nflux'
        pid = where(plot_vars eq var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            xrange = time_range
            yrange = get_var_setting(var, 'yrange')
            ylog = 1
            plot, xrange, yrange, $
                xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                position=tpos, nodata=1, noerase=1, ylog=1
                
            msg = 'Flux norm. to 100 km'
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*0.8
            xyouts, tx,ty,normal=1, msg, charsize=label_size
            
            msg = 'PA ['+strjoin(string(pa_range,format='(I0)'),',')+']'+'!U'+tex2str('circ')+'!N'
            ty = tpos[3]-ychsz*1.6
            xyouts, tx,ty,normal=1, msg, charsize=label_size

            msg = 'EN ['+strjoin(string(en_range,format='(I0)'),',')+']'+' eV'
            ty = tpos[3]-ychsz*2.4
            xyouts, tx,ty,normal=1, msg, charsize=label_size
        endif
    endforeach

    if keyword_set(test) then stop
    sgclose

    return, plot_file

end

print, fig_rbsp_eflux_v01(event_info=event_info, test=0)
end