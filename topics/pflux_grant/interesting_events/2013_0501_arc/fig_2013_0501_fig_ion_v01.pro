;+
; To plot H and O outflow and pitch angle distributions.
;-

function fig_2013_0501_fig_ion_v01, event_info=event_info

test = 0

;---Load data.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
    plasma_param = event_info['plasma_param']

;---Settings.
    prefix = event_info['prefix']
    probe = event_info['probe']
    time_range = event_info['time_range']
    time_range = event_info['time_range']+[0,60*20]
    snapshot_time = event_info['snapshot_time']
    fac_labels = event_info['fac_labels']
    bar_times = make_bins(time_range,600, inner=1)
    
    bar_thick = (keyword_set(test))? 1: 4
    label_size = 0.8

;    beam_dis = 2.1
    beam_dis = 3.0
    psym = -1
    symsz = 0.3
    test_particle_color = sgcolor('red')
    color_beam = sgcolor('red')
;    color_conics = sgcolor(['deep_pink','purple'])
;    test_pitch_angles = [162d,144]
    color_conics = sgcolor(['orange'])
    test_pitch_angles = [162d]
    
    ct_oxygen = 64
    ct_proton = 63
    zrange_proton = [1e4,1e6]
    zrange_oxygen = [1e4,1e6]
    pa_times_proton = time_double('2013-05-01/'+$
        ['07:40:08','07:40:53','07:41:16'])
    pa_times_oxygen = time_double('2013-05-01/'+$
        ['07:40:53','07:41:16','07:42:24'])
    npa_panel = n_elements(pa_times_proton)
    
;---Trace O+ data.
    model_time = time_double('2013-05-01/07:38:03')

    test_info_list = list()
    test_info_list.add, dictionary($
        'species', 'o', $
        'trs', time_double('2013-05-01/'+['07:40:53','07:43:32']), $
        'ens', [6200,900] )
    test_info_list.add, dictionary($
        'species', 'p', $
        'trs', time_double('2013-05-01/'+['07:40:08','07:42:01']), $
        'ens', [7000,200] )


    trace_input_list = list()
    foreach info, test_info_list do begin
        species = info['species']
        trs = info['trs']
        ens = info['ens']
        
        test_time_range = trs
        
    ;---To get the trace info.
        get_data, prefix+species+'_en_spec', times, data
        index = lazy_where(times, '[]', trs)
        test_times = times[index]
        log_ens = alog10(ens)
        foreach time, test_times, time_id do begin
            energy = 10.^(log_ens[0]+(log_ens[1]-log_ens[0])/(trs[1]-trs[0])*(time-trs[0]))
            
            trace_input_list.add, dictionary($
                'species', species, $
                'time', time, $
                'energy', energy, $
                'pitch_angles', test_pitch_angles, $
                'beam_dis', beam_dis, $
                'mod_time', model_time, $
                'model', 't89', $
                'igrf', 0 )
        endforeach
    endforeach



;---Trace and plot output.
    r_gsm_var = prefix+'r_gsm'
    trace_output_list = list()
    foreach info, trace_input_list do begin        
        time = time_double(info['time'])
        info['r_gsm'] = get_var_data(r_gsm_var, at=time)
        info['time'] = time
        par_var = info['model']+'_var'
        info['par'] = reform(get_var_data(par_var, at=time))
    
        trace_output = trace_ion_to_ionosphere(time, _extra=info.tostruct())
        trace_output['species'] = info['species']
        trace_output['energy'] = info['energy']
        trace_output_list.add, trace_output
    endforeach


    ; Collect plot_info.
    plot_info = dictionary()
    dirs = ['anti','para']
    foreach info, trace_output_list do begin
        species = info['species']
        foreach dir, dirs do begin
            the_var = prefix+species+'_en_spec_'+dir+'_plot'
            if ~plot_info.haskey(the_var) then begin
                plot_info[the_var] = dictionary($
                    'conic_times', list(), $
                    'beam_times', list(), $
                    'energys', list(), $
                    'conic_dis', list(), $
                    'beam_times_2', list(), $
                    'conic_times_2', list(), $
                    'beam_times_3', list(), $
                    'conic_times_3', list(), $
                    'beam_times_0', list(), $
                    'test_times', list() )
            endif

            key = (dir eq 'anti')? 'north': 'south'
            the_output = info[key]
            (plot_info[the_var])['conic_times'].add, the_output['conic_time']
            (plot_info[the_var])['beam_times'].add, the_output['beam_time']
            (plot_info[the_var])['energys'].add, info['energy']
            (plot_info[the_var])['conic_dis'].add, the_output['conic_dis']
            
            ; Get the info for the original test particle and those mirrored twice.
            if dir eq 'anti' then begin
                north_output = info['north']
                south_output = info['south']
                t0 = north_output['beam_time']
                t_north_sc = north_output['beam_dtime']
                t_south_sc = south_output['beam_dtime']
                t1 = t0+t_north_sc
                
                t2 = t0+t_south_sc+2*t_north_sc
                (plot_info[the_var])['beam_times_2'].add, t2
                t3 = t1+2*(t_south_sc+t_north_sc)
                (plot_info[the_var])['beam_times_3'].add, t3
                
                t0 = north_output['conic_time']
                t_north_sc = north_output['conic_dtime']
                t_south_sc = south_output['conic_dtime']
                t2 = t0+t_south_sc+2*t_north_sc
                (plot_info[the_var])['conic_times_2'].add, t2
                t3 = t1+2*(t_south_sc+t_north_sc)
                (plot_info[the_var])['conic_times_3'].add, t3
                
                ; t1 is the time of the test O+.
                (plot_info[the_var])['test_times'].add, t1
            endif
        endforeach
    endforeach


;---Plot data.
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,'fig_2013_0501_fig_ion_v01.pdf'])
    if keyword_set(test) then plot_file = 0


    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz
    pa_xpad = 0.8
    pa_size = 1.3
    xpan_size = npa_panel*pa_size+(npa_panel-1)*abs_xchsz*pa_xpad
    ypan_size = pa_size
    
    ypads = [0.4+fltarr(2),4]
    spec_pan_ratio = 0.6
    ypans = [spec_pan_ratio+fltarr(3),1]*ypan_size
    
    margins = [8,4,1.5,4]
    poss = panel_pos(plot_file, pansize=[xpan_size,ypan_size*spec_pan_ratio], $
        ypans=ypans, $
        xpads=2, ypads=ypads, margins=margins, fig_size=fig_size, $
        nxpan=2, nypan=4 )

    
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    uniform_ticklen = -ychsz*0.15*fig_size[0]
    
    pa_oxygen_poss = sgcalcpos(1,npa_panel, xpad=pa_xpad, position=poss[*,0,3])
    pa_proton_poss = sgcalcpos(1,npa_panel, xpad=pa_xpad, position=poss[*,1,3])
    
    p_vars = prefix+['p_en_spec','p_en_spec_'+['para','anti']]
    o_vars = prefix+['o_en_spec','o_en_spec_'+['para','anti']]
    orig_vars = [p_vars, o_vars]
    spec_vars = orig_vars+'_plot'
    p_vars = p_vars+'_plot'
    o_vars = o_vars+'_plot'
    foreach var, spec_vars, var_id do begin
        var = duplicate_var(orig_vars[var_id], output=var)
        get_data, var, times, data, val
        index = where(data eq 0, count)
        if count ne 0 then begin
            data[index] = 0.01
            store_data, var, times, data, val
        endif
    endforeach
    nspec_var = n_elements(spec_vars)
    
    
    yrange = [5,5e4]
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
    foreach var, spec_vars do begin
        options, var, 'yrange', yrange
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'yminor', yminor
        options, var, 'ytickname', ytickn
        options, var, 'constant', ytickv
        options, var, 'ytitle', 'Energy!C(eV)'
    endforeach
    
    
    o_spec_poss = reform(poss[*,0,0:2])
    p_spec_poss = reform(poss[*,1,0:2])
    foreach species, ['o','p'] do begin
        if species eq 'o' then begin
            vars = o_vars
            ct = ct_oxygen
            zrange = zrange_oxygen
            poss = o_spec_poss
            species_str = 'O+'
            novtitle = 0
            fig_letters = ['a','b','c']
        endif else if species eq 'p' then begin
            vars = p_vars
            ct = ct_proton
            zrange = zrange_proton
            poss = p_spec_poss
            species_str = 'H+'
            novtitle = 1
            fig_letters = ['e','f','g']
        endif
        options, vars, 'color_table', ct
        options, vars, 'zrange', zrange
        options, vars, 'no_color_scale', 1
        tplot_options, 'num_lab_min', 5

        tpos = poss[*,0]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, vars, 'xticklen', xticklen
        options, vars, 'yticklen', yticklen
        

        if species eq 'p' then begin
            options, vars, 'ytitle', ' '
            options, vars, 'ytickformat', '(A1)'
        endif
        
        tplot, vars, trange=time_range, position=poss, noerase=1, novtitle=novtitle
        timebar, snapshot_time, color=sgcolor('red')
        
        msgs = 'PA = ['+['0,180','0,45','135,180']+'] deg'
        foreach var, vars, var_id do begin
            tpos = poss[*,var_id]
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            msg = fig_letters[var_id]+')'
            xyouts, tx,ty,normal=1, msg
            tx = tpos[2]-xchsz*0.5
            ty = tpos[3]-ychsz*1
            msg = msgs[var_id]
            xyouts, tx,ty,normal=1, msg, alignment=1
        endforeach
        
        tpos = poss[*,0]
        cbpos = tpos
        cbpos[1] = tpos[3]+ychsz*0.4
        cbpos[3] = cbpos[1]+ychsz*0.5
        
        ncolor = 256
        colors = findgen(ncolor)
        ztitle = species_str+' '+get_setting(vars[0],'ztitle')
        zlog = 1
        zrange = get_setting(vars[0],'zrange')
        zticklen = uniform_ticklen/(cbpos[3]-cbpos[1])/fig_size[1]
        sgcolorbar, colors, horizontal=1, ct=ct, position=cbpos, $
            log=zlog, ztitle=ztitle, zrange=zrange, zticklen=zticklen, zminor=10, zcharsize=1
    endforeach



    
;---Plot plot_info.
    spec_poss = [[p_spec_poss],[o_spec_poss]]

    foreach var, plot_info.keys() do begin
        pid = where(spec_vars eq var)
        tpos = spec_poss[*,pid]
        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, ylog=1, $
            nodata=1, noerase=1, position=tpos
    
        the_output = plot_info[var]
        conic_times = (the_output['conic_times']).toarray()
        beam_times = (the_output['beam_times']).toarray()
        energys = (the_output['energys']).toarray()
    
        foreach color_conic, color_conics, conic_id do begin
            plots, conic_times[*,conic_id], energys, data=1, color=color_conic, psym=psym, symsize=symsz
        endforeach
    
        plots, beam_times, energys, data=1, color=color_beam, psym=psym, symsize=symsz
    
        ; Add labeling.
        center_time = mean([max(conic_times),max(beam_times)])
        ty = 200d
        tmp = convert_coord(center_time, ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        msg = 'P0'
        if strpos(var, 'para') ne -1 then msg = 'P2'
        xyouts, tx,ty,normal=1, msg, alignment=0.5, charsize=label_size

        
        ; Add test times.
        test_times = (the_output['test_times']).toarray()
        ntest_time = n_elements(test_times)
        if ntest_time eq 0 then continue
        plots, test_times, energys, data=1, color=test_particle_color, psym=psym, symsize=symsz
    
        center_time = max(test_times)
        ty = 200d
        tmp = convert_coord(center_time, ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        msg = 'P1'
        xyouts, tx,ty,normal=1, msg, alignment=0, charsize=label_size
    
        ; Add more beam tracings.
        beam_times_2 = (the_output['beam_times_2']).toarray()
        plots, beam_times_2, energys, data=1, color=test_particle_color, psym=psym, symsize=symsz
        
        beam_times_3 = (the_output['beam_times_3']).toarray()
        plots, beam_times_3, energys, data=1, color=test_particle_color, psym=psym, symsize=symsz
        
        ; Add more conic tracings.
        conic_times_2 = (the_output['conic_times_2']).toarray()
        foreach color_conic, color_conics, conic_id do begin
            plots, conic_times_2[*,conic_id], energys, data=1, color=color_conic, psym=psym, symsize=symsz
        endforeach

        conic_times_3 = (the_output['conic_times_3']).toarray()
        foreach color_conic, color_conics, conic_id do begin
            plots, conic_times_3[*,conic_id], energys, data=1, color=color_conic, psym=psym, symsize=symsz
        endforeach
        
        ; Add labels for P3 and P4.
        center_time = mean([max(conic_times_2),max(beam_times_2)])
        ty = 200d
        tmp = convert_coord(center_time, ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        msg = 'P3'
        xyouts, tx,ty,normal=1, msg, alignment=0.5, charsize=label_size

        center_time = mean([max(conic_times_3),max(beam_times_3)])
        ty = 200d
        tmp = convert_coord(center_time, ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        msg = 'P4'
        xyouts, tx,ty,normal=1, msg, alignment=0.5, charsize=label_size
    endforeach
    
    
    var = prefix+'o_en_spec_para_plot'
    pid = where(spec_vars eq var)
    tpos = spec_poss[*,pid]
    xrange = time_range
    yrange = get_setting(var, 'yrange')
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, ylog=1, $
        nodata=1, noerase=1, position=tpos
    tx = snapshot_time
    ty = yrange[0]
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]+xchsz*0.5
    tx = tpos[0]+xchsz*0.5
    the_info = plot_info[var]
    conic_dis = the_info['conic_dis'].toarray()
    foreach color_conic, color_conics, conic_id do begin
        ty = tpos[1]+(0.2+conic_id)*ychsz
        the_dis = mean(conic_dis[*,conic_id])
        msg = 'Conic, PA '+string(test_pitch_angles[conic_id],format='(I0)')+' deg, '+string(the_dis-1,format='(F3.1)')+' Re'
        xyouts, tx,ty,normal=1, msg, charsize=label_size, color=color_conic
    endforeach
    
    ty = tpos[1]+(0.2+n_elements(color_conics))*ychsz
    msg = 'Beam, '+string(beam_dis-1,format='(F3.1)')+' Re'
    xyouts, tx,ty,normal=1, msg, charsize=label_size, color=color_beam



;---Plot pitch angle distribution.
    pa_unit = 'energy'
    foreach species, ['o','p'] do begin
        if species eq 'o' then begin
            pa_times = pa_times_oxygen
            pa_poss = pa_oxygen_poss
            zrange = zrange_oxygen
            spec_pos = o_spec_poss[*,-1]
            pa_letter = 'd'
        endif else if species eq 'p' then begin
            pa_times = pa_times_proton
            pa_poss = pa_proton_poss
            zrange = zrange_proton
            spec_pos = p_spec_poss[*,-1]
            pa_letter = 'h'
        endif
        
        pa_var = rbsp_plot_pa2d_read_data(time_range, probe=probe, species=species)
        foreach time, pa_times, time_id do begin
            tpos = pa_poss[*,time_id]
            is_last_panel = time_id eq n_elements(pa_times)-1
            is_first_panel = time_id eq 0
            
            no_colorbar = 1
            ztitle = (is_last_panel)? ' ': get_setting(prefix+species+'_en_spec_plot', 'ztitle')
            ytickformat = (is_first_panel)? '': '(A1)'
            ytitle = (is_first_panel)? !null: ' '
            if species eq 'p' then begin
                ytitle = ' '
                ytickformat = '(A1)'
            endif
            tmp = rbsp_plot_pa2d(time=time, pa2d_var=pa_var, time_range, probe=probe, $
                species=species, unit=pa_unit, zrange=zrange, position=tpos, $
                no_colorbar=no_colorbar, title=' ', ztitle=ztitle, $
                ytickformat=ytickformat, ytitle=ytitle )
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            msg = pa_letter+'-'+string(time_id+1,format='(I0)')+')'
            xyouts, tx,ty, normal=1, msg
            msg = time_string(time)+' UT'
            tx = tpos[0]+xchsz*0.5
            ty = tpos[1]+ychsz*0.2
            xyouts, tx,ty,normal=1, msg, charsize=label_size
        endforeach
        
        ; Add bars for each time
        tpos = spec_pos
        xrange = time_range
        yrange = [0,1]
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach time, pa_times, time_id do begin
            tmp = convert_coord(time, yrange[0], data=1, to_normal=1)
            txs = tmp[0]+[0,0]
            tys = tmp[1]+ychsz*[-0.6,-0.4]
            plots, txs, tys, normal=1, thick=bar_thick
        endforeach
    endforeach
    
    
    if keyword_set(test) then stop
    sgclose

    return, plot_file

end


print, fig_2013_0501_fig_ion_v01(event_info=event_info)
end
