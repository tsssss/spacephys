;+
; Trace O+ twice.
;-



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
    

    xticklen_chsz = -0.2
    yticklen_chsz = -0.35
    psym = -1
    symsz = 0.5


;---Trace O+ data.
    model_time = time_double('2013-05-01/07:38:03')
    beam_dis = 2.1
    beam_dis = 3.0
    color_beam = sgcolor('red')
    color_conics = sgcolor(['deep_pink','purple'])
    test_pitch_angles = [162d,144]
    color_conics = sgcolor(['deep_pink'])
    test_pitch_angles = [162d]

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
        index = where_pro(times, '[]', trs)
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
            the_var = prefix+species+'_en_spec_'+dir
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


    ; Plot data.
    spec_vars = prefix+[$
        ['p_en_spec','p_en_spec_'+['para','anti']], $
        ['o_en_spec','o_en_spec_'+['para','anti']] ]
    nspec_var = n_elements(spec_vars)
    
    
    yrange = [5,5e4]
    foreach var, spec_vars do begin
        options, var, 'yrange', yrange
    endforeach

    margins = [15,4,10,2]
    spec_poss = sgcalcpos(nspec_var, margins=margins)
    sgopen, 0, size=[8,8], xchsz=xchsz, ychsz=ychsz
    tplot, spec_vars, trange=time_range, position=spec_poss
    timebar, snapshot_time, color=sgcolor('red')
    
    ; Plot plot_info.
    test_color = sgcolor('red')
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
    
        
        ; Add test times.
        test_times = (the_output['test_times']).toarray()
        ntest_time = n_elements(test_times)
        if ntest_time eq 0 then continue
        plots, test_times, energys, data=1, color=test_color, psym=psym, symsize=symsz
    
        beam_times_2 = (the_output['beam_times_2']).toarray()
        plots, beam_times_2, energys, data=1, color=test_color, psym=psym, symsize=symsz
        
        beam_times_3 = (the_output['beam_times_3']).toarray()
        plots, beam_times_3, energys, data=1, color=test_color, psym=psym, symsize=symsz
        
        
        conic_times_2 = (the_output['conic_times_2']).toarray()
        foreach color_conic, color_conics, conic_id do begin
            plots, conic_times_2[*,conic_id], energys, data=1, color=color_conic, psym=psym, symsize=symsz
        endforeach

        conic_times_3 = (the_output['conic_times_3']).toarray()
        foreach color_conic, color_conics, conic_id do begin
            plots, conic_times_3[*,conic_id], energys, data=1, color=color_conic, psym=psym, symsize=symsz
        endforeach
    endforeach
    
    
    var = prefix+'o_en_spec_para'
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
    
    ty = tpos[1]+(0.2+2)*ychsz
    msg = 'Beam, '+string(beam_dis-1,format='(F3.1)')+' Re'
    xyouts, tx,ty,normal=1, msg, charsize=label_size, color=color_beam



end