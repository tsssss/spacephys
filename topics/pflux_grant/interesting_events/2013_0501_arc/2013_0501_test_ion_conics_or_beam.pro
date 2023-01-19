;+
; Trace H+ and O+ distribution backward in time to check if they are ion conics or beams.
;-


;---Settings.
    time_range = time_double(['2013-05-01/07:20','2013-05-01/07:50'])
    probe = 'b'
    model_time = time_double('2013-05-01/07:38:03')
    beam_dis = 2.1
    
    
    vars = list()
    foreach species, ['p','o'] do begin
        var = rbsp_read_en_spec(time_range, probe=probe, species=species, pitch_angle_range=[0,45])
        vars.add, rename_var(var, output=var+'_para')
        var = rbsp_read_en_spec(time_range, probe=probe, species=species, pitch_angle_range=[45,135])
        vars.add, rename_var(var, output=var+'_perp')
        var = rbsp_read_en_spec(time_range, probe=probe, species=species, pitch_angle_range=[135,180])
        vars.add, rename_var(var, output=var+'_anti')
    endforeach
    vars = vars.toarray()
    options, vars, 'zrange', [1e4, 1e6]
    options, vars, 'color_table', 40


;---Plot H+ and O+ spec.
    sgopen, 0, xsize=8, ysize=10
    nvar = n_elements(vars)
    margins = [15,4,12,1]
    poss = sgcalcpos(nvar, margins=margins)
    tplot, vars, trange=time_range, position=poss
    times = make_bins(time_range, 5*60, inner=1)
    timebar, times, color=sgcolor('white'), linestyle=1
    constants = [1e1,1e2,1e3,1e4]
    yrange = [1,5e4]
    xrange = time_range
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach ty, constants do begin
            oplot, xrange, ty+[0,0], linestyle=1, color=sgcolor('white')
        endforeach
    endfor
    event_time = time_double('2013-05-01/07:38:03')
    timebar, event_time, color=sgcolor('red')



;---Select trace info.
    test_info_list = list()
    test_info_list.add, dictionary($
        'species', 'o', $
    ;    'trs', time_double('2013-05-01/'+['07:40:53','07:43:09']), $
    ;    'ens', [6200,1200] )
        'trs', time_double('2013-05-01/'+['07:40:53','07:43:09']), $
        'ens', [6200,1300] )
    test_info_list.add, dictionary($
        'species', 'p', $
    ;    'trs', time_double('2013-05-01/'+['07:39:22','07:40:30']), $
    ;    'ens', [6000,1000] )
        'trs', time_double('2013-05-01/'+['07:39:22','07:40:53']), $
        'ens', [6000,300] )
    
    foreach info, test_info_list do begin
        species = info['species']
        trs = info['trs']
        ens = info['ens']
        
        trace_input_list = list()
        
        prefix = 'rbsp'+probe+'_'
        north_var = prefix+species+'_en_spec_anti'
        south_var = prefix+species+'_en_spec_para'
        test_time_range = trs
        energy_range = [100,50000]
        pitch_angle_range = [135,180]
        
    ;---To get the trace info.
        files = rbsp_load_hope(test_time_range, id='l3%pa', probe=probe, errmsg=errmsg)
        var_list = list()
        suffix = (species eq 'e')? '_Ele': '_Ion'
        time_var = 'Epoch'+suffix
        energy_var = 'HOPE_ENERGY'+suffix
        flux_var = strupcase('f'+species+'du')
        var_list.add, dictionary($
            'in_vars', [energy_var,flux_var], $
            'time_var_name', time_var, $
            'time_var_type', 'Epoch' )
        read_vars, test_time_range, files=files, var_list=var_list, errmsg=errmsg
            
        pitch_angles = cdf_read_var('PITCH_ANGLE', filename=files[0])
        get_data, flux_var, test_times, fluxs
        energys = get_var_data(energy_var)
        index = where(abs(fluxs) ge 1e30, count)
        if count ne 0 then fluxs[index] = !values.f_nan
    
    
        log_ens = alog10(ens)
        pa_index = lazy_where(pitch_angles, '[]', pitch_angle_range)
        foreach time, test_times, time_id do begin
            the_fluxs = reform(fluxs[time_id,*,*])
            the_energys = energys[time_id,*]
            index = lazy_where(the_energys, '[]', energy_range)
            the_fluxs = (the_fluxs[index,*])[*,pa_index]
            the_energys = the_energys[index]
            the_pitch_angles = pitch_angles[pa_index]
            max = max(the_fluxs, index, nan=1)
            tmp = array_indices(the_fluxs, index)
            energy = the_energys[tmp[0]]
            energy = 10.^(log_ens[0]+(log_ens[1]-log_ens[0])/(trs[1]-trs[0])*(time-trs[0]))
            pitch_angle = the_pitch_angles[tmp[1]]
            
            trace_input_list.add, dictionary($
                'species', species, $
                'time', time, $
                'energy', energy, $
                'pitch_angles', !null, $
                'beam_dis', beam_dis, $
                'mod_time', model_time, $
                'model', 't01', $
                'igrf', 0 )
        endforeach

    
    
    
    ;---Plot selected trace info.
        index = where(vars eq north_var)
        tpos = poss[*,index]
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach trace_input, trace_input_list do begin
            time = trace_input['time']
            energy = trace_input['energy']
            plots, time, energy, data=1, psym=6, symsize=0.5, color=sgcolor('red')
        endforeach
    
    
    ;---Trace and plot output.
        pinfo = _2013_0501_load_data()
        probe = pinfo['probe']
        prefix = pinfo['prefix']
        r_gsm_var = prefix+'r_gsm'
        
        trace_output_list = list()
        foreach trace_input, trace_input_list do begin
            time = time_double(trace_input['time'])
            trace_input['r_gsm'] = get_var_data(r_gsm_var, at=time)
            trace_input['time'] = time
            par_var = trace_input['model']+'_var'
            trace_input['par'] = reform(get_var_data(par_var, at=time))
        
            trace_output = trace_ion_to_ionosphere(time, _extra=trace_input.tostruct())
            trace_output_list.add, trace_output
        endforeach
        
        foreach trace_output, trace_output_list, trace_id do begin
            trace_input = trace_input_list[trace_id]
            
            keys = ['north','south']
            foreach key, keys do begin
                var = (key eq 'north')? north_var: south_var
                index = where(vars eq var)
                tpos = poss[*,index]
                plot, xrange, yrange, $
                    xstyle=5, ystyle=5, ylog=1, $
                    nodata=1, noerase=1, position=tpos


                the_output = trace_output[key]

                conic_times = the_output['conic_time']
                beam_time = the_output['beam_time']
                nconic = n_elements(conic_times)
                colors = get_color(nconic+1)
                colors = sgcolor(['red','orange','yellow'])
                color_conics = colors[1:*]
                color_beam = colors[0]
                energy = trace_input['energy']


                foreach conic_time, conic_times, conic_id do begin
                    color_conic = color_conics[conic_id]
                    plots, conic_time, energy, data=1, psym=1, symsize=0.5, color=color_conic
                endforeach
                plots, beam_time, energy, data=1, psym=1, symsize=0.5, color=color_beam
            endforeach
            
        endforeach
    
    
    
        foreach trace_output, trace_output_list do begin
            foreach key, ['north','south'] do begin
                the_output = trace_output[key]
                print, the_output['conic_dis'], the_output['beam_dis']
            endforeach
        endforeach
    
    endforeach

end