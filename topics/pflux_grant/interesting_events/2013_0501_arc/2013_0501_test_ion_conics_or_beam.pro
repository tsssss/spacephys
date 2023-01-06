;+
; Trace H+ and O+ distribution backward in time to check if they are ion conics or beams.
;-

function trace_ion_to_ionosphere, time, species=species, model=model, mod_time=model_time, $
    energy=energy, pitch_angles=pitch_angles, r_gsm=r_gsm, igrf=igrf, $
    beam_dis=beam_dis, _extra=ex
    
    if n_elements(pitch_angles) eq 0 then pitch_angles = [175d,162]
    if n_elements(beam_dis) eq 0 then beam_dis = 2.1


    rad = constant('rad')
    deg = constant('deg')
    
    
    tmp = geopack_resolve_model(model)
    t89 = tmp.t89
    t96 = tmp.t96
    t01 = tmp.t01
    ts04 = tmp.ts04
    storm = tmp.storm
    
    if n_elements(model_time) eq 0 then model_time = time
    
    trace_result = dictionary()
    foreach trace_dir, [-1,1] do begin
        msg = (trace_dir eq -1)? 'north': 'south'
        
        if n_elements(par) eq 0 then begin
            time_range = model_time+[-1,1]*600
            if keyword_set(t89) then t89_par = 2
            par_var = geopack_read_par(time_range, model=model, t89_par=t89_par, _extra=ex)
            par = get_var_data(par_var, at=model_time)
        endif


        ps = geopack_recalc(time)
        xp = r_gsm[0]
        yp = r_gsm[1]
        zp = r_gsm[2]

        geopack_trace, xp,yp,zp, trace_dir, par, $
            xf,yf,zf, r0=r0, refine=1, ionosphere=1, $
            t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm, igrf=igrf, $
            fline=fline

        ntrace = n_elements(fline[*,0])
        ndim = 3
        bline = fltarr(ntrace,ndim)
        for ii=0,ntrace-1 do begin
            rx = fline[ii,0]
            ry = fline[ii,1]
            rz = fline[ii,2]

            if keyword_set(igrf) then begin
                geopack_igrf_gsm, rx,ry,rz, bx,by,bz
            endif else begin
                geopack_dip, rx,ry,rz, bx,by,bz
            endelse

            if model eq 'dip' or model eq 'dipole' or model eq 'igrf' then begin
                dbx = 0
                dby = 0
                dbz = 0
            endif else begin
                routine = 'geopack_'+model
                if model eq 't04s' then routine = 'geopack_ts04'
                call_procedure, routine, par, rx,ry,rz, dbx,dby,dbz
            endelse

            bline[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endfor


        bmag = snorm(bline)
        case species of
            'e': mass0 = 1d/1836
            'p': mass0 = 1d
            'o': mass0 = 16d
            'he': mass0 = 4d
        endcase
        mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.


    ;---Conics.
        re = constant('re')
        conic_times = []
        conic_dis = []
        fline_dis = snorm(fline)

;        if n_elements(pitch_angles) eq 0 then pitch_angles = [175d,162]
        foreach pitch_angle, pitch_angles do begin
            pa = pitch_angle*rad
            mu = sin(pa)^2/bmag[0]
            b_mirror = 1/mu
            r_mirror = sinterpol(fline, bmag, b_mirror)

            dis_mirror = (snorm(r_mirror))[0]
            index = where(fline_dis gt dis_mirror, nstep)
            v_para = fltarr(nstep)
            bmag_step = [bmag[0:nstep-1],b_mirror]
            bmag_step = (bmag_step[0:nstep-1]+bmag_step[1:nstep])*0.5
            r_step = [fline[0:nstep-1,*],r_mirror]
            dr_step = snorm(r_step[0:nstep-1,*]-r_step[1:nstep,*])

            vmag = sqrt(2*energy/mass0)*1e-3
            pitch_angle_step = asin(sqrt(mu*bmag_step))*deg
            vpara_step = vmag*cos(pitch_angle_step*rad)

            dt_step = dr_step*re/vpara_step
            dt_conic = total(dt_step)

            conic_times = [conic_times,dt_conic]
            conic_dis = [conic_dis,dis_mirror]
        endforeach


    ;---Beam.
        index = where(fline_dis gt beam_dis, nstep)
        v_para = fltarr(nstep)
        b_mirror = sinterpol(bmag, fline_dis, beam_dis)
        bmag_step = [bmag[0:nstep-1],b_mirror]
        bmag_step = (bmag_step[0:nstep-1]+bmag_step[1:nstep])*0.5
        r_mirror = sinterpol(fline, fline_dis, beam_dis)
        r_step = [fline[0:nstep-1,*],r_mirror]
        dr_step = snorm(r_step[0:nstep-1,*]-r_step[1:nstep,*])
        dt_beam = total(dr_step*re/vmag)
        beam_time = dt_beam

        trace_result[msg] = dictionary($
            'conic_dtime', conic_times, $
            'conic_dis', conic_dis, $
            'beam_dtime', beam_time, $
            'beam_dis', beam_dis )
    endforeach
    
    
    info = trace_result['north']
    info['conic_time'] = time-info['conic_dtime']
    info['beam_time'] = time-info['beam_dtime']
    
    info2 = trace_result['south']
    info2['conic_time'] = info['conic_time']+info2['conic_dtime']
    info2['beam_time'] = info['beam_time']+info2['beam_dtime']
    
    return, trace_result
    
end


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