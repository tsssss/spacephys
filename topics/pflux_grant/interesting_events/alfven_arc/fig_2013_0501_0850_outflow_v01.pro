;+
; To plot H and O outflow and pitch angle distributions.
;-

function fig_2013_0501_0850_outflow_v01, event_info=event_info

test = 1

;---Load data.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_0850_load_data()


;---Settings.
    prefix = event_info['prefix']
    probe = event_info['probe']
    time_range = event_info['time_range']+[0,1800]
    
    bar_thick = (keyword_set(test))? 2: 6
    label_size = 0.8
    symsize = 0.5

    ct_oxygen = 64
    ct_proton = 63
    zrange_proton = [1e3,1e6]
    zrange_oxygen = [1e3,1e6]

    o_spec_var = prefix+'o_en_spec_'+['anti','perp','para']
    vars = o_spec_var
    options, vars, 'color_table', ct_oxygen
    options, vars, 'zrange', zrange_oxygen
    options, vars, 'yrange', [10,5e4]
    foreach var, vars do begin
        get_data, var, times, data, val
        index = where(data eq 0, count)
        if count ne 0 then begin
            data[index] = 0.01
            store_data, var, times, data, val
        endif
    endforeach

    
    p_spec_var = prefix+'p_en_spec_'+['anti','perp','para']
    vars = p_spec_var
    options, vars, 'color_table', ct_proton
    options, vars, 'zrange', zrange_proton
    options, vars, 'yrange', [10,5e4]
    foreach var, vars do begin
        get_data, var, times, data, val
        index = where(data eq 0, count)
        if count ne 0 then begin
            data[index] = 0.01
            store_data, var, times, data, val
        endif
    endforeach

    

    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,'fig_2013_0501_0850_outflow_v01.pdf'])
    if keyword_set(test) then plot_file = 0
    
    
    components = ['anti']
    plot_vars = prefix+['o_en_spec_'+components,'p_en_spec_'+components]
    foreach var, plot_vars do copy_data, var, var+'_clean'
    plot_vars = [plot_vars+'_clean',plot_vars]
    nplot_var = n_elements(plot_vars)
    margins = [12,4,12,1]
    poss = sgcalcpos(nplot_var, margins=margins)
    sgopen, plot_file, size=[6,6]
    tplot, plot_vars, trange=time_range, position=poss

    bar_times = time_double('2013-05-01/'+['07:38','08:44'])
    timebar, bar_times, color=sgcolor('red'), linestyle=1



;---Identify the energy and time of the outflow.
    obs_info = dictionary()
    obs_info['o1'] = dictionary($
        'species', 'o', $
        'spec_var', prefix+'o_en_spec_anti', $
        'nanti_appearance', 2, $
        'npara_appearance', 1, $
        'trs', time_double('2013-05-01/'+['07:40:53','07:43:32']), $
        'ens', [6200,900] )
    obs_info['o2'] = dictionary($
        'species', 'o', $
        'nanti_appearance', 1, $
        'npara_appearance', 0, $
        'spec_var', prefix+'o_en_spec_anti', $
        'trs', time_double('2013-05-01/'+['08:45:57','08:51:37']), $
        'ens', [4000,300] )
;    obs_info['o3'] = dictionary($
;        'species', 'o', $
;        'spec_var', prefix+'o_en_spec_anti', $
;        'trs', time_double('2013-05-01/'+['08::49:44','08:52:00']), $
;        'ens', [900,120] )



    foreach the_info, obs_info do begin

        get_data, the_info['spec_var'], times, data, vals, limits=lim
        time_index = lazy_where(times, '[]', the_info['trs'], count=ntest_time)
        if ntest_time eq 0 then message, 'Inconsistency ...'

        test_times = times[time_index]
        test_energys = 10.^interpol(alog10(the_info['ens']), the_info['trs'], test_times)
        the_info['test_times'] = test_times
        the_info['test_energys'] = test_energys
        
;        energy_range = minmax(the_info['ens'])
;        test_energys = fltarr(ntest_time)
;        foreach ii, time_index, test_id do begin
;            energy_index = lazy_where(reform(vals[ii,*]), '[]', energy_range, count=count)
;            if count eq 0 then message, 'Inconsistency ...'
;            max_flux = max(data[ii,energy_index], index)
;            test_energys[test_id] = vals[ii,energy_index[index]]
;        endforeach

        pid = where(plot_vars eq the_info['spec_var'])
        tpos = poss[*,pid]
        yrange = lim.yrange
        xrange = time_range
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach test_time, test_times, ii do begin
            plots, test_times[ii], test_energys[ii], psym=1, symsize=symsize, color=sgcolor('red')
        endforeach
    endforeach
    
    
;---Solve for beam generation altitude and time.
    r_gsm_var = prefix+'r_gsm'
    foreach the_info, obs_info do begin
        model_time = (the_info['trs'])[0]
        external_model = 't89'
        internal_model = 'dipole'
        r_gsm = get_var_data(r_gsm_var, at=model_time)
        par = 2d

    ;---Get the model field line.
        time = model_time
        model = external_model
        fline_info = dictionary()
        foreach trace_dir, [-1,1] do begin
            msg = (trace_dir eq -1)? 'north': 'south'

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
                    dbx = -1
                    dby = -1
                    dbz = -1
                endif else begin
                    routine = 'geopack_'+model
                    if model eq 't03s' then routine = 'geopack_ts04'
                    call_procedure, routine, par, rx,ry,rz, dbx,dby,dbz
                endelse

                bline[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
            endfor

            fline_info[msg] = dictionary($
                'fline', fline, $
                'bline', bline )
        endforeach
        
    ;---Solve each energy bin.
        species = the_info.species
        case species of
            'e': mass0 = 1d/1836
            'p': mass0 = 1d
            'o': mass0 = 16d
            'he': mass0 = 4d
        endcase
        mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.

        test_energys = the_info.test_energys
        test_times = the_info.test_times
        ntest_energy = n_elements(test_energys)

    ;---Solve northern hemisphere.
        flines = fline_info.north.fline
        blines = fline_info.north.bline
        bmags = snorm(blines)
        nfline = n_elements(blines[*,0])
        trace_times = dblarr(nfline,ntest_energy)

        foreach energy, test_energys, test_id do begin
            vmag = sqrt(2*energy/mass0)*1e-3
            drs = snorm(flines[1:nfline-1,*]-flines[0:nfline-2,*])*constant('re')
            dts = drs/vmag
            trace_times[0,test_id] = test_times[test_id]
            for ii=1,nfline-1 do trace_times[ii,test_id] = trace_times[ii-1,test_id]-dts[ii-1]
        endforeach

        fline_diss = snorm(flines)
        dis_step = 0.1
        uniform_diss = make_bins(minmax(fline_diss), dis_step)
        nuniform_dis = n_elements(uniform_diss)
        uniform_times = sinterpol(trace_times, fline_diss, uniform_diss)
        time_diffs = dblarr(nuniform_dis)
        for ii=0,nuniform_dis-1 do begin
            time_diffs[ii] = stddev(uniform_times[ii,*])
        endfor
        beam_time_error = min(time_diffs, index)
        beam_dis = uniform_diss[index]
        beam_time = mean(uniform_times[index,*])
        the_info['beam_time'] = beam_time
        the_info['beam_dis'] = beam_dis
        the_info['beam_time_error'] = beam_time_error
        
        plots, uniform_times[index,*], test_energys, psym=1, symsize=symsize, color=sgcolor('blue')
        print, beam_dis
        print, time_string(beam_time)
        print, beam_time_error


    ;---Basic travel time.
        travel_info = dictionary()

        ; North.
        test_times = the_info.test_times
        beam_time = the_info.beam_time
        travel_times = abs(test_times-beam_time)
        travel_info['north'] = travel_times

        ; South.
        flines = fline_info.south.fline
        blines = fline_info.south.bline
        bmags = snorm(blines)
        nfline = n_elements(blines[*,0])
        trace_times = dblarr(nfline,ntest_energy)

        foreach energy, test_energys, test_id do begin
            vmag = sqrt(2*energy/mass0)*1e-3
            drs = snorm(flines[1:nfline-1,*]-flines[0:nfline-2,*])*constant('re')
            dts = drs/vmag
            trace_times[0,test_id] = test_times[test_id]
            for ii=1,nfline-1 do trace_times[ii,test_id] = trace_times[ii-1,test_id]-dts[ii-1]
        endforeach

        fline_diss = snorm(flines)
        beam_times = dblarr(ntest_energy)
        for ii=0,ntest_energy-1 do beam_times[ii] = sinterpol(trace_times[*,ii], fline_diss, beam_dis)
        travel_times = abs(test_times-beam_times)
        travel_info['south'] = travel_times
        the_info['travel_info'] = travel_info
    endforeach



;---Trace multiple times.
    south_color = sgcolor('magenta')
    north_color = sgcolor('red')
    coef = 1.0
    L_phi = 0.5
    dphi = 5
    dt_phi = 117d*L_phi/sqrt(dphi)    ; time to go through potential drop one time, t = 117*L/sqrt(phi), L in Re, phi in kV.
    

    pid = where(plot_vars eq prefix+'o_en_spec_anti', count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        yrange = lim.yrange
        xrange = time_range
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach the_info, obs_info do begin
            beam_time = the_info.beam_time
            test_energys = the_info.test_energys
            travel_info = the_info.travel_info
            nanti_appearance = the_info['nanti_appearance']
            
            ; North.
            hem_color = north_color
            t0 = beam_time
            ; 1st appearance.
            if nanti_appearance ge 1 then begin            
                t1 = t0+travel_info.north
                plots, t1, test_energys, psym=1, symsize=symsize, color=hem_color
            endif
            ; 2nd appearnce.
            if nanti_appearance ge 2 then begin
                t2 = t1+travel_info.south+dt_phi+$ ; before reaching southern hemisphere.
                    (travel_info.south+travel_info.north*2+dt_phi*3)/sqrt(coef)  ; after reaching southern hemisphere.
                plots, t2, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
            
            ; South.
            hem_color = south_color
            t0 = beam_time
            ; 1st appearance.
            t1 = t0+(travel_info.south+travel_info.north*2+dt_phi*2)/sqrt(coef)
            plots, t1, test_energys*coef, psym=1, symsize=symsize, color=hem_color
        endforeach
    endif

    
;---Predict H+.
    pid = where(plot_vars eq prefix+'p_en_spec_anti', count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        yrange = lim.yrange
        xrange = time_range
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach the_info, obs_info do begin
            beam_time = the_info.beam_time
            test_energys = the_info.test_energys
            travel_info = the_info.travel_info
            nanti_appearance = the_info['nanti_appearance']
    
            ; North.
            hem_color = north_color
            t0 = beam_time
            ; 1st appearance.
            if nanti_appearance ge 1 then begin
                t1 = t0+(travel_info.north)/4
                plots, t1, test_energys, psym=1, symsize=symsize, color=hem_color
            endif
            ; 2nd appearnce.
            if nanti_appearance ge 2 then begin
                t2 = t1+(travel_info.south+$ ; before reaching southern hemisphere.
                    (travel_info.south+travel_info.north*2+dt_phi*3)/sqrt(coef))/4  ; after reaching southern hemisphere.
                plots, t2, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
            
            ; South.
            hem_color = south_color
            t0 = beam_time
            ; 1st appearance.
            t1 = t0+(travel_info.south+travel_info.north*2+dt_phi*2)/sqrt(coef)
            plots, t1, test_energys*coef, psym=1, symsize=symsize, color=hem_color
        endforeach
    endif
    
    
    ; parallel.
    pid = where(plot_vars eq prefix+'o_en_spec_para', count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        yrange = lim.yrange
        xrange = time_range
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach the_info, obs_info do begin
            beam_time = the_info.beam_time
            test_energys = the_info.test_energys
            travel_info = the_info.travel_info
            npara_appearance = the_info['npara_appearance']

            ; South.
            hem_color = south_color
            t0 = beam_time
            ; 1st appearance.
            if nanti_appearance ge 1 then begin            
                t1 = t0+(travel_info.south)/sqrt(coef)
                plots, t1, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
            ; 2nd appearance.
            if nanti_appearance ge 2 then begin
                t2 = t1+(travel_info.north+travel_info.south+dt_phi*2)*2/sqrt(coef)
                plots, t2, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
            
            ; North.
            hem_color = north_color
            t0 = beam_time
            ; 1st appearance.
            if nanti_appearance ge 1 then begin
                t1 = t0+travel_info.north+travel_info.south+dt_phi+(travel_info.south+dt_phi)/sqrt(coef)
                plots, t1, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
            ; 2nd appearance.
            if nanti_appearance ge 2 then begin
                t2 = t1+(travel_info.north+travel_info.south+dt_phi*2)*2/sqrt(coef)
                plots, t2, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
        endforeach
    endif
    
    ; parallel.
    npara_appearance = 1
    pid = where(plot_vars eq prefix+'p_en_spec_para', count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        yrange = lim.yrange
        xrange = time_range
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach the_info, obs_info do begin
            beam_time = the_info.beam_time
            test_energys = the_info.test_energys
            travel_info = the_info.travel_info
            npara_appearance = the_info['npara_appearance']

            ; South.
            hem_color = south_color
            t0 = beam_time
            ; 1st appearance.
            if npara_appearance ge 1 then begin
                t1 = t0+(travel_info.south)/sqrt(coef)/4
                plots, t1, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
            ; 2nd appearance.
            if npara_appearance ge 2 then begin
                t2 = t1+(travel_info.north+travel_info.south)*2/sqrt(coef)/4
                plots, t2, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
            

            ; North.
            hem_color = north_color
            t0 = beam_time
            ; 1st appearance.
            if npara_appearance ge 1 then begin
                t1 = t0+(travel_info.north+travel_info.south+(travel_info.south)/sqrt(coef))/4
                plots, t1, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
            ; 2nd appearance.
            if npara_appearance ge 2 then begin
                t2 = t1+((travel_info.north+travel_info.south)*2/sqrt(coef))/4
                plots, t2, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
        endforeach
    endif


    if keyword_set(test) then stop
    sgclose

    return, plot_file

end


print, fig_2013_0501_0850_outflow_v01(event_info=event_info)
end
