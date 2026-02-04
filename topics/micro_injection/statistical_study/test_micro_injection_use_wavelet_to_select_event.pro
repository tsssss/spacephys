;+
; Check to use wavelet to select micro-injection events.
;-

    test = 0
    mission_probe = ['mms','4']
    time_range = ['2015-09-01','2015-09-02']
    plot_tr = time_double(['2015-09-01/10:00','2015-09-02'])
    event_id = time_string(time_double(time_range[0]), tformat='YYYY_MMDD')
    plot_dir = join_path([srootdir(),'plot',event_id])

    errmsg = ''
    retval = !null
    version = 'v01'

    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = mission+probe+'_'
    routine = mission+'_read_kev_electron'
    energy = 1.04e5 ; eV.
    scale_info = {s0:40d, s1:4000, dj:1d/8, ns:0d }

;---Load data.
    common_times = make_bins(time_double(time_range), 1d)

    flux_var = call_function(routine, time_range, probe=probe, energy=energy, errmsg=errmsg)
    ; Work in log scale.
    flux_var_2 = flux_var+'_log'
    fluxs = get_var_data(flux_var, times=times, energy_bins, limits=lim)
    store_data, flux_var_2, times, alog10(fluxs), energy_bins, limits=lim
    options, flux_var_2, ylog=0, ytitle='Log!D10!N('+lim.unit+')'    
    interp_time, flux_var_2, common_times

    ; B field.
    routine = mission+'_read_bfield'
    b_var = call_function(routine, time_range, probe=probe, errmsg=errmsg)
    interp_time, b_var, common_times
    ; E field.
    routine = mission+'_read_efield'
    e_var = call_function(routine, time_range, probe=probe, errmsg=errmsg)
    e_gse = get_var_data(e_var, times=times)
    e_gsm = cotran_pro(e_gse, times, coord_msg=['gse','gsm'])
    e_var = prefix+'e_gsm'
    store_data, e_var, times, e_gsm
    interp_time, e_var, common_times
    add_setting, e_var, smart=1, dictionary($
        'display_type', 'vector', $
        'coord', 'gsm', $
        'unit', 'mV/m', $
        'short_name', 'E' )
    
    
    ; Ion vel.
    routine = mission+'_read_ion_vel'
    u_var = call_function(routine, time_range, probe=probe, errmsg=errmsg)
    interp_time, u_var, common_times
    
    ; Orbit.
    routine = mission+'_read_orbit'
    orbit_var = call_function(routine, time_range, probe=probe, errmsg=errmsg)
    interp_time, orbit_var, common_times


    ; Convert to FAC.
    external_model = 't89'
    internal_model = 'dipole'
    b0_window = 1200d
    bmod_var = lets_read_geopack_bfield(orbit_var=orbit_var, external_model=external_model, internal_model=internal_model)
    b_vars = lets_decompose_bfield(b_var=b_var, b0_window=b0_window, bmod_var=bmod_var)
    b0_var = b_vars['b0']
    b1_var = b_vars['b1']
    ; Delete B field data close to Earth.
    min_dis = 4.
    foreach var, [b1_var] do begin
        data = get_var_data(var, times=times)
        dis = snorm(get_var_data(orbit_var, at=times))
        index = where(dis le min_dis, count)
        if count ne 0 then begin
            data[index,*] = !values.f_nan
            store_data, var, times, data
        endif
    endforeach

    ; Convert to FAC.
    default_coord = 'gsm'
    fac_coord = mission+'_fac'
    q_fac_var = lets_define_fac(r_var=orbit_var, b_var=b0_var, fac_coord=fac_coord)
    coord_msgs = [default_coord,fac_coord]
    fac_vars = list()
    foreach var, [b1_var,e_var,u_var] do begin
        options, var, mission='mms'
        out_var = streplace(var,default_coord,fac_coord)
        fac_vars.add, lets_cotran(coord_msgs, input=var, q_var=q_fac_var, output=out_var)
    endforeach
    fac_vars = fac_vars.toarray()


;---Calc wavelet.
    fac_labels = ['b','w','o']

    spec_var = stplot_mor_new(flux_var_2, scale_info=scale_info)
    spec = get_var_data(spec_var, times=times, freqs, limits=lim)
    spec_var_2 = spec_var+'_plot'
    store_data, spec_var_2, times, spec, freqs*1e3, limits=lim
    options, spec_var_2, ytitle='Freq (mHz)', yrange=lim.yrange*1e3, ztitle='PSD Ele flux'

    vars = prefix+['b1','e','u']+'_'+fac_coord
    wanted_comps = list(['b','w','o'],['w','o'],['w','o'])
    ct = !null
    zrange = !null
    foreach var, vars, vid do begin
        vec = get_var_data(var, times=times, settings=settings)
        wanted_comp = wanted_comps[vid]
        foreach comp, wanted_comp, cid do begin
            wanted_index = where_pro(fac_labels, 'eq', comp)
            field_var = var+'_'+comp
            mor_var = field_var+'_mor'
            if tnames(mor_var) ne '' then continue
            
            dat = vec[*,wanted_index]
            store_data, field_var, times, dat
            spec_var = stplot_mor_new(field_var, scale_info=scale_info)
            get_data, spec_var, times, data, freqs
            store_data, spec_var, times, data, freqs*1e3
            unit = settings.unit
            short_name = settings.short_name
            add_setting, spec_var, smart=1, dictionary($
                'requested_time_range', time_range, $
                'no_interp', 1, $
                'display_type', 'spec', $
                'unit', unit, $
                'ytitle', 'Freq (mHz)', $
                'yrange', minmax(freqs*1e3), $
                'ylog', 1, $
                'zlog', 1, $
                'short_name', short_name )
            ztitle = short_name+'('+unit+')!U2!N'
            options, spec_var, ztitle=ztitle, color_table=ct, zrange=zrange
        endforeach
    endforeach

    

;---Calc PSD.
    suffix = ['','_?_mor']
    plot_vars = [flux_var,spec_var_2,$
        prefix+'b1_mms_fac'+suffix,prefix+'e_mms_fac'+suffix,prefix+'u_mms_fac'+suffix]
    plot_file = join_path([plot_dir,'test_micro_injection_use_wavelet_to_select_event_'+event_id+'_overview_v01.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, size=[6,8]
    tplot, plot_vars, trange=plot_tr
    if keyword_set(test) then stop
    sgclose
    
    
    window = 1200d   ; sec.
    psd_times = make_bins(time_double(time_range), window)
    nsector = n_elements(psd_times)-1

    var_list = list()
    var_list.add, prefix+'kev_e_flux_log_mor_plot'
    var_list.add, prefix+'b1_mms_fac_'+['b','w','o']+'_mor'
    var_list.add, prefix+'e_mms_fac_'+['w','o']+'_mor'
    var_list.add, prefix+'u_mms_fac_'+['w','o']+'_mor'
    var_strs = ['kev_ele','b1','e','ion_vel']
    foreach vars, var_list, vid do begin
        nvar = n_elements(vars)
        specs = get_var_data(vars[0], times=times, freqs)
        ntime = n_elements(times)
        nfreq = n_elements(freqs)
        all_specs = fltarr(ntime,nfreq,nvar)
        foreach var, vars, ii do all_specs[*,*,ii] = get_var_data(var)
        var_str = var_strs[vid]

        psds = fltarr(nsector,nfreq,nvar)
        flags = fltarr(nsector)
        for ii=0,nsector-1 do begin
            tr = psd_times[ii:ii+1]
            index = where_pro(times,'[]',tr, count=count)
            if count eq 0 then continue
            
            psds[ii,*,*] = mean(all_specs[index,*,*],dimension=1)
            
            plot_file = join_path([plot_dir,var_str+'_psd','test_micro_injection_use_wavelet_to_select_event_'+event_id+'_'+var_str+'_psd_'+time_string(tr[0],tformat='YYYY_MMDD_hhmm')+'_v01.pdf'])
            if keyword_set(test) then plot_file = 0
            sgopen, plot_file, size=[6,6]
            margins = [10,4,12,1]
            poss = sgcalcpos(2, margins=margins, ypad=4)
            tpos = poss[*,0]
            var = vars[0]
            line_var = strmid(var,0,strpos(var,'fac')+3)
            if nvar eq 1 then line_var = strmid(var,0,strpos(var,'_mor'))
            tplot, line_var, trange=tr, position=tpos
            tpos = poss[*,1]

            xtitle = 'Freq (mHz)'
            ytitle = 'PSD'
            xxs = freqs
            xrange = minmax(xxs)
            yrange = [0.1,1e6]

            if nvar eq 1 then begin
                yys = reform(psds[ii,*])
                ;yys *= 1e-3
            endif else begin
                yys = mean(reform(psds[ii,*,*]),dimension=2)
            endelse
            plot, xxs, yys, $
                xlog=1, ylog=1, position=tpos, noerase=1, $
                xtitle=xtitle, xstyle=1, xrange=xrange, $
                ytitle=ytitle, ystyle=1, yrange=yrange
            
            
            ; Check original data to see if the peak is real.
            is_real = 1
            orig_data = get_var_data(line_var, in=tr, times=uts)
            index = where(finite(orig_data,nan=1), count)
            ndata = n_elements(uts)
            if count gt 0.1*ndata then is_real = 0
            
            ; Check PSD max value within the interested range.
            freq_range = [1,10.]   ; mHz
            index = where_pro(freqs,'[]', freq_range, count=count)
            if count eq 0 then is_real = 0
            the_psd = yys
            psd_max = max(the_psd[index], max_index)
            index_lim = 3
            if max_index lt index_lim then is_real = 0 ; Too close to the low freq.
            if max_index gt n_elements(index)-index_lim then is_real = 0 ; Too close to the high freq.
            
            the_index = where(the_psd eq psd_max)
            if is_real then begin
                plots, xxs[the_index], yys[the_index], psym=1, data=1, color=sgcolor('red')
                flags[ii] = 1
            endif else begin
                plots, xxs[the_index], yys[the_index], psym=1, data=1, color=sgcolor('green')
                flags[ii] = 0
            endelse
            if keyword_set(test) then stop
            sgclose
        endfor

        
        plot_file = join_path([plot_dir,var_str+'_psd','test_micro_injection_use_wavelet_to_select_event_'+event_id+'_'+var_str+'_result_v01.pdf'])
        if keyword_set(test) then plot_file = 0
        sgopen, plot_file, size=[6,4]
        thick = 2
        plot_vars = [line_var,vars]
        nplot_var = n_elements(plot_vars)
        poss = sgcalcpos(nplot_var, margins=margins)
        
        tplot, plot_vars, trange=plot_tr, position=poss
        tpos = poss[*,0]
        set_axis, line_var, position=tpos, xrange=plot_tr, yrange=[0,1]
        for ii=0,nsector-1 do begin
            if flags[ii] eq 0 then continue
            tr = psd_times[ii:ii+1]
            oplot, tr, [0,0], color=sgcolor('red'), thick=thick;, data=1
        endfor
        if keyword_set(test) then stop
        sgclose

    endforeach


end