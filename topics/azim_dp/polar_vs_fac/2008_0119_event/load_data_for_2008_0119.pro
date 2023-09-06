;+
; Load data for the 2008_0119 event.
;-

pro load_data_for_2008_0119_trace, default_settings

    trace_var = 'tracing_infos'
    if ~check_if_update(trace_var) then return

    tracing_infos = hash()

;---LANL 97A.
    probe = 'LANL-97A'
    injections = orderedhash()
    aurora_time = '2008-01-19/06:28'
    injections[aurora_time] = list($
        list(140, '2008-01-19/06:44:15'), $
        list(206, '2008-01-19/06:43:15'), $
        list(316, '2008-01-19/06:42:15') )
    aurora_time = '2008-01-19/07:01'
    injections[aurora_time] = list($
        list(140, '2008-01-19/07:09:15'), $
        list(206, '2008-01-19/07:08:15') )
    aurora_time = '2008-01-19/07:14'
    injections[aurora_time] = list($
        list(140, '2008-01-19/07:26:05'), $
        list(206, '2008-01-19/07:24:50'), $
        list(316, '2008-01-19/07:23:55') )
    tracing_infos[probe] = dictionary($
        'probe', probe, $
        'species', 'p', $
        'injections', injections)


;---LANL 1989-046.
    probe = '1989-046'
    injections = orderedhash()
    aurora_time = '2008-01-19/06:28'
    injections[aurora_time] = list($
        list(!null, '2008-01-19/06:34:00'), $
        list(!null, '2008-01-19/06:39:42'))
    aurora_time = '2008-01-19/07:01'
    injections[aurora_time] = list($
        list(176, '2008-01-19/07:06:00'), $
        list( 60, '2008-01-19/07:08:20'))
    aurora_time = '2008-01-19/07:14'
    injections[aurora_time] = list($
        list(!null, '2008-01-19/07:21:30'))
    tracing_infos[probe] = dictionary($
        'probe', probe, $
        'species', 'e', $
        'injections', injections)

;---THA.
    probe = 'tha'
    injections = orderedhash()
    aurora_time = '2008-01-19/06:28'
    injections[aurora_time] = list($
        list(140, '2008-01-19/06:37:39'), $
        list( 66, '2008-01-19/06:41:09'))
    aurora_time = '2008-01-19/07:01'
    injections[aurora_time] = list($
        list(204, '2008-01-19/07:04:40'), $
        list( 66, '2008-01-19/07:09:00'))
    aurora_time = '2008-01-19/07:14'
    injections[aurora_time] = list($
        list(204, '2008-01-19/07:22:46'), $
        list( 66, '2008-01-19/07:26:46'))
    aurora_time = '2008-01-19/07:39'
    injections[aurora_time] = list($
        list(293, '2008-01-19/07:47:00'), $
        list( 66, '2008-01-19/07:53:00'))
    tracing_infos[probe] = dictionary($
        'probe', probe, $
        'species', 'e', $
        'injections', injections)
        
;---THE.
    probe = 'the'
    injections = orderedhash()
    aurora_time = '2008-01-19/06:28'
    injections[aurora_time] = list($
        list(293, '2008-01-19/06:35:30'), $
        list( 93, '2008-01-19/06:37:20'))
    aurora_time = '2008-01-19/07:01'
    injections[aurora_time] = list($
        list(204, '2008-01-19/07:04:00'), $
        list( 83, '2008-01-19/07:04:09'))
    aurora_time = '2008-01-19/07:14'
    injections[aurora_time] = list($
        list( 41, '2008-01-19/07:20:21'), $
        list(204, '2008-01-19/07:18:51'))
    aurora_time = '2008-01-19/07:39'
    injections[aurora_time] = list($
        list(204, '2008-01-19/07:44:21'), $
        list( 52, '2008-01-19/07:46:00'))
    tracing_infos[probe] = dictionary($
        'probe', probe, $
        'species', 'e', $
        'injections', injections)

;---Trace injections backward in time.
    time_step = 4d
    ndim = 3
    foreach tracing_info, tracing_infos do begin
        probe = tracing_info['probe']
        species = tracing_info['species']
        is_ion = (species eq 'e')? 0: 1
        prefix = probe+'_'
        
        flux_var = prefix+species+'_flux'
        get_data, flux_var, times, fluxs, energys
        energy_ratio = sqrt(mean(energys[1:-1]/energys[0:-2]))
        
        r_gsm_var = prefix+'r_gsm'
        r_sm_var = prefix+'r_sm'
        if check_if_update(r_gsm_var) then begin
            get_data, r_sm_var, times, r_sm
            r_gsm = cotran(r_sm, times, 'sm2gsm')
            store_data, r_gsm_var, times, r_gsm
        endif
        
        injection_infos = tracing_info['injections']
        start_times = injection_infos.keys()
        trace_result = orderedhash()
        foreach start_time_str, start_times do begin
            injections = injection_infos[start_time_str]

            ; dispersionless injection.
            if n_elements((injections[0])[0]) eq 0 then continue

            the_trace = dictionary()
            the_trace['energy_ratio'] = energy_ratio

            ninjection = n_elements(injections)
            injection_time_strs = strarr(ninjection)
            trace_info = list()
            foreach injection, injections, injection_id do begin
                injection_time_strs[injection_id] = injection[1]
                trace_info.add, list()
            endforeach
            
            ; Trace each injection to start time and save the traced positions.
            start_time = time_double(start_time_str)
            foreach injection, injections, injection_id do begin
                energy = injection[0]
                if n_elements(energy) eq 0 then break   ; dispersionless.
                
                energy *= energy_ratio  ; use the maximum energy of the energy bin.
                end_time = time_double(injection[1])
                trace_times = reverse(smkarthm(start_time, end_time, time_step, 'dx'))
                r_gsm = get_var_data(r_gsm_var, at=end_time)
                trace_r_gsms = trace_injection_backward(r_gsm, trace_times, energy, ion=is_ion)
                trace_mlts = pseudo_mlt(cotran(trace_r_gsms, trace_times, 'gsm2sm'))
                ; This updates injections and even tracing_infos and tracing_info.
                injection.add, trace_times      ; injection[2].
                injection.add, trace_r_gsms     ; injection[3].
                injection.add, trace_mlts       ; injection[4].
            endforeach
            
            ; Solve for the injection time and mlt.
            injection_times = time_double(injection_time_strs)
            common_times = smkarthm(start_time,min(injection_times), time_step, 'dx')
            ncommon_time = n_elements(common_times)
            mlts = fltarr(ncommon_time,ninjection)
            foreach injection, injections, injection_id do begin
                energy = injection[0]
                if n_elements(energy) eq 0 then continue    ; dispersionless.
                trace_times = injection[2]
                trace_mlts = injection[4]
                mlts[*,injection_id] = interpol(trace_mlts, trace_times, common_times)
            endforeach
            mlt_devs = fltarr(ncommon_time)
            for ii=0, ncommon_time-1 do begin
                mlt_devs[ii] = stddev(mlts[ii,*])
            endfor
            min_mlt_dev = min(mlt_devs, index)
            injection_time = common_times[index]
            injection_mlt = mean(mlts[index,*])
            print, time_string(injection_time)
            print, injection_mlt
            print, energy_ratio
            the_trace['injection_time'] = injection_time
            the_trace['injection_mlt'] = injection_mlt

            trace_result[start_time_str] = the_trace
        endforeach
        tracing_info['trace_result'] = trace_result
    endforeach
    


    store_data, trace_var, 0, tracing_infos



end

pro load_data_for_2008_0119_injection, default_settings

    injection_var = 'injection_result'
    if ~check_if_update(injection_var) then return

    tracing_infos = get_var_data('tracing_infos')
    injection_result = hash()
    foreach probe, tracing_infos.keys() do begin
        r_gsm_var = probe+'_r_gsm'
        tracing_info = tracing_infos[probe]
        injection_info = orderedhash()

        trace_result = tracing_info['trace_result']
        raw_injections = tracing_info['injections']
        aurora_time_strs = raw_injections.keys()
        foreach aurora_time_str, aurora_time_strs do begin
            injections = raw_injections[aurora_time_str]

            foreach injection, injections do begin
                if n_elements(injection[0]) eq 0 then begin
                    injection_time_str = injection[1]
                    injection_time = time_double(injection_time_str)
                    r_gsm = get_var_data(r_gsm_var, at=injection_time)
                    injection_mlt = pseudo_mlt(cotran(r_gsm, injection_time, 'gsm2sm'))
                    injection_info[injection_time_str] = dictionary($
                        'type', 'dispersionless', $
                        'time', injection_time, $
                        'mlt', injection_mlt, $
                        'injections', injection )
                endif else begin
                    trace = trace_result[aurora_time_str]
                    injection_time = trace['injection_time']
                    injection_time_str = time_string(injection_time)
                    injection_mlt = trace['injection_mlt']
                    if injection_info.haskey(injection_time_str) then continue
                    injection_info[injection_time_str] = dictionary($
                        'type', 'dispersive', $
                        'time', injection_time, $
                        'mlt', injection_mlt, $
                        'injections', injections )
                endelse
            endforeach
            injection_result[probe] = injection_info
        endforeach
    endforeach

    store_data, injection_var, 0, injection_result

end



pro load_data_for_2008_0119_footpoint, default_settings


;---Settings.
    time_range = default_settings.event_time_range+[-1,1]*1800d
    probes = [['1989-046','1994-084','LANL-97A'],'th'+['a','d','e'],'g'+['11','12']]
    models = ['t89','t96','t01','t04s']
    default_settings['mapping'] = dictionary($
        'probes', probes, $
        'models', models )
    data_file = default_settings.data_file
    
    foreach probe, probes do begin
        prefix = probe+'_'
        r_sm_var = prefix+'r_sm'
        get_data, r_sm_var, times, r_sm, limits=lim
        r_gsm = cotran(r_sm, times, 'sm2gsm')
        r_gsm_var = prefix+'r_gsm'        
        store_data, r_gsm_var, times, r_gsm
        add_setting, r_gsm_var, smart=1, dictionary($
            'short_name', 'R', $
            'unit', 'Re', $
            'display_type', 'vector', $
            'coord', 'GSM' )

        time_var = 'time_mapping_'+probe
        types = ['fmlt','fmlat','fmlon']
        fvars = prefix+types
        load = 0
        foreach var, prefix+types do begin
            if ~cdf_has_var(var, filename=data_file) then begin
                load = 1
                break
            endif
        endforeach
        
        if load then begin
            colors = sgcolor(['black','red','green','blue'])

            foreach model, models do begin
                read_geopack_info, r_gsm_var, model=model, prefix=prefix, suffix='_'+model, direction=-1
            endforeach
            foreach var, prefix+types do begin
                stplot_merge, var+'_'+models, newname=var, labels=strupcase(models), colors=colors
                options, var, 'ynozero', 1
                options, var, 'labflag', -1
            endforeach

            var = prefix+types[0]
            get_data, var, times
            time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
            cdf_save_var, time_var, value=times, filename=data_file
            cdf_save_setting, time_settings, varname=time_var, filename=data_file

            foreach var, prefix+types do begin
                get_data, var, times, data
                settings = dictionary(lim)
                settings['depend_0'] = time_var
                cdf_save_var, var, value=data, filename=data_file
                cdf_save_setting, settings, varname=var, filename=data_file
            endforeach
        endif
        foreach var, prefix+types do begin
            cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
        endforeach

        
    endforeach
    

end

pro load_data_for_2008_0119_themis_uniform_mlt_image_uniform, default_settings

    return
;---Settings.
    time_range = default_settings.event_time_range
    sites = ['fsim','gill','inuv']
    min_elev = 2.5  ; deg
    mlat_range = [60d,80]   ; deg.
    mlt_range = [-6d,3]     ; hr.
    merge_method = 'merge_elev'
    var = 'thg_asf_mlt_image_uniform'
    time_var = 'time_'+var
    themis_asi_settings = dictionary($
        'time_range', time_range, $
        'sites', sites, $
        'min_elev', min_elev, $
        'mlat_range', mlat_range, $
        'mlt_range', mlt_range, $
        'var', var )
    default_settings[var] = themis_asi_settings
    if tnames(var) ne '' then return


;---Load data.
    data_file = default_settings.data_file
    if ~cdf_has_var(var, filename=data_file) then begin
        get_data, 'thg_asf_mlt_image', times, mlt_images, limits=lim
        
        ntime = n_elements(times)
        mlt_step = 0.1
        mlt_bins = make_bins(mlt_range, mlt_step)
        nmlt_bin = n_elements(mlt_bins)-1
        mlt_bin_centers = mlt_bins[0:nmlt_bin-1]+mlt_step*0.5
        
        fillval = 1
        ewogram = fltarr(ntime,nmlt_bin)+fillval
        
        npixel = product(lim.image_size)
        mlt_images = reform(mlt_images,[ntime,npixel])>0
        index = where(mlt_images eq 0, count)
        if count ne 0 then mlt_images[index] = !values.f_nan

        pixel_mlt = lim.pixel_mlt
        pixel_mlat = lim.pixel_mlat
        asi_min_count = 5e2
        asi_max_mlt = -2.5
        omega = 24/86400d   ; hour/sec.
        for ii=0, nmlt_bin-1 do begin
            index = where($
                pixel_mlat ge mlat_range[0] and $
                pixel_mlat le mlat_range[1] and $
                pixel_mlt ge mlt_bins[ii] and $
                pixel_mlt le mlt_bins[ii+1], count)
            if count eq 0 then continue
            the_mlts = pixel_mlt[index]

            foreach time, times, time_id do begin
                tmp = mlt_images[time_id,index]
                tindex = where(tmp le asi_min_count and the_mlts le asi_max_mlt+(time-times[0])*omega, count)
                if count ne 0 then tmp[tindex] = fillval
                ewogram[time_id,ii] = mean(tmp,nan=1)
;                ewogram[time_id,ii] = mean(tmp,nan=1)
            endforeach
        endfor
        index = where(finite(ewogram,nan=1), count)
        if count ne 0 then ewogram[index] = fillval

        ewogram_var = 'thg_asf_ewogram'
        store_data, ewogram_var, times, ewogram, mlt_bin_centers
        add_setting, ewogram_var, smart=1, dictionary($
            'display_type', 'spec', $
            'yrange', mlt_range, $
            'ytitle', 'MLT!C(h)', $
            'color_table', 49, $
            'unit', 'Count #', $
            'ztitle', '(Count #)' )
        
        stop
        site_infos = themis_read_mlonimg_default_site_info(sites)
        foreach site, sites, ii do site_infos[ii].min_elev = min_elev
        themis_read_mltimg, time_range, sites=sites, site_infos=site_infos, $
            mlon_range=!null, mlat_range=mlat_range, merge_method=merge_method
        
        stop
        var = themis_calc_mlt_image(time_range, sites=sites, smooth_window=smooth_window, $
            min_elev=min_elev, merge_method=merge_method, mlat_range=mlat_range)
        get_data, var, times, mlt_images, limits=lim
        settings = dictionary(lim)
        time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
        cdf_save_var, time_var, value=times, filename=data_file
        cdf_save_setting, time_settings, varname=time_var, filename=data_file
        settings['depend_0'] = time_var
        cdf_save_var, var, value=mlt_images, filename=data_file
        cdf_save_setting, settings, varname=var, filename=data_file
    endif

    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
    get_data, var, limits=lim
    settings = dictionary(lim)
    image_size = settings.image_size
    npixel = product(image_size)
    foreach key, settings.keys() do begin
        val = settings[key]
        if n_elements(val) eq npixel then val = reform(val,image_size)
        settings[key] = val
    endforeach
    lim = settings.tostruct()
    store_data, var, limits=lim

end

pro load_data_for_2008_0119_lanl, default_settings

;---Settings.
    time_range = default_settings.event_time_range
    probes = ['LANL-'+['01A','02A','04A','97A'],'1989-046','1994-084']
    coord = 'sm'
    
    var = 'lanl'
    if n_elements(e_energy_range) eq 0 then e_energy_range = [0,300d]
    if n_elements(p_energy_range) eq 0 then p_energy_range = [100,1000d]
    default_settings[var] = dictionary($
        'probes', probes, $
        'e_energy_range', e_energy_range, $
        'p_energy_range', p_energy_range, $
        'time_range', time_range )
    data_file = default_settings.data_file

    foreach probe, probes do begin
        prefix = probe+'_'
        r_var = prefix+'r_sm'
        mlt_var = prefix+'mlt'
        e_flux_var = prefix+'kev_e_flux'
        p_flux_var = prefix+'kev_h_flux'
        time_var = 'time_'+probe

    ;---Orbit.
        if ~cdf_has_var(r_var, filename=data_file) then begin
            lanl_read_orbit, time_range, probe=probe, coord=coord
            get_data, r_var, times, r_sm, limits=lim
            settings = dictionary(lim)
            time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
            cdf_save_var, time_var, value=times, filename=data_file
            cdf_save_setting, time_settings, varname=time_var, filename=data_file
            settings['depend_0'] = time_var
            var = r_var
            cdf_save_var, var, value=r_sm, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, r_var, filename=data_file, time_var=time_var, time_type='unix'

        
    ;---MLT.
        if ~cdf_has_var(mlt_var, filename=data_file) then begin
            get_data, r_var, times, r_sm
            mlt = pseudo_mlt(r_sm)
            settings = dictionary($
                'display_type', 'scalar', $
                'depend_0', time_var, $
                'unit', 'hr', $
                'short_name', 'MLT' )
            var = mlt_var
            cdf_save_var, var, value=mlt, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, mlt_var, filename=data_file, time_var=time_var, time_type='unix'

            
    ;---e- and H+ flux.
        e_energy_var = 'e_energy_'+probe
        if ~cdf_has_var(e_flux_var, filename=data_file) then begin
            lanl_read_kev_electron, time_range, probe=probe
            get_data, e_flux_var, times, flux, energys, limits=lim

            var = e_energy_var
            settings = dictionary($
                'unit', 'keV', $
                'data_type', 'metadata' )
            cdf_save_var, var, value=energys, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
            
            var = e_flux_var
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            settings['depend_1'] = e_energy_var
            cdf_save_var, var, value=flux, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, e_flux_var, filename=data_file, time_var=time_var, time_type='unix'
        energys = cdf_read_var(e_energy_var, filename=data_file)
        get_data, e_flux_var, times, fluxs
        store_data, e_flux_var, times, fluxs, energys

        p_energy_var = 'h_energy_'+probe
        if ~cdf_has_var(p_flux_var, filename=data_file) then begin
            lanl_read_kev_proton, time_range, probe=probe
            get_data, p_flux_var, times, flux, energys, limits=lim

            var = p_energy_var
            settings = dictionary($
                'unit', 'keV', $
                'data_type', 'metadata' )
            cdf_save_var, var, value=energys, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
            
            var = p_flux_var
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            settings['depend_1'] = p_energy_var
            cdf_save_var, var, value=flux, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, p_flux_var, filename=data_file, time_var=time_var, time_type='unix'
        energys = cdf_read_var(p_energy_var, filename=data_file)
        get_data, p_flux_var, times, fluxs
        store_data, p_flux_var, times, fluxs, energys

    ;---Make some adjustments.
        rename_var, e_flux_var, to=prefix+'e_flux'
        rename_var, p_flux_var, to=prefix+'p_flux'
        flux_vars = prefix+['e','p']+'_flux'
        time_step = 10d
        common_times = make_bins(time_range, time_step)
        foreach flux_var, flux_vars do begin
            is_eflux = (strsplit(flux_var,'_',extract=1))[-2] eq 'e'
            
            if is_eflux then interp_time, flux_var, common_times
            get_data, flux_var, times, fluxes, energys, limits=lim

            if is_eflux then begin
                energy_range = e_energy_range
            endif else begin
                energy_range = p_energy_range
            endelse
            
            ; Remove some bad H+ channels.
            sopa_p6_energy_range = [400,670d]
            sopa_p12_energy_range = [50,113]
            if not is_eflux then begin
                index = where_pro(energys, ')(', sopa_p6_energy_range, count=count)
                if count ne 0 then begin
                    energys = energys[index]
                    fluxes = fluxes[*,index]
                endif
                index = where_pro(energys, ')(', sopa_p12_energy_range, count=count)
                if count ne 0 then begin
                    energys = energys[index]
                    fluxes = fluxes[*,index]
                endif
            endif
            
            index = where_pro(energys, '[]', energy_range, count=count)
            if count eq 0 then message, 'Inconsistency...'
            store_data, flux_var, times, fluxes[*,index], energys[index], limits=lim
            options, flux_var, 'labels', lim.labels[index]
            options, flux_var, 'colors', lim.colors[index]
            
            yrange = minmax(fluxes[*,index])
            log_yrange = [floor(alog10(yrange[0])),ceil(alog10(yrange[1]))]
            yrange = 10d^log_yrange
            species_str = (is_eflux)? 'e!U-!N': 'H!U+!N'
            add_setting, flux_var, /smart, {$
                display_type: 'list', $
                ylog: 1, $
                yrange: yrange, $
                color_table: 40, $
                unit: '#/cm!U2!N-s-sr-keV', $
                value_unit: 'keV', $
                short_name: species_str+' flux'}
        endforeach
    endforeach

    

end

pro load_data_for_2008_0119_themis_data, default_settings

;---Settings.
    time_range = default_settings.event_time_range+[-1,1]*1800d
    probes = themis_probes()
    coord = 'sm'
    var = 'themis_data'
    energy_range = [0,400d]
    themis_settings = dictionary($
        'time_range', time_range, $
        'energy_range', energy_range, $
        'var', var )
    default_settings[var] = themis_settings


    time_var = 'time_themis_data'
    data_file = default_settings.data_file
    if ~cdf_has_var(time_var, filename=data_file) then begin
        time_step = 3d
        common_times = make_bins(time_range, time_step)
        time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
        cdf_save_var, time_var, value=common_times, filename=data_file
        cdf_save_setting, time_settings, varname=time_var, filename=data_file
    endif
    common_times = cdf_read_var(time_var, filename=data_file)
    

;---Loop through each probe.
    foreach probe, probes do begin
        prefix = 'th'+probe+'_'

    ;---Orbit.
        r_var = themis_read_orbit(time_range, probe=probe, coord=coord, get_name=1)
        if ~cdf_has_var(r_var, filename=data_file) then begin
            r_var = themis_read_orbit(time_range, probe=probe, coord=coord)
            interp_time, r_var, common_times
            get_data, r_var, times, r_sm, limits=lim
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            var = r_var
            cdf_save_var, var, value=r_sm, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, r_var, filename=data_file, time_var=time_var, time_type='unix'

    ;---MLT.
        mlt_var = prefix+'mlt'
        if ~cdf_has_var(mlt_var, filename=data_file) then begin
            get_data, r_var, times, r_sm
            mlt = pseudo_mlt(r_sm)
            settings = dictionary($
                'display_type', 'scalar', $
                'depend_0', time_var, $
                'unit', 'hr', $
                'short_name', 'MLT' )
            var = mlt_var
            cdf_save_var, var, value=mlt, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, mlt_var, filename=data_file, time_var=time_var, time_type='unix'
        

    ;---B field.
        b_var = themis_read_bfield(time_range, probe=probe, coord=coord, get_name=1)
        if ~cdf_has_var(b_var, filename=data_file) then begin
            b_var = themis_read_bfield(time_range, probe=probe, coord=coord)
            interp_time, b_var, common_times
            get_data, b_var, times, b_sm, limits=lim
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            var = b_var
            cdf_save_var, var, value=b_sm, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, b_var, filename=data_file, time_var=time_var, time_type='unix'
        
    ;---keV flux.
        if probe eq 'b' then continue
        e_energy_var = 'e_energy_'+probe
        e_flux_var = themis_read_kev_flux(time_range, probe=probe, no_spec=1, id='e', get_name=1)
;cdf_del_var, e_flux_var, filename=data_file
        if ~cdf_has_var(e_flux_var, filename=data_file) then begin 
            e_flux_var = themis_read_kev_flux(time_range, probe=probe, no_spec=1, id='e')
            interp_time, e_flux_var, common_times
            get_data, e_flux_var, times, flux, energys, limits=lim

            var = e_energy_var
            settings = dictionary($
                'unit', 'keV', $
                'data_type', 'metadata' )
            cdf_save_var, var, value=energys, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
            
            var = e_flux_var
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            settings['depend_1'] = e_energy_var
            cdf_save_var, var, value=flux, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, e_flux_var, filename=data_file, time_var=time_var, time_type='unix'


        p_energy_var = 'p_energy_'+probe
        p_flux_var = themis_read_kev_flux(time_range, probe=probe, no_spec=1, id='p', get_name=1)
;cdf_del_var, p_flux_var, filename=data_file
        if ~cdf_has_var(p_flux_var, filename=data_file) then begin 
            p_flux_var = themis_read_kev_flux(time_range, probe=probe, no_spec=1, id='p')
            interp_time, p_flux_var, common_times
            get_data, p_flux_var, times, flux, energys, limits=lim

            var = p_energy_var
            settings = dictionary($
                'unit', 'keV', $
                'data_type', 'metadata' )
            cdf_save_var, var, value=energys, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
            
            var = p_flux_var
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            settings['depend_1'] = p_energy_var
            cdf_save_var, var, value=flux, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, p_flux_var, filename=data_file, time_var=time_var, time_type='unix'


    ;---Make some adjustments.
        flux_vars = prefix+['e','p']+'_flux'
        foreach flux_var, flux_vars do begin
            get_data, flux_var, times, fluxes, energys, limits=lim
            index = where_pro(energys, '[]', energy_range, count=count)
            if count eq 0 then message, 'Inconsistency...'
            store_data, flux_var, times, fluxes[*,index], energys[index], limits=lim
            options, flux_var, 'labels', lim.labels[index]
            options, flux_var, 'colors', lim.colors[index]
        endforeach
    endforeach

end

pro load_data_for_2008_0119_j_ver_mlt_image, default_settings
    
;---Settings.
    time_range = default_settings.event_time_range
    mlat_range = [50d,90]   ; deg
    var = 'thg_j_ver_mlt_image'
    time_var = 'time_weygand_j'
    j_ver_settings = dictionary($
        'time_range', time_range, $
        'mlat_range', mlat_range, $
        'var', var )
    default_settings[var] = j_ver_settings
    if tnames(var) ne '' then return
    
;---Load data.
    data_file = default_settings.data_file
    if ~cdf_has_var(var, filename=data_file) then begin
        var = themis_read_j_ver_mlt_image(time_range)
        get_data, var, times, mlt_images, limits=lim
        settings = dictionary(lim)
        time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
        cdf_save_var, time_var, value=times, filename=data_file
        cdf_save_setting, time_settings, varname=time_var, filename=data_file
        settings['depend_0'] = time_var
        cdf_save_var, var, value=mlt_images, filename=data_file
        cdf_save_setting, settings, varname=var, filename=data_file
    endif
    
    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
    get_data, var, limits=lim
    settings = dictionary(lim)
    image_size = settings.image_size
    npixel = product(image_size)
    foreach key, settings.keys() do begin
        val = settings[key]
        if n_elements(val) eq npixel then val = reform(val,image_size)
        settings[key] = val
    endforeach
    lim = settings.tostruct()
    store_data, var, limits=lim
end

pro load_data_for_2008_0119_thg_mlt_image, default_settings

;---Settings.
    time_range = default_settings.event_time_range
    sites = ['fsim','gill','inuv']
    min_elev = 2.5  ; deg
    mlat_range = [50d,90]   ; deg
    merge_method = 'merge_elev'
    smooth_window = 2.5*60  ; sec
    var = 'thg_asf_mlt_image'
    time_var = 'time_'+var
    themis_asi_settings = dictionary($
        'time_range', time_range, $
        'sites', sites, $
        'smooth_window', smooth_window, $
        'min_elev', min_elev, $
        'mlat_range', mlat_range, $
        'merge_method', merge_method, $
        'var', var )
    default_settings[var] = themis_asi_settings
    if tnames(var) ne '' then return

    
;---Load data.
    data_file = default_settings.data_file
    if ~cdf_has_var(var, filename=data_file) then begin
        var = themis_calc_mlt_image(time_range, sites=sites, smooth_window=smooth_window, $
            min_elev=min_elev, merge_method=merge_method, mlat_range=mlat_range)
        get_data, var, times, mlt_images, limits=lim
        settings = dictionary(lim)
        time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
        cdf_save_var, time_var, value=times, filename=data_file
        cdf_save_setting, time_settings, varname=time_var, filename=data_file
        settings['depend_0'] = time_var
        cdf_save_var, var, value=mlt_images, filename=data_file
        cdf_save_setting, settings, varname=var, filename=data_file
    endif

    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
    get_data, var, limits=lim
    settings = dictionary(lim)
    image_size = settings.image_size
    npixel = product(image_size)
    foreach key, settings.keys() do begin
        val = settings[key]
        if n_elements(val) eq npixel then val = reform(val,image_size)
        settings[key] = val
    endforeach
    lim = settings.tostruct()
    store_data, var, limits=lim

end


pro load_data_for_2008_0119_po_mlt_image, default_settings

;---Settings.
    time_range = default_settings.event_time_range
    mlat_range = [50d,90]   ; deg
    var = 'po_mlt_image'
    time_var = 'time_'+var
    themis_asi_settings = dictionary($
        'time_range', time_range, $
        'sites', sites, $
        'mlat_range', mlat_range, $
        'var', var )
    default_settings[var] = themis_asi_settings
    if tnames(var) ne '' then return

;---Load data.
    data_file = default_settings.data_file
    if ~cdf_has_var(var, filename=data_file) then begin
        polar_calc_mlt_image, time_range
        get_data, var, times, mlt_images, limits=lim
        settings = dictionary(lim)
        time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
        cdf_save_var, time_var, value=times, filename=data_file
        cdf_save_setting, time_settings, varname=time_var, filename=data_file
        settings['depend_0'] = time_var
        cdf_save_var, var, value=mlt_images, filename=data_file
        cdf_save_setting, settings, varname=var, filename=data_file
    endif

    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
    get_data, var, limits=lim
    settings = dictionary(lim)
    image_size = settings.image_size
    npixel = product(image_size)
    foreach key, settings.keys() do begin
        val = settings[key]
        if n_elements(val) eq npixel then val = reform(val,image_size)
        settings[key] = val
    endforeach
    lim = settings.tostruct()
    store_data, var, limits=lim

end



pro load_data_for_2008_0110_goes, default_setting

    probes = ['11','12']
    coord = 'sm'
    time_range = default_setting.event_time_range
    foreach probe, probes do begin
        prefix = 'g'+probe+'_'
        goes_read_bfield, time_range, probe=probe, coord=coord
        b_var = prefix+'b_'+coord
        get_data, b_var, times, b_sm
        tilt = asin(b_sm[*,2]/snorm(b_sm))*constant('deg')
        tilt_var = prefix+'theta'
        store_data, tilt_var, times, tilt
        add_setting, tilt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'Tilt', $
            'unit', 'deg' )
        
        r_var = goes_read_orbit(time_range, probe=probe, coord=coord)

        mlt_var = prefix+'mlt'
        get_data, r_var, times, r_sm
        mlt = pseudo_mlt(r_sm)
        store_data, mlt_var, times, mlt
        add_setting, mlt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'MLT', $
            'unit', 'hr' )
    endforeach
    
    
end



function load_data_for_2008_0119, energy_range=energy_range

    root_dir = join_path([googledir(),'works','polar_vs_fac'])
    data_dir = join_path([root_dir,'data'])
    data_file = join_path([data_dir,'2008_0119_data.cdf'])
    plot_dir = join_path([root_dir,'plot'])
    if file_test(data_file) eq 0 then cdf_touch, data_file
    if n_elements(energy_range) eq 0 then energy_range = [0d,400]
    cd, root_dir


;---Global attributes.
    event_time_range = time_double(['2008-01-19/06:00','2008-01-19/09:00'])
    aurora_times = time_double('2008-01-19/'+['06:08','06:19','06:28','06:50','07:01','07:14','07:39','08:06'])
    default_settings = dictionary($
        'event_time_range', event_time_range, $
        'aurora_times', aurora_times, $
        'root_dir', root_dir, $
        'data_dir', data_dir, $
        'plot_dir', plot_dir, $
        'data_file', data_file, $
        'energy_range', energy_range )
        

;---LANL.
    load_data_for_2008_0119_lanl, default_settings
    
;---GOES.
    load_data_for_2008_0110_goes, default_settings
    
;---THEMIS data.
    load_data_for_2008_0119_themis_data, default_settings

;---Polar UVI.
    load_data_for_2008_0119_po_mlt_image, default_settings

;---THEMIS ASI.
    load_data_for_2008_0119_thg_mlt_image, default_settings

;---THEMIS current.
    load_data_for_2008_0119_j_ver_mlt_image, default_settings
    
;---Load footpoints.
    load_data_for_2008_0119_footpoint, default_settings
    
;---Load injection data.
    load_data_for_2008_0119_trace, default_settings
    load_data_for_2008_0119_injection, default_settings

;---THEMIS ASI MLT image uniform.
;    load_data_for_2008_0119_themis_uniform_mlt_image_uniform, default_settings


    return, default_settings
    
end

tmp = load_data_for_2008_0119()
stop
vars = ['1989-046','thd','tha','1994-084']+'_mlt'
time_step = 10d
time_range = tmp.event_time_range
common_times = make_bins(time_range, time_step)
foreach var, vars do begin
    get_data, var, times, mlt
    mlt = interpol(mlt, times, common_times)
    store_data, var, common_times, mlt
endforeach
stplot_merge, vars, newname='combo_mlt', $
    colors=sgcolor(['black','blue','red','green']), $
    labels=strupcase(['1989-046','th-d','th-a','1994-084']), $
    ytitle='(hr)'
vars = ['1989-046','thd','tha','1994-084']+'_e_flux'
tplot, [vars,'combo_mlt'], trange=time_range
timebar, tmp.aurora_times, linestyle=1
end