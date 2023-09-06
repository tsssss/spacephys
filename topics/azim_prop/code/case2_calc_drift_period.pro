;+
; Calculate the drift period as a function of time and energy.
;-

id = '2014_0828_10'
info_var = id+'_event_info'
if tnames(info_var) eq '' then _2014_0828_10_load_data
time_range = time_double(['2014-08-28/10:00','2014-08-28/10:40'])

foreach species, ['e'] do begin
;foreach species, ['e','h'] do begin
    energy_range = (species eq 'e')? [50,500]: [70,500]
    probes = (species eq 'e')? ['thd','g15','1991-080','g13','1994-084']: $
        ['thd','LANL-01A','LANL-02A','LANL-04A','LANL-97A']
    foreach probe, probes do begin
        prefix = probe+'_'
        prefix2 = prefix+species+'_'
        drift_period_var = prefix2+'drift_period'
        if tnames(drift_period_var) ne '' then continue
        
        the_var = prefix+'kev_'+species+'_flux'
        get_data, the_var, times, fluxes, energies
        ; Filter energy.
        index = where_pro(energies,'[]', energy_range)
        energies = energies[index]
        fluxes = fluxes[*,index]
        ; Filter time.
        index = where_pro(times, '[]', time_range)
        fluxes = fluxes[index,*]
        times = times[index]
        ; Select times to calc drift period.
        times = time_range
        ntime = n_elements(times)
        nenergy = n_elements(energies)
        drift_periods = fltarr(ntime,nenergy)
        bounce_periods = fltarr(ntime,nenergy)
        lshells = fltarr(ntime,nenergy)
        
        pos_var = prefix+'r_gsm'
        r_gsms = get_var_data(pos_var, at=times)
        foreach time, times, ii do begin
            foreach energy, energies, jj do begin
                drift_periods[ii,jj] = calc_drift_period(rgsm=r_gsms[ii,*], species=species, $
                    time=time, energy=energy, bounce_period=bounce_period, lshell=lshell)
                bounce_periods[ii,jj] = bounce_period
                lshells[ii,jj] = lshell
            endforeach
        endforeach
        
        store_data, drift_period_var, times, drift_periods, energies
        add_setting, drift_period_var, /smart, {$
            display_type: 'list', $
            color_table: 52, $
            unit: 'sec', $
            value_unit: 'keV', $
            short_name: 'T drift'}
        store_data, prefix2+'bounce_period', times, bounce_periods, energies
        add_setting, prefix2+'bounce_period', /smart, {$
            display_type: 'list', $
            color_table: 52, $
            unit: 'sec', $
            value_unit: 'keV', $
            short_name: 'T bounce'}
        store_data, prefix2+'lshell', times, lshells, energies
        add_setting, prefix2+'lshell', /smart, {$
            display_type: 'list', $
            color_table: 52, $
            unit: '#', $
            value_unit: 'keV', $
            short_name: 'L'}
    endforeach
    
    test_time_range = time_double(['2014-08-28/10:00','2014-08-28/10:10'])
    nprobe = n_elements(probes)
    mlts = fltarr(nprobe)
    dmlt = 0.1
    test_mlts = make_bins([-0.5,4.5],dmlt)
    test_xcors = dblarr(n_elements(test_mlts))
    dmlts = fltarr(nprobe)
    foreach probe, probes, probe_id do begin
        foreach test_mlt, test_mlts, jj do begin
            prefix = probe+'_'
            prefix2 = prefix+species+'_'
            drift_period_var = prefix2+'drift_period'
            
            ; Get the MLT.
            r_gsm = get_var_data(prefix+'r_gsm', at=test_time_range)
            r_sm = cotran(r_gsm, test_time_range, 'gsm2sm')
            mlts[probe_id] = mean(pseudo_mlt(r_sm))

            the_var = prefix+'kev_'+species+'_flux'
            get_data, the_var, times, fluxes, energies, limits=lim
            ; Filter energy.
            index = where_pro(energies,'[]', energy_range)
            energies = energies[index]
            fluxes = fluxes[*,index]
            ; Filter time.
            index = where_pro(times, '[]', time_range)
            fluxes = fluxes[index,*]
            times = times[index]
            
            drift_periods = get_var_data(drift_period_var, at=times)
            new_fluxes = fluxes
            foreach energy, energies, ii do begin
                dtimes = test_mlt/24*drift_periods[*,ii]
                new_times = times-dtimes
                new_fluxes[*,ii] = interpol(fluxes[*,ii], new_times, times)
            endforeach
            new_var = prefix+'kev_'+species+'_new_flux'
            store_data, new_var, times, new_fluxes, energies
            
            yrange = alog10(minmax(new_fluxes))
            yrange = 10.^[floor(yrange[0]),ceil(yrange[1])]
            add_setting, new_var, /smart, {$
                display_type: 'list', $
                ylog: 1, $
                yrange: yrange, $
                color_table: 52, $
                unit: '#/cm!U2!N-s-sr-keV', $
                value_unit: 'keV', $
                short_name: 'e!U-!N flux'}
            
            tplot, [the_var, new_var]
            
            new_fluxes = get_var_data(new_var, in=test_time_range)
            xcor = 0
            foreach energy, energies, ii do begin
                if ii eq 0 then continue
                xcor += c_correlate(new_fluxes[*,ii-1],new_fluxes[*,ii],0)
            endforeach
            test_xcors[jj] = xcor/(n_elements(energies)-1)
        endforeach
        max_xcor = max(test_xcors, index, /nan)
        test_mlt = test_mlts[index]
        
        prefix = probe+'_'
        prefix2 = prefix+species+'_'
        drift_period_var = prefix2+'drift_period'

        the_var = prefix+'kev_'+species+'_flux'
        get_data, the_var, times, fluxes, energies, limits=lim
        ; Filter energy.
        index = where_pro(energies,'[]', energy_range)
        energies = energies[index]
        fluxes = fluxes[*,index]
        ; Filter time.
        index = where_pro(times, '[]', time_range)
        fluxes = fluxes[index,*]
        times = times[index]

        drift_periods = get_var_data(drift_period_var, at=times)
        new_fluxes = fluxes
        foreach energy, energies, ii do begin
            dtimes = test_mlt/24*drift_periods[*,ii]
            new_times = times-dtimes
            new_fluxes[*,ii] = interpol(fluxes[*,ii], new_times, times)
        endforeach
        new_var = prefix+'kev_'+species+'_new_flux'
        store_data, new_var, times, new_fluxes, energies
        
        yrange = alog10(minmax(new_fluxes))
        yrange = 10.^[floor(yrange[0]),ceil(yrange[1])]
        add_setting, new_var, /smart, {$
            display_type: 'list', $
            ylog: 1, $
            yrange: yrange, $
            color_table: 52, $
            unit: '#/cm!U2!N-s-sr-keV', $
            value_unit: 'keV', $
            short_name: 'e!U-!N flux'}

        tplot, [the_var, new_var]
        dmlts[probe_id] = test_mlt
    endforeach
    
    
    test_time = time_double('2014-08-28/10:08')
    new_mlts = fltarr(nprobe)
    foreach probe, probes, probe_id do begin
        prefix = probe+'_'
        prefix2 = prefix+species+'_'
        
        new_var = prefix+'kev_'+species+'_new_flux'
        get_data, new_var, new_times, new_fluxes, energies
        ; Get drift period.
        drift_period_var = prefix2+'drift_period'
        drift_periods = get_var_data(drift_period_var, at=new_times)
        tmp = min(new_times-test_time, time_index, /absolute)
        new_time = new_times[time_index]
        
        nenergy = n_elements(energies)
        the_mlts = fltarr(nenergy)
        pos_var = prefix+'r_gsm'
        dmlt = dmlts[probe_id]
        foreach energy, energies, energy_id do begin
            dtime = dmlt/24*drift_periods[time_index,energy_id]
            the_time = new_time+dtime
            r_gsm = get_var_data(pos_var, at=the_time)
            r_sm = cotran(r_gsm, the_time, 'gsm2sm')
            the_mlts[energy_id] = pseudo_mlt(r_sm)
        endforeach
        new_mlts[probe_id] = mean(the_mlts)-dmlt
    endforeach
    
    stop
    
endforeach

end