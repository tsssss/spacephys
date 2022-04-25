;+
; Generate 1-min resolution HOPE integrated flux for e, p, o, he.
; Perigee data are excluded. Dates when HOPE data are missing are filled with NaNs.
;-
pro ml_rbsp_hope_load_en_spec_gen_file, input_date, probe=probe, $
    filename=out_file, resolution=time_step, downsample_method=downsample_method

test = 0

    if n_elements(input_date) eq 0 then begin
        errmsg = 'No input time ...'
        return
    endif
    if n_elements(probe) eq 0 then begin
        errmsg = 'No input probe ...'
        return
    endif
    if n_elements(time_step) eq 0 then time_step = ml_time_step()
    ; can be mean, median, log_mean.
    if n_elements(downsample_method) eq 0 then downsample_method = 'mean'

    year_str = time_string(input_date[0],tformat='YYYY')
    year = float(year_str)
    time_range = time_double([year_str,string(year+1,format='(I04)')])
    secofday = constant('secofday')
if keyword_set(test) then time_range = time_range[0]+[0,30*secofday]
if keyword_set(test) then time_range = time_double(['2015-03-15','2015-03-18'])
    days = make_bins(time_range+[0,-1]*secofday,secofday)
    common_times = make_bins(time_range+[0,-1]*time_step,time_step)+time_step*0.5

    prefix = 'rbsp'+probe+'_'
    all_species = rbsp_hope_species()

;---Init file.
    if file_test(out_file) eq 0 then cdf_touch, out_file

;---Add time.
    time_var = 'unix_time'
    if ~cdf_has_var(time_var, filename=out_file) then begin
        common_times = make_bins(time_range+[0,-1]*time_step, time_step)+time_step*0.5

        cdf_save_var, time_var, value=common_times, filename=out_file
        settings = dictionary($
            'FIELDNAM', 'Time', $
            'UNITS', 'sec', $
            'TEXT', 'center time of the sampling interval', $
            'VAR_TYPE', 'support_data' )
        cdf_save_setting, settings, var=time_var, filename=out_file
    endif
    common_times = cdf_read_var(time_var, filename=out_file)
    ntime = n_elements(common_times)
    
;---Add energy bins.
    foreach species, all_species do begin
        energy_var = 'energy_bins_'+species
        energy_bins = ml_rbsp_hope_read_energy_bins(probe=probe, species=species)
        if ~cdf_has_var(energy_var, filename=out_file) then begin
            cdf_save_var, energy_var, value=energy_bins, filename=out_file, save_as_one=1
            settings = dictionary($
                'FIELDNAM', 'Energy', $
                'UNITS', 'eV', $
                'VAR_TYPE', 'support_data')
            cdf_save_setting, settings, var=energy_var, filename=out_file
        endif
    endforeach

;---Add energy spectrogram.
    fillval = 0.001
    files = rbsp_load_hope(time_range, probe=probe, id='l3%pa')
    foreach species, all_species do begin
        energy_var = 'energy_bins_'+species
        energy_bins = cdf_read_var(energy_var, filename=out_file)
        nenergy = n_elements(energy_bins)
        energy_range = minmax(energy_bins)
        energy_step = mean(energy_bins[1:nenergy-1]/energy_bins[0:nenergy-2])
        energy_edges = sqrt(energy_bins[1:nenergy-1]*energy_bins[0:nenergy-2])
        energy_edges = [energy_edges[0]/energy_step,energy_edges,energy_edges[-1]*energy_step]
        en_spec_var = prefix+'hope_flux_'+species
        flux_var = strupcase('f'+species+'do')
        flux_lowres = fltarr(ntime,nenergy)+fillval

        settings = dictionary($
            'FIELDNAM', species+' flux integrated over pitch angle and gyro angle', $
            'UNITS', '#/s-cm!U-2!N-sr-keV', $
            'VAR_TYPE', 'data', $
            'DEPEND_0', time_var, $
            'DEPEND_1', energy_var, $
            'DISPLAY_TYPE', 'spectrogram')

        var_list = list()
        the_time_var = (species eq 'e')? 'Epoch_Ele': 'Epoch_Ion'
        energy_var = (species eq 'e')? 'HOPE_ENERGY_Ele': 'HOPE_ENERGY_Ion'
        var_list.add, dictionary($
            'in_vars', flux_var, $
            'out_vars', en_spec_var, $
            'time_var_name', the_time_var, $
            'time_var_type', 'epoch' )
        var_list.add, dictionary($
            'in_vars', energy_var, $
            'time_var_name', the_time_var, $
            'time_var_type', 'epoch' )
        read_vars, time_range, files=files, var_list=var_list
        get_data, en_spec_var, times, flux, limits=lim, dlimits=dlim
        get_data, energy_var, times, energys
        index = where(flux le 0 or flux ge 1e30, count)
        if count ne 0 then flux[index] = fillval

        ; Calculate the lowres data.
        foreach time, common_times, time_id do begin
            time_index = where(times ge time-0.5*time_step and $
                times le time+0.5*time_step, count)
            if count eq 0 then continue
            the_energy = energys[time_index,*]
            the_flux = flux[time_index,*]
            foreach energy, energy_bins, energy_id do begin
                index = where(the_energy ge energy_edges[energy_id] and $
                    the_energy le energy_edges[energy_id+1], count)
                if count eq 0 then continue
                if downsample_method eq 'mean' then begin
                    flux_lowres[time_id,energy_id] = mean(the_flux[index], nan=1)
                endif else if downsample_method eq 'median' then begin
                    flux_lowres[time_id,energy_id] = median(the_flux[index], nan=1)
                endif else if downsample_method eq 'log_mean' then begin
                    flux_lowres[time_id,energy_id] = 10^mean(alog10(the_flux[index]), nan=1)
                endif
            endforeach
        endforeach
                
        if keyword_set(test) then begin
            lowres_en_spec_var = en_spec_var+'_lowres'
            store_data, lowres_en_spec_var, common_times, flux_lowres, energy_bins, limits=lim, dlimits=dlim
            store_data, en_spec_var, times, flux, energys
            cdf2tplot, file
            vars = [en_spec_var, lowres_en_spec_var]
            zlim, vars, 1e4, 1e8, 1
            ylim, vars, 4, 4e4, 1
            options, vars, 'spec', 1
            options, vars, 'no_interp', 1
            sgopen, 0, xsize=800, ysize=450
            tplot, vars, trange=time_range[0]+[0,secofday]
            stop
        endif

        cdf_save_var, en_spec_var, value=flux_lowres, filename=out_file
        cdf_save_setting, settings, var=en_spec_var, filename=out_file
    endforeach

end


;time_range = time_double(['2013-01-01','2019-11-01'])
;;time_range = time_double(['2015-01-01','2020-01-01'])
;;time_range = time_double(['2015-01-01','2019-11-01'])
;secofday = constant('secofday')
;dates = make_bins(time_range+[0,-1]*secofday, secofday)
;foreach date, dates do begin
;    foreach probe, ['a','b'] do begin
;        rbspx = 'rbsp'+probe
;        year = time_string(date, tformat='YYYY')
;        time_str = time_string(date, tformat='YYYY_MMDD')
;        file = join_path([diskdir('data'),'sdata','rbsp','ml_hope_flux',$
;            rbspx,year,rbspx+'_ml_hope_flux_'+time_str+'_v01.cdf'])
;if file_test(file) eq 1 then file_delete, file
;        ml_rbsp_hope_load_en_spec_gen_file, date, probe=probe, filename=file
;    endforeach
;endforeach
;
;stop

date = '2015-03-17'
probe = 'a'
file = join_path([homedir(),'test_ml_hope.cdf'])
file_delete, file, allow_nonexistent=1
ml_rbsp_hope_load_en_spec_gen_file, date, probe=probe, filename=file


end