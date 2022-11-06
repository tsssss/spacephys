;+
; Read training data for HOPE ML.
; 
; nan_duration=. Default is 20 min. Smaller gaps will be eliminated.
;-
pro ml_hope_load_training_data_gen_file, input_date, filename=out_file, resolution=time_step, nan_duration=nan_duration

    if n_elements(nan_duration) eq 0 then nan_duration = 20*60    ; 20 min.

    if n_elements(input_date) eq 0 then begin
        errmsg = 'No input time ...'
        return
    endif
    if n_elements(time_step) eq 0 then time_step = ml_time_step()

    year_str = time_string(input_date[0],tformat='YYYY')
    year = float(year_str)
    time_range = time_double([year_str,string(year+1,format='(I04)')])
    secofday = constant('secofday')
if keyword_set(test) then time_range = time_range[0]+[0,30*secofday]

;---Init file.
    if file_test(out_file) eq 0 then cdf_touch, out_file

;---Add time.
    time_var = 'unix_time'
    if ~cdf_has_var(time_var, filename=out_file) then begin
        common_times = make_bins(time_range+[0,-1]*time_step,time_step)+time_step*0.5
        cdf_save_var, time_var, value=common_times, filename=out_file
        settings = dictionary('UNITS', 'sec', 'VAR_TYPE', 'meta_data')
        cdf_save_setting, settings, varname=time_var, filename=out_file
    endif
    common_times = cdf_read_var(time_var, filename=out_file)
    
    
;---Add RBSP data.
    probes = ['a','b']
    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
    
    ;---Flux.
        all_species = rbsp_hope_species()
        foreach species, all_species do begin
            var = ml_rbsp_hope_read_en_spec(time_range, probe=probe, species=species)
            if cdf_has_var(var, filename=out_file) then continue
            data = get_var_data(var, energys, limits=lim)
            energy_var = species+'_energy'
            if ~cdf_has_var(energy_var, filename=out_file) then begin
                cdf_save_var, energy_var, value=energys, filename=out_file
                settings = dictionary('UNITS', 'eV', 'VAR_TYPE', 'meta_data')
                cdf_save_setting, settings, filename=out_file, varname=energy_var
            endif
            cdf_save_var, var, value=data, filename=out_file
            settings = dictionary(lim)
            settings['DEPEND_0'] = time_var
            settings['DEPEND_1'] = energy_var
            cdf_save_setting, settings, filename=out_file, varname=var
        endforeach
        
    ;---Orbit vars.
        routines = 'ml_rbsp_read_'+['lshell','mlt','mlat']
        foreach routine, routines do begin
            var = call_function(routine, time_range, probe=probe)
            if cdf_has_var(var, filename=out_file) then continue
            data = get_var_data(var, limits=lim)
            cdf_save_var, var, value=data, filename=out_file
            settings = dictionary(lim)
            settings['DEPEND_0'] = time_var
            cdf_save_setting, settings, filename=out_file, varname=var
        endforeach
    endforeach

;---Add OMNI data.
    file = ml_omni_load_param(time_range)
    omni_time_var = 'unix_time'
    times = cdf_read_var(omni_time_var, filename=file)
    foreach var, cdf_vars(file) do begin
        if var eq omni_time_var then continue
        if cdf_has_var(var, filename=out_file) then continue
        data = cdf_read_var(var, filename=file)
        ndim = size(data,n_dimensions=1)
        data_mag = (ndim eq 1)? data: snorm(data)
        nan_index = where(finite(data_mag,nan=1), count, complement=good_index)
        if count ne 0 then begin
            ; Interpolate to remove nan.
            data = sinterpol(data[good_index,*], times[good_index], times)
            ; Keep nan sections that are long.
            nan_time_ranges = times[time_to_range(nan_index, time_step=1)]
            nnan_time_range = n_elements(nan_time_ranges)*0.5
            durations = nan_time_ranges[*,1]-nan_time_ranges[*,0]
            index = where(durations ge nan_duration, nnan_section)
            if nnan_section eq 0 then begin
                nan_sections = !null
            endif else begin
                nan_section_times = nan_time_ranges[index,*]
                nan_sections = (nan_section_times-times[0])/time_step
            endelse
            for ii=0,nnan_section-1 do begin
                i0 = nan_sections[ii,0]
                i1 = nan_sections[ii,1]
                data[i0:i1,*] = !values.f_nan
            endfor
        endif        
        cdf_save_var, var, value=data, filename=out_file
        settings = cdf_read_setting(var, filename=file)
        settings['DEPEND_0'] = time_var
        cdf_save_setting, settings, filename=out_file, varname=var
    endforeach

end