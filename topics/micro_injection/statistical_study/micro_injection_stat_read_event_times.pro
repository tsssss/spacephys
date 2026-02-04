;+
; Read event time ranges for a given time_range.
;-

function micro_injection_stat_read_event_times_gen_file_v01, input_time_range, $
    mission_probe=mission_probe, filename=data_file, $
    errmsg=errmsg, local_root=local_root

    compile_opt idl2
    on_error, 0
    errmsg = ''
    retval = !null

    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = mission+probe+'_'
    search_var = prefix+'mi_search'
    date = time_double(input_time_range[0])
    secofday = constant('secofday')
    date = date-(date mod secofday)
    time_range = date+[0,secofday]
    if check_if_update(search_var, time_range) then begin
    ;---Load data.
        log_flux_vars = list()
        routine = 'micro_injection_read_cwt_kev_electron'
        flux_vars = call_function(routine, time_range, mission_probe=mission_probe, errmsg=errmsg, log=1)
        spec_vars = flux_vars+'_mor'

    ;---Find micro injections.
        freq_range = minmax(1d/([1d,10]*60))    ; in Hz.
        significant_flux_levels = [1.5,0.3]
        significant_psd_levels = 1e-1*[1,1]

    ;---Collect info.
        search_info = orderedhash()
        energy_strs = list()
        foreach flux_var, flux_vars, vid do begin
            energy_str = get_var_setting(flux_var, 'labels')
            energy_strs.add, energy_str
            spec_var = spec_vars[vid]
            search_info[energy_str] = dictionary($
                'flux_var', flux_var, $
                'spec_var', spec_var, $
                'significant_flux_level', significant_flux_levels[vid], $
                'significant_psd_level', significant_psd_levels[vid] )
            options, [flux_var,spec_var], energy_str=energy_str
        endforeach
        energy_strs = energy_strs.toarray()
        flux_unit = get_var_setting(flux_vars[0],'unit')
        spec_unit = get_var_setting(spec_vars[0],'unit')
        search_settings = dictionary($
            'requested_time_range', time_range, $
            'mission', mission, $
            'probe', probe, $
            'time_range', time_range, $
            'energy_strs', energy_strs, $
            'freq_range', freq_range, $
            'freq_range_unit', 'Hz', $
            'spec_unit', spec_unit, $
            'flux_unit', flux_unit )


    ;---Find microinjections.
        foreach energy_str, search_settings['energy_strs'] do begin
            my_info = search_info[energy_str]
            flux_var = my_info['flux_var']
            spec_var = my_info['spec_var']
            significant_flux_level = my_info['significant_flux_level']
            significant_psd_level = my_info['significant_psd_level']
            options, spec_var, constant=freq_range
            options, flux_var, constant=significant_flux_level
            mi_search = micro_injection_use_wavelet_to_find_event(spec_var, flux_var, $
                significant_flux_level, significant_psd_level, freq_range, time_range, $
                test=test, plot_dir=plot_dir, gen_plot=gen_psd_plot)

            mi_times = list()
            my_info['mi_search'] = mi_search
            if n_elements(mi_search) ne 0 then begin
                foreach time, mi_search.keys() do begin
                    info = mi_search[time]
                    if not info['is_mi'] then continue
                    mi_times.add, time
                endforeach
                my_info['mi_times'] = mi_times.toarray()
            endif else begin
                my_info['mi_times'] = []
            endelse
        endforeach

        date = time_range[0]
        search_var = var_store(search_var, search_info, date, search_settings)
        options, search_var, 'requested_time_range', time_range
    endif else begin
        search_info = get_var_data(search_var, search_settings, times=date)
        search_settings = search_settings[0]
    endelse
    
;---Save to file.
    cdf_save_setting, filename=data_file, search_settings    
    foreach energy_str, search_info.keys() do begin
        prefix = energy_str+'_'
        my_info = search_info[energy_str]
        settings = dictionary()
        foreach key, my_info.keys() do begin
            val = my_info[key]
            if size(val,type=1) eq 11 or n_elements(val) gt 1 then continue
            settings[key] = val
        endforeach
        mi_search = my_info['mi_search']
        if size(mi_search,type=1) ne 11 then begin
            ; No time of interest.
            times = !null
        endif else begin
            times = (mi_search.keys()).toarray()
        endelse

        ntime = n_elements(times)        
        if ntime eq 0 then begin
            times = date
            flags = [!values.f_nan]
            msgs = ['']
        endif else begin
            flags = intarr(ntime)
            msgs = strarr(ntime)
            foreach time, times, tid do begin
                tmp = mi_search[time]
                flags[tid] = tmp['is_mi']
                msgs[tid] = strjoin((tmp['msgs']).toarray(),'%')
            endforeach
        endelse
        searched_time_var = var_store(prefix+'searched_time_flags', flags, times)
        
        mi_index = where(flags eq 1, nmi_time)
        if nmi_time eq 0 then begin
            mi_times = date
            mi_freqs = [!values.f_nan]
            mi_psds = [!values.f_nan]
            mi_psd_ratios = [!values.f_nan]
        endif else begin
            mi_times = times[mi_index]
            mi_freqs = fltarr(nmi_time)
            mi_psds = fltarr(nmi_time)
            mi_psd_ratios = fltarr(nmi_time)
            foreach time, mi_times, tid do begin
                tmp = mi_search[time]
                mi_freqs[tid] = tmp['freq']
                mi_psds[tid] = tmp['psd']
                mi_psd_ratios[tid] = tmp['psd_ratio']
            endforeach
        endelse
        mi_freq_var = var_store(prefix+'mi_freq', mi_freqs, mi_times, settings=dictionary('unit','Hz'))
        mi_psd_var = var_store(prefix+'mi_psd', mi_psds, mi_times, settings=dictionary('unit',search_settings['spec_unit']))
        mi_psd_ratio_var = var_store(prefix+'mi_psd_ratio', mi_psd_ratios, mi_times)
         
        vars = [mi_freq_var,mi_psd_var,mi_psd_ratio_var]
        print, vars
        if n_elements(vars) eq 0 then stop
        stplot2cdf, vars, time_var=prefix+'mi_times', filename=data_file
        
        vars = [searched_time_var]
        stplot2cdf, vars, time_var=prefix+'searched_times', filename=data_file
    endforeach

    
    return, data_file

end


function micro_injection_stat_read_event_times, time_range, mission_probe=mission_probe, id=method_id, version=version

    if n_elements(version) eq 0 then version = 'v01'
    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = mission+probe+'_'
    mission_probe_str = mission+'_'+probe

    project_info = micro_injection_stat_load_project()
    data_dir = project_info['data_dir']

;---Load files.
    base_name = mission_probe_str+'_event_times_%Y_%m%d_'+version+'.cdf'
    local_root = data_dir
    local_path = [local_root,'event_times',mission_probe_str,'%Y']
    request = dictionary($
        'pattern', dictionary($
        'local_file', join_path([local_path,base_name]), $
        'local_index_file', join_path([local_path,default_index_file()])), $
        'cadence', 'day' )

    files = prepare_files(request=request, errmsg=errmsg, $
        file_times=file_times, time=time_range, nonexist_files=nonexist_files)
    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            routine = 'micro_injection_stat_read_event_times_gen_file_'+version
            file = call_function(routine, file_time, mission_probe=mission_probe, filename=local_file)
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif

;---Read mi_times, the simplest is to combine both energy bins' mi_times.
    gatt = cdf_read_setting(filename=files[0])
    energy_strs = gatt['energy_strs']
    vars = energy_strs+'_mi_times'
    mi_times = []
    foreach var, vars do begin
        times = []
        foreach file, files do begin
            times = [times, cdf_read_var(var, filename=file)]
        endforeach
        mi_times = [mi_times, times]
    endforeach
    if n_elements(mi_times) eq 0 then return, []
    index = where_pro(mi_times, '[]', time_range, count=count)
    if count eq 0 then return, []
    time_step = project_info['common_time_step']
    mi_trs = time_to_range(sort_uniq(mi_times[index]), time_step=time_step)
    
    return, mi_trs
    
end

test = 0
missions = ['mms']

search_trs = micro_injection_load_search_time_range()
nsearch_tr = n_elements(search_trs[*,0])
dates = list()
for ii=0,nsearch_tr-1 do begin
    tr = reform(search_trs[ii,*])
    nday = total(tr*[-1,1])/constant('secofday')
    dates.add, findgen(nday)*constant('secofday')+tr[0], extract=1
endfor

foreach date, dates do begin
if time_double(date) lt time_double('2016-01-21') then continue
    time_range = date+[0,constant('secofday')]
    foreach mission, missions do begin
        probes = call_function(mission+'_probes')
        foreach probe, probes do begin
            mission_probe = [mission,probe]
            mission_probe_str = mission+'_'+probe
            del_data, '*'
            mi_trs = micro_injection_stat_read_event_times(time_range, mission_probe=mission_probe)
            print, 'Processed '+mission_probe_str+' for '+time_string(date)+' ...'
        endforeach
    endforeach
endforeach

stop



tr = ['2015-10-01','2015-10-02']
mission_probe = ['mms','1']
file = join_path([homedir(),'test_mi_search.cdf'])
if file_test(file) then file_delete, file
print, micro_injection_stat_read_event_times_gen_file_v01(tr, mission_probe=mission_probe, filename=file)
stop
tmp = micro_injection_stat_read_event_times(tr, mission_probe=mission_probe)
end