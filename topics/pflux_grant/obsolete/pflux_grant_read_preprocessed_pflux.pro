;+
; Read preprocessed dE and dB.
; Save output data in rbspx_de_fac, rbspx_dedot0_fac, rbspx_db_fac.
;-
pro pflux_grant_read_preprocessed_pflux, time, probe=probe, filename=data_file, errmsg=errmsg, local_root=local_root, project=project

    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','rbsp'])
    if n_elements(probe) ne 1 then message, 'Invalid probe ...'
    if n_elements(version) eq 0 then version = 'v01'
    index = where(['a','b'] eq probe, count)
    if count eq 0 then message, 'Invalid probe ...'
    rbspx = 'rbsp'+probe
    prefix = rbspx+'_'
    local_path = [local_root,rbspx,'ebfield','%Y']
    base_name = rbspx+'_ebfield_%Y_%m%d_'+version+'.cdf'
    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project.pflux_calc_setting
    valid_range = pflux_calc_setting[rbspx].time_range+[-1,1]*constant('secofday')

    request = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', valid_range, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', prefix+['pf_fac_norm','pfdot0_fac_norm','bw_ratio'], $
                'time_var_name', 'ut_sec', $
                'time_var_type', 'unix'), $
            dictionary($
                'in_vars', prefix+['cnorm'], $
                'time_var_name', 'ut_orbit', $
                'time_var_type', 'unix')))

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)


;---Read data from files and save to memory.
    read_files, time, files=files, request=request
    
    
;---Remove uncertain pfdot0 data using bw_ratio.
    min_bw_ratio = 0.2
    bw_ratio_var = prefix+'bw_ratio'
    get_data, bw_ratio_var, times, bw_ratio
    index = where(abs(bw_ratio) le min_bw_ratio, count)
    if count ne 0 then begin    
        if n_elements(project) eq 0 then project = pflux_grant_load_project()
        pflux_calc_setting = project['pflux_calc_setting']
        pad_time = max(pflux_calc_setting.period_range)
           
        efield_time_step = 1d/16
        time_ranges = time_to_range(times[index], time_step=efield_time_step, pad_time=pad_time)
        ntime_range = n_elements(time_ranges)/2
        
        fillval = !values.f_nan
        pfdot0_var = prefix+'pfdot0_fac_norm'
        get_data, pfdot0_var, times, pfdot0
        for ii=0, ntime_range-1 do begin
            index = where_pro(times, '[]', time_ranges[ii,*], count=count)
            if count eq 0 then continue
            pfdot0[index,*] = fillval
        endfor
        store_data, pfdot0_var, times, pfdot0
    endif
    
    
;---Add settings.
    fac_labels = ['earth','west','out']
    unit = 'mW/m!U2!N'
    var = prefix+'pf_fac_norm'
    add_setting, var, /smart, {$
        display_type: 'vector', $
        unit: unit, $
        short_name: 'S', $
        coord: '', $
        coord_labels: fac_labels}
        
    var = prefix+'pfdot0_fac_norm'
    add_setting, var, /smart, {$
        display_type: 'vector', $
        unit: unit, $
        short_name: 'S!S!Udot0!R!N', $
        coord: '', $
        coord_labels: fac_labels}

end

time_range = time_double(['2013-06-07','2013-06-08'])
time_range = time_double(['2013-05-01','2013-05-02'])
probe = 'b'
pflux_grant_read_preprocessed_pflux, time_range, probe=probe
end