;+
; Read DF for given probe, time_range, mlt_range, and rxy_range.
;
; time_range.
; probe=.
; mlt_range=.
; rxy_range=.
; project=.
;-

function azim_df_read_df, time_range, probe=probe, project=project, $
    mlt_range=mlt_range, rxy_range=rxy_range

    retval = list()
    if n_elements(probe) eq 0 then return, retval
    if n_elements(project) eq 0 then project = azim_df_load_project()
    cdf_file = join_path([project.data_dir, 'azim_df_detected_df_'+probe+'.cdf'])
    log_file = join_path([project.data_dir, 'azim_df_detected_df_'+probe+'.log'])
    txt_file = join_path([project.data_dir, 'azim_df_detected_df_'+probe+'.txt'])
    if file_test(cdf_file) eq 0 then begin
        azim_df_detect_df, probe, log_file=log_file, data_file=cdf_file, txt_file=txt_file
    endif

    in_vars = ['obs_rxy','obs_mlt','obs_r_sm','height', $
        'scaled_height','width','time_range','theta_range']
    prefix = probe+'_'
    out_vars = prefix+in_vars
    del_data, out_vars
    request = dictionary($
        'var_list', list($
            dictionary($
                'in_vars', in_vars, $
                'out_vars', out_vars, $
                'time_var_name', 'obs_time', $
                'time_var_type', 'unix')))
    read_files, time_range, files=cdf_file, request=request, errmsg=errmsg

    retval = list()
    the_var = out_vars[0]
    data = get_var_data(the_var)
    if n_elements(data) eq 0 then return, retval

    if n_elements(mlt_range) eq 2 then begin
        the_range = minmax(mlt_range)
        the_var = prefix+'obs_mlt'
        data = get_var_data(the_var)
        index = where_pro(data, '()', the_range, count=count)
        if count eq 0 then return, retval
        foreach var, out_vars do begin
            get_data, var, time, data
            store_data, var, time[index], data[index,*]
        endforeach
    endif

    if n_elements(rxy_range) eq 2 then begin
        the_range = minmax(rxy_range)
        the_var = prefix+'obs_rxy'
        data = get_var_data(the_var)
        index = where_pro(data, '()', the_range, count=count)
        if count eq 0 then return, retval
        foreach var, out_vars do begin
            get_data, var, time, data
            store_data, var, time[index], data[index,*]
        endforeach
    endif


    df_list = list()
    get_data, out_vars[0], obs_times
    if n_elements(obs_times) eq 1 then begin
        if obs_times[0] eq 0 then return, df_list
        the_df = dictionary($
            'probe', probe, $
            'obs_time', obs_times[0])
        foreach key, in_vars, var_id do begin
            the_data = get_var_data(out_vars[var_id])
            the_df[key] = (n_elements(the_data) eq 1)? the_data[0]: reform(the_data)
        endforeach
        df_list.add, the_df
    endif else begin
        data_list = list()
        data_dims = list()
        foreach var, out_vars do begin
            data = get_var_data(var)
            data_list.add, data
            data_dims.add, size(data,/n_dimension)
        endforeach
        foreach obs_time, obs_times, ii do begin
            the_df = dictionary($
                'probe', probe, $
                'obs_time', obs_time)
            foreach key, in_vars, var_id do begin
                the_df[key] = (data_dims[var_id] eq 1)? (data_list[var_id])[ii]: reform((data_list[var_id])[ii,*])
            endforeach
            df_list.add, the_df
        endforeach
    endelse

    return, df_list

end

;time_range = time_double(['2014-08-28/09:30','2014-08-28/11:00'])
time_range = time_double(['2009','2014'])
project = azim_df_load_project()
foreach probe, project.all_probes do df_list = azim_df_read_df(time_range, probe=probe)
end
