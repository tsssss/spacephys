;+
; Read ramps for given probe, time_range, mlt_range, and rxy_range.
;
; time_range.
; probe=.
; mlt_range=.
; rxy_range=.
; project=.
;-

pro azim_df_read_ramp_gen_cdf, filename=cdf_file, txt_files=txt_files

    lprmsg, 'Saving ramps to CDF'
    lprmsg, '    From: '+txt_files
    lprmsg, '    To: '+cdf_file


;---Read files to get all ramps.
    ramp_list = list()
    foreach file, txt_files do begin
        lines = read_all_lines(file)
        nline = n_elements(lines)
        foreach line, lines do begin
            the_ramp = azim_df_ramp_read(line)
            ramp_list.add, the_ramp
        endforeach
    endforeach

;---Sort by obs_time.
    nramp = ramp_list.length
    df_obs_times = dblarr(nramp)
    foreach ramp, ramp_list, df_id do df_obs_times[df_id] = ramp.obs_time
    index = sort(df_obs_times)
    ramp_list = ramp_list[index]

;---Collect data.
    df_time_ranges = dblarr(nramp,2)
    df_theta_ranges = fltarr(nramp,2)
    df_widths = fltarr(nramp)
    df_heights = fltarr(nramp)
    df_obs_r_sms = fltarr(nramp,3)
    foreach ramp, ramp_list, df_id do begin
        df_obs_times[df_id] = ramp.obs_time
        df_time_ranges[df_id,*] = ramp.time_range
        df_theta_ranges[df_id,*] = ramp.value_range
        df_widths[df_id] = ramp.width
        df_heights[df_id] = ramp.height
        df_obs_r_sms[df_id,*] = ramp.obs_r_sm
    endforeach

    df_obs_mlts = pseudo_mlt(df_obs_r_sms)
    df_obs_rxys = snorm(df_obs_r_sms[*,0:1])
    project = azim_df_load_project()
    scale_width = project.scale_width
    df_scaled_heights = azim_df_scale_theta(df_heights, df_obs_mlts, width=scale_width)

;---Save data to cdf.
    cdfid = cdf_create(cdf_file)
    time_var = 'obs_time'
    settings = dictionary($
        'unit', 'sec', $
        'fieldnam', 'observed time')
    cdf_save_var, time_var, value=df_obs_times, filename=cdfid
    cdf_save_setting, settings, filename=cdfid, varname=time_var

    the_var = 'time_range'
    settings = dictionary($
        'depend_0', time_var, $
        'unit', 'sec', $
        'fieldnam', 'time range')
    cdf_save_var, the_var, value=df_time_ranges, filename=cdfid
    cdf_save_setting, settings, filename=cdfid, varname=the_var

    the_var = 'theta_range'
    settings = dictionary($
        'depend_0', time_var, $
        'unit', 'deg', $
        'fieldnam', 'theta range')
    cdf_save_var, the_var, value=df_theta_ranges, filename=cdfid
    cdf_save_setting, settings, filename=cdfid, varname=the_var

    the_var = 'width'
    settings = dictionary($
        'depend_0', time_var, $
        'unit', 'sec', $
        'display_type', 'scalar', $
        'short_name', 'dT', $
        'fieldnam', 'dtime')
    cdf_save_var, the_var, value=df_widths, filename=cdfid
    cdf_save_setting, settings, filename=cdfid, varname=the_var

    the_var = 'height'
    settings = dictionary($
        'depend_0', time_var, $
        'unit', 'deg', $
        'display_type', 'scalar', $
        'short_name', 'd'+tex2str('theta'), $
        'fieldnam', 'dtheta')
    cdf_save_var, the_var, value=df_heights, filename=cdfid
    cdf_save_setting, settings, filename=cdfid, varname=the_var

    the_var = 'scaled_height'
    settings = dictionary($
        'depend_0', time_var, $
        'scale_width', scale_width, $
        'unit', 'deg', $
        'display_type', 'scalar', $
        'short_name', 'd'+tex2str('theta'), $
        'fieldnam', 'scaled dtheta')
    cdf_save_var, the_var, value=df_scaled_heights, filename=cdfid
    cdf_save_setting, settings, filename=cdfid, varname=the_var

    the_var = 'obs_r_sm'
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'unit', 'Re', $
        'coord', 'SM', $
        'coord_labels', ['x','y','z'], $
        'short_name', 'R', $
        'fieldnam', 'R SM at observed time')
    cdf_save_var, the_var, value=df_obs_r_sms, filename=cdfid
    cdf_save_setting, settings, filename=cdfid, varname=the_var

    the_var = 'obs_mlt'
    settings = dictionary($
        'depend_0', time_var, $
        'unit', 'hr', $
        'display_type', 'scalar', $
        'short_name', 'MLT', $
        'fieldnam', 'MLT at observed time')
    cdf_save_var, the_var, value=df_obs_mlts, filename=cdfid
    cdf_save_setting, settings, filename=cdfid, varname=the_var

    the_var = 'obs_rxy'
    settings = dictionary($
        'depend_0', time_var, $
        'unit', 'Re', $
        'display_type', 'scalar', $
        'short_name', 'R!Dxy!N', $
        'fieldnam', 'Rxy at observed time')
    cdf_save_var, the_var, value=df_obs_rxys, filename=cdfid
    cdf_save_setting, settings, filename=cdfid, varname=the_var

    cdf_close, cdfid
end


function azim_df_read_ramp, time_range, probe=probe, project=project, $
    mlt_range=mlt_range, rxy_range=rxy_range


    retval = list()
    if n_elements(probe) eq 0 then return, retval
    if n_elements(project) eq 0 then project = azim_df_load_project()
    cdf_file = join_path([project.data_dir, 'azim_df_detected_ramp_'+probe+'.cdf'])
    if file_test(cdf_file) eq 0 then begin
        search_types = project.search_types
        txt_files = list()
        foreach search_type, search_types do begin
            index = where(search_type.probes eq probe, count)
            if count eq 0 then continue
            file_suffix = 'azim_df_detect_ramp_'+probe+'_'+search_type.name+'.txt'
            txt_files.add, join_path([project.data_dir,file_suffix])
        endforeach
        txt_files = txt_files.toarray()
        foreach txt_file, txt_files do begin
            if file_test(txt_file) eq 0 then begin
                stop
                azim_df_detect_ramp, project=project
            endif
        endforeach
        azim_df_read_ramp_gen_cdf, filename=cdf_file, txt_files=txt_files
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


    ramp_list = list()
    get_data, out_vars[0], obs_times
    if n_elements(obs_times) eq 1 then begin
        if obs_times[0] eq 0 then return, ramp_list
        the_ramp = dictionary($
            'probe', probe, $
            'obs_time', obs_times[0])
        foreach key, in_vars, var_id do begin
            the_data = get_var_data(out_vars[var_id])
            the_ramp[key] = (n_elements(the_data) eq 1)? the_data[0]: reform(the_data)
        endforeach
        ramp_list.add, the_ramp
    endif else begin
        data_list = list()
        data_dims = list()
        foreach var, out_vars do begin
            data = get_var_data(var)
            data_list.add, data
            data_dims.add, size(data,/n_dimension)
        endforeach
        foreach obs_time, obs_times, ii do begin
            the_ramp = dictionary($
                'probe', probe, $
                'obs_time', obs_time)
            foreach key, in_vars, var_id do begin
                the_ramp[key] = (data_dims[var_id] eq 1)? (data_list[var_id])[ii]: reform((data_list[var_id])[ii,*])
            endforeach
            ramp_list.add, the_ramp
        endforeach
    endelse

    return, ramp_list

end

;time_range = time_double(['2014-08-28/09:30','2014-08-28/11:00'])
time_range = time_double(['2009','2014'])
project = azim_df_load_project()
foreach probe, project.all_probes do ramp_list = azim_df_read_ramp(time_range, probe=probe)
end
