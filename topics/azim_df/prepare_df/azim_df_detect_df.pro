;+
; Detect DF from all ramps of a given probe.
;-

pro azim_df_detect_df, probe, log_file=log_file, data_file=data_file, txt_file=txt_file

;---Settings.
    project = azim_df_load_project()
    clear_window_size = 10.*60  ; sec.
    file_suffix = project.name+'_detected_df_'+probe
    if n_elements(log_file) eq 0 then log_file = join_path([project.data_dir,file_suffix+'.log'])
    if n_elements(data_file) eq 0 then data_file = join_path([project.data_dir,file_suffix+'.cdf'])
    if file_test(data_file) eq 1 then return
    if file_test(log_file) eq 0 then ftouch, log_file
    if n_elements(txt_file) eq 0 then txt_file = join_path([project.data_dir,file_suffix+'.txt'])
    if file_test(txt_file) eq 0 then ftouch, txt_file


;---Load all ramps.
    prefix = probe+'_'
    cdf_file = join_path([project.data_dir, 'azim_df_detected_ramp_'+probe+'.cdf'])
    if file_test(cdf_file) eq 0 then begin
        message, 'CDF file does not exist ...'
        tmp = azim_df_read_ramp(time_double(['2009','2014']), probe=probe, project=project)
    endif
    cdf_id = cdf_open(cdf_file)
    data_dict = dictionary()
    vars = ['obs_time','obs_r_sm','obs_rxy','obs_mlt','width','height','scaled_height','time_range','theta_range']
    foreach var, vars do data_dict[var] = cdf_read_var(var, filename=cdf_id)
    cdf_close, cdf_id


;---Select isolated ramps.
    obs_times = data_dict['obs_time']
    nramp = n_elements(obs_times)
    is_isolated_ramp = bytarr(nramp)+1
    dtimes = obs_times[1:nramp-1]-obs_times[0:nramp-2]
    index = where(dtimes lt clear_window_size, count)
    if count ne 0 then is_isolated_ramp[index] = 0
    lprmsg, '', log_file
    lprmsg, 'Found '+string(count,format='(I0)')+' ramps that are not isolated ...', log_file
    foreach ramp_id, index do begin
        ramp = dictionary('probe', probe)
        foreach key, data_dict.keys() do begin
            tmp = (data_dict[key])[ramp_id,*]
            ramp[key] = (n_elements(tmp) eq 1)? tmp[0]: tmp[*]
        endforeach
        azim_df_vertex_write, ramp, filename=log_file
    endforeach


    index = where(is_isolated_ramp eq 1, nramp)
    if nramp eq 0 then return
    foreach key, data_dict.keys() do data_dict[key] = (data_dict[key])[index,*]


;---Select ramps with proper shape.
    df_index = list()
    for ramp_id=0,nramp-1 do begin
        ramp = dictionary('probe', probe)
        foreach key, data_dict.keys() do begin
            tmp = (data_dict[key])[ramp_id,*]
            ramp[key] = (n_elements(tmp) eq 1)? tmp[0]: tmp[*]
        endforeach
        df = azim_df_filter_vertex(ramp, project=project, log_file=log_file)
        if n_elements(df) eq 0 then continue
        df_index.add, ramp_id
        azim_df_vertex_write, df, filename=txt_file
    endfor

    if df_index.length eq 0 then return
    index = df_index.toarray()
    foreach key, data_dict.keys() do data_dict[key] = (data_dict[key])[index,*]


;---Save data to cdf.
    df_obs_times = data_dict.obs_time
    df_time_ranges = data_dict.time_range
    df_theta_ranges = data_dict.theta_range
    df_widths = data_dict.width
    df_heights = data_dict.height
    df_scaled_heights = data_dict.scaled_height
    df_obs_r_sms = data_dict.obs_r_sm
    df_obs_mlts = data_dict.obs_mlt
    df_obs_rxys = data_dict.obs_rxy
    scale_width = project.scale_width

    cdfid = cdf_create(data_file)
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

project = azim_df_load_project()
foreach probe, project.all_probes do azim_df_detect_df, probe
end
