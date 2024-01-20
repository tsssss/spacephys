;+
; event_info needs to have time_range, sites.
;-

function event_study_read_asi_mlt_image, event_info, time_var=time_var, get_name=get_name, update=update

;---Settings.
    mlt_image_var = themis_asf_read_mlt_image(get_name=1)
    if keyword_set(get_name) then return, mlt_image_var
    if keyword_set(update) then del_data, mlt_image_var
    time_range = event_info['time_range']
    if ~check_if_update(mlt_image_var, time_range) then return, mlt_image_var


    data_file = event_info['data_file']
    if keyword_set(update) then cdf_del_var, mlt_image_var, filename=data_file
    if ~cdf_has_var(mlt_image_var, filename=data_file) then begin
        sites = event_info['sites']

        if ~event_info.haskey('min_elevs') then event_info['min_elevs'] = themis_asf_get_default_min_elevs(sites)

        if ~event_info.haskey('emission_height') then event_info['emission_height'] = 110d
        emission_height = event_info['emission_height']

        if ~event_info.haskey('merge_method') then event_info['merge_method'] = 'merge_elev'
        merge_method = event_info['merge_method']

        if ~event_info.haskey('calibration_method') then event_info['calibration_method'] = 'simple'
        calibration_method = event_info['calibration_method']

        mlt_image_var = themis_asf_read_mlt_image(time_range, sites=sites, min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method)

        times = get_var_time(mlt_image_var)
        if n_elements(time_var) eq 0 then time_var = 'asi_time'
        cdf_save_var, time_var, value=times, filename=data_file
        time_step = 3d
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', time_range )
        cdf_save_setting, settings, filename=data_file, varname=time_var

        data = get_var_data(mlt_image_var, limits=limits)
        cdf_save_var, mlt_image_var, value=data, filename=data_file

        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=mlt_image_var
    endif


    if check_if_update(mlt_image_var, time_range) then begin
        cdf_load_var, mlt_image_var, filename=data_file
    endif

    return, mlt_image_var


end