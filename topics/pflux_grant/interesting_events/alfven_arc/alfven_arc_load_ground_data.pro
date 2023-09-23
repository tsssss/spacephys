;+
;-

function alfven_arc_load_ground_data_weygand_ut, event_info, time_var=var, get_name=get_name

    if n_elements(var) eq 0 then var = 'weygand_ut'
    if keyword_set(get_name) then return, var

    data_file = event_info['data_file']
    if ~cdf_has_var(var, filename=data_file) then begin
        time_range = event_info['weygand_time_range']
        time_step = event_info['weygand_time_step']
        times = make_bins(time_range, time_step)
        cdf_save_var, var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', time_range )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    return, cdf_read_var(var, filename=data_file)

end

function alfven_arc_load_ground_data_weygand, event_info, time_var=time_var

    data_file = event_info['data_file']
    time_range = event_info['weygand_time_range']
    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_ground_data_weygand_ut(get_name=1)

    types = ['hor','ver']
    j_vars = 'thg_j_'+types
    foreach type, types do begin
        var = 'thg_j_'+type
        if ~cdf_has_var(var, filename=data_file) then begin
            var = themis_read_weygand_j(time_range, id='j_'+type)

            data = get_var_data(var, limits=limits)
            cdf_save_var, var, value=data, filename=data_file
            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=var
        endif

        if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    endforeach


    mlat_range = [50d,90]   ; deg
    var = 'thg_j_ver_mlt_image'
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

    if check_if_update(var, time_range) then begin
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
    endif

    j_vars = [j_vars,var]
    return, j_vars
    
end


function alfven_arc_load_ground_data_asi_ut, event_info, time_var=var, get_name=get_name

    if n_elements(var) eq 0 then var = 'asi_ut'
    if keyword_set(get_name) then return, var

    data_file = event_info['data_file']
    if ~cdf_has_var(var, filename=data_file) then begin
        time_range = event_info['asi_time_range']
        time_step = event_info['asi_time_step']
        times = make_bins(time_range+[0,-1]*time_step, time_step)
        cdf_save_var, var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', time_range )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    return, cdf_read_var(var, filename=data_file)

end


function alfven_arc_load_ground_data_asi, event_info, filename=data_file, time_var=time_var, asi_setting=asi_setting

    event_info['asi_setting'] = asi_setting
    data_file = event_info['data_file']
    time_range = event_info['asi_time_range']
    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_ground_data_asi_ut(event_info, get_name=1)
    
    min_elevs = asi_setting['min_elevs']
    merge_method = asi_setting['merge_method']
    calibration_method = asi_setting['calibration_method']

    mlt_image_var = themis_asf_read_mlt_image(time_range, get_name=1)
    if ~cdf_has_var(mlt_image_var, filename=data_file) then begin
        sites = asi_setting['sites']
        mlt_image_var = themis_asf_read_mlt_image(time_range, sites=sites, min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method)
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

    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
    if ~cdf_has_var(mlt_image_rect_var, filename=data_file) then begin
        sites = asi_setting['sites']
        mlt_range = asi_setting['mlt_range']
        mlat_range = asi_setting['mlat_range']
        mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, sites=sites, $
            min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method, $
            mlt_range=mlt_range, mlat_range=mlat_range)
        data = get_var_data(mlt_image_rect_var, limits=limits)
        cdf_save_var, mlt_image_rect_var, value=data, filename=data_file

        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=mlt_image_rect_var
    endif
    if check_if_update(mlt_image_rect_var, time_range) then begin
        cdf_load_var, mlt_image_rect_var, filename=data_file
    endif

;;---Ewogram.
;;---Keogram.
;    fmlt_vars = prefix+'fmlt_'+['dipole','igrf']
;;    mlt_range = []
;;    foreach var, fmlt_vars do begin
;;        mlt_range = [mlt_range,minmax(get_var_data(var,in=time_range))]
;;    endforeach
;;    mlt_range = minmax(mlt_range)
;    mlt_range = []
;    model_setting = event_info['model_setting']
;    internal_model = model_setting['internal_model']
;    mlt_var = prefix+'fmlt_'+internal_model+'_north'
;    models = model_setting['models']
;    external_model = model_setting['external_model']
;    model_index = where(models eq external_model)
;    mlt_range = reform((get_var_data(mlt_var,at=time_range))[*,model_index])
;    mlt_range = round(mlt_range*10)*0.1
;    mlat_range = [60,67]
;    keo_var = themis_read_mlt_image_rect_keo(mlt_image_var=mlt_image_rect_var, $
;        mlt_range=mlt_range, mlat_range=mlat_range)
;    zlim, keo_var, 1,8e3, 0

    return, [mlt_image_var,mlt_image_rect_var]

end


function alfven_arc_load_ground_data, input_time_range, filename=data_file, asi_setting=asi_setting, _extra=ex

    time_range = time_double(input_time_range)

    if n_elements(data_file) eq 0 then begin
        base = time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_ground_data_v01.cdf'
        data_file = join_path([googledir(),'works','pflux_grant','alfven_arc','data',base])
    endif
    

    if file_test(data_file) eq 0 then begin
        cdf_save_setting, 'time_range', time_range, filename=data_file
    endif
    event_info = cdf_read_setting(filename=data_file)
    event_info['data_file'] = data_file

    time_range = event_info['time_range']

;---SECS.
    weygand_time_var = 'weygand_ut'
    weygand_time_step = 10d
    event_info['weygand_time_var'] = weygand_time_var
    event_info['weygand_time_step'] = weygand_time_step
    event_info['weygand_time_range'] = time_range
    times = alfven_arc_load_ground_data_weygand_ut(event_info, time_var=weygand_time_var)
    j_vars = alfven_arc_load_ground_data_weygand(event_info, time_var=weygand_time_var)

;---ASI.
    asi_time_var = 'asi_ut'
    asi_time_step = 3d
    event_info['asi_time_var'] = asi_time_var
    event_info['asi_time_step'] = asi_time_step
    event_info['asi_time_range'] = time_range
    times = alfven_arc_load_ground_data_asi_ut(event_info, time_var=asi_time_var)

    if n_elements(asi_setting) eq 0 then message, 'No input asi_setting ...'
;    asi_setting = dictionary($
;        'sites', ['fykn','mcgr'], $
;        'min_elevs', float([5,10]), $
;        'best_site', 'fykn', $
;        'mlt_range', [-6d,0], $
;        'mlat_range', [55d,70] )
    asi_vars = alfven_arc_load_ground_data_asi(event_info, filename=data_file, time_var=time_var, asi_setting=asi_setting)
    
    return, event_info

end