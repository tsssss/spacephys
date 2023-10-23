;+
; 
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





function themis_asf_read_mlt_image_2015_0416_calibrate_mlon_image, mlon_image_var, time_range

    ; customized background removal for mlon_image.
    get_data, mlon_image_var, times, mlon_images, limits=lim
    ntime = n_elements(mlon_images[*,0,0])
    npixel = product(lim.image_size)
    mlon_images = reform(mlon_images, [ntime,npixel])
    illuminated_indexs = where(total(mlon_images,1) ne 0, nilluminated_index)
    smooth_window = 2*60d
    time_step = 3d
    smooth_width = smooth_window/time_step
    foreach index, illuminated_indexs do begin
        data = mlon_images[*,index]
        time_index = where(data ne 0, count)
        if count eq 0 then continue
        if count lt smooth_width*3 then begin
            mlon_images[time_index,index] -= mean(data[time_index])
        endif else begin
    ;        mlt_images[time_index,index] -= smooth(data[time_index], smooth_width, edge_mirror=1)
            bg = calc_baseline_smooth(data[time_index], smooth_width)
            for ii=0,5 do bg = calc_baseline_smooth(bg, smooth_width)
            mlon_images[time_index,index] -= bg
        endelse
    endforeach

    mlon_images = reform(mlon_images, [ntime,lim.image_size])

    time_index = where_pro(times, '[)', time_range)
    mlon_images = mlon_images[time_index,*,*]
    times = times[time_index]
    store_data, mlon_image_var, times, mlon_images
    return, mlon_image_var

    ;tmp = mlon_images[where(mlon_images ne 0)]
    ;min_bg = (sort_uniq(tmp))[n_elements(tmp)*0.01]
    ;print, min_bg
    ;print, min(tmp)
    ;mlt_images += min_bg

end

function themis_asf_read_mlt_image_2015_0416, input_time_range, sites=sites, $
    errmsg=errmsg, get_name=get_name, update=update, $
    calibration_method=calibration_method, min_elev=min_elev, merge_method=merge_method, _extra=extra

    errmsg = ''
    retval = ''
    mlt_image_var = 'thg_asf_mlt_image'
    if keyword_set(get_name) then return, mlt_image_var
    if keyword_set(update) then del_data, mlt_image_var
    time_range = time_double(input_time_range)
    if ~check_if_update(mlt_image_var, time_range) then return, mlt_image_var

    ; get mlon_image.
    pad_time = 30*60d
    data_time_range = time_range+[-1,1]*pad_time
    mlon_image_var = themis_read_mlon_image(data_time_range, sites=sites, $
        min_elevs=min_elevs, resolutions=resolutions, $
        merge_method=merge_method, calibration_method=calibration_method)
    
    mlon_image_var = themis_asf_read_mlt_image_2015_0416_calibrate_mlon_image(mlon_image_var, time_range)

    ; convert to mlt_image.
    mlt_image_var = 'thg_asf_mlt_image'
    mlt_image_var = mlon_image_to_mlt_image(mlon_image_var, output=mlt_image_var)
    options, mlt_image_var, 'requested_time_range', time_range

    return, mlt_image_var

end


function themis_asf_read_mlt_image_rect_2015_0416, input_time_range, $
    errmsg=errmsg, get_name=get_name, $
    mlat_range=mlat_range, mlt_range=mlt_range, $
    output=mlt_image_var, calibration_method=calibration_method, _extra=ex

    if n_elements(mlt_image_var) eq 0 then mlt_image_var = 'thg_asf_mlt_image_rect'
    if keyword_set(get_name) then return, mlt_image_var
    time_range = time_double(input_time_range)
    if ~check_if_update(mlt_image_var, time_range) then return, mlt_image_var
    
    
;---Check input.
    if n_elements(mlat_range) eq 0 then mlat_range = [55.,85]

;---Load MLon image.
    pad_time = 30*60d
    data_time_range = time_range+[-1,1]*pad_time
    mlon_image_var = themis_asf_read_mlon_image_rect(data_time_range, calibration_method=calibration_method, _extra=ex)
    get_data, mlon_image_var, times, mlon_images, limits=lim
    index = where(finite(mlon_images,nan=1), count)
    if count ne 0 then mlon_images[index] = 0
    mlon_image_var = themis_asf_read_mlt_image_2015_0416_calibrate_mlon_image(mlon_image_var, time_range)
    mlon_images = get_var_data(mlon_image_var, times=times)
    
;---Prepare.
    illuminated_index = where(lim.illuminated_pixels eq 1)
    if n_elements(mlat_range) ne 2 then begin
        pixel_mlat = lim.pixel_mlat
        mlat_range = minmax(pixel_mlat[illuminated_index])
    endif
    if n_elements(mlt_range) ne 2 then begin
        pixel_mlon = lim.pixel_mlon
        mlon_range = minmax(pixel_mlon[illuminated_index])

        mlt_range = []
        foreach tmp, mlon_range do begin
            mlt_ranges = mlon2mlt(tmp, times)
            dmlt = mlt_ranges[1:-1]-mlt_ranges[0:-2]
            index = where(dmlt lt 0, count)
            for ii=0, count-1 do mlt_ranges[index[ii]+1:*] += 24
            mlt_range = [mlt_range, mlt_ranges]
        endforeach
        mlt_range = minmax(mlt_range)
    endif

    mlt_binsize = lim.dmlon/15.
    mlt_bins = make_bins(mlt_range, mlt_binsize)
    nmlt_bin = n_elements(mlt_bins)
    mlat_bins = lim.mlat_bins
    mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nmlat_bin)

;---Convert to MLT image.
    mltimg_size = [nmlt_bin,nmlat_bin]
    ntime = n_elements(times)
    fillval = !values.f_nan
    mlt_images = fltarr([ntime,mltimg_size])+fillval
    mlon_bins = lim.mlon_bins
    foreach time, times, time_id do begin
        the_mlt = mlon2mlt(mlon_bins, time)
        the_dmlt = the_mlt[1:-1]-the_mlt[0:-2]
        index = where(the_dmlt lt 0, count)
        for ii=0, count-1 do the_mlt[index[ii]+1:*] += 24
        index = where(the_mlt ge 12, count)
        if count ne 0 then the_mlt[index] -= 24 ; make it in [12,12]
        mlt_index = where_pro(the_mlt, '[]', minmax(mlt_bins), count=count)
        if count eq 0 then continue
        the_mlt = the_mlt[mlt_index]
        mlon_image = reform(mlon_images[time_id,mlt_index,mlat_index])
        mlt_image = sinterpol(mlon_image, the_mlt, mlt_bins)
        index = where_pro(mlt_bins,'][', minmax(the_mlt), count=count)
        if count ne 0 then mlt_image[index,*] = 0
        mlt_images[time_id,*,*] = mlt_image
    endforeach

    store_data, mlt_image_var, times, mlt_images, limits={$
        requested_time_range: time_range, $
        unit: '(#)', $
        image_size: mltimg_size, $
        mlt_range: mlt_range, $
        mlat_range: mlat_range, $
        mlt_bins: mlt_bins, $
        mlat_bins: mlat_bins[mlat_index] }
        
    return, mlt_image_var
    
end


function alfven_arc_load_ground_data_asi_2015_0416, event_info, filename=data_file, time_var=time_var, asi_setting=asi_setting

    event_info['asi_setting'] = asi_setting
    data_file = event_info['data_file']
    time_range = event_info['asi_time_range']
    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_ground_data_asi_ut(event_info, get_name=1)
    
    min_elevs = asi_setting['min_elevs']
    merge_method = asi_setting['merge_method']
    calibration_method = asi_setting['calibration_method']

    mlt_image_var = themis_asf_read_mlt_image_2015_0416(time_range, get_name=1)
    if ~cdf_has_var(mlt_image_var, filename=data_file) then begin
        sites = asi_setting['sites']
        mlt_image_var = themis_asf_read_mlt_image_2015_0416(time_range, sites=sites, min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method)
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

    mlt_image_rect_var = themis_asf_read_mlt_image_rect_2015_0416(time_range, get_name=1)
    if ~cdf_has_var(mlt_image_rect_var, filename=data_file) then begin
        sites = asi_setting['sites']
        mlt_range = asi_setting['mlt_range']
        mlat_range = asi_setting['mlat_range']
        mlt_image_rect_var = themis_asf_read_mlt_image_rect_2015_0416(time_range, sites=sites, $
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

    return, [mlt_image_var,mlt_image_rect_var]

end




function alfven_arc_load_ground_data_2015_0416, input_time_range, filename=data_file, asi_setting=asi_setting, version=version, _extra=ex

    time_range = time_double(input_time_range)
    if n_elements(version) eq 0 then version = 'v01'
    

    if n_elements(data_file) eq 0 then begin
        base = time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_ground_data_'+version+'.cdf'
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
    asi_vars = alfven_arc_load_ground_data_asi_2015_0416(event_info, filename=data_file, time_var=time_var, asi_setting=asi_setting)
    
    return, event_info

end



;get_data, mlt_image_var, times, mlt_images
;
;zlog = 1
;if zlog then begin
;    zrange = [1e1,1e4]
;    zzs = alog10(mlt_images)
;    zr = alog10(zrange)
;    ct = 49
;endif else begin
;    zrange = [-1,1]*1e3
;    zzs = mlt_images
;    zr = zrange
;    ct = 70
;endelse
;
;
;stop
;sgopen, 0, size=[1,1]*6
;angles = smkarthm(0,2*!dpi,40,'n')
;foreach time, times, time_id do begin
;    sgtv, bytscl(reform(zzs[time_id,*,*]),top=254, min=zr[0], max=zr[1]), position=[0,0,1,1], ct=ct
;    plot, [-1,1],[-1,1], nodata=1, noerase=1, position=[0,0,1,1], xstyle=5, ystyle=5
;    foreach rr, smkarthm(0,1,0.125,'dx') do plots, rr*cos(angles),rr*sin(angles), color=sgcolor('silver'), linestyle=1
;    foreach rr, smkarthm(0,1,0.25,'dx') do plots, rr*cos(angles),rr*sin(angles), color=sgcolor('silver'), linestyle=0
;    foreach tt, smkarthm(0,2*!dpi,25,'n') do plots, [0,1]*cos(tt),[0,1]*sin(tt), color=sgcolor('silver'), linestyle=1
;    foreach tt, smkarthm(0,2*!dpi,13,'n') do plots, [0,1]*cos(tt),[0,1]*sin(tt), color=sgcolor('silver'), linestyle=0
;    xyouts, 10,10, device=1, time_string(time)
;    if time eq time_double('2015-04-16/08:03') then stop
;    if time eq time_double('2015-04-16/08:06') then stop
;    if time eq time_double('2015-04-16/08:09') then stop
;endforeach




;    mlt_image_var = 'thg_asf_mlt_image'
;    if keyword_set(get_name) then return, mlt_image_var
;    time_range = time_double(input_time_range)
;    if ~check_if_update(mlt_image_var, time_range) then return, mlt_image_var
;
;    mlon_image_var = themis_asf_read_mlon_image(input_time_range, sites=sites, $
;        min_elev=min_elev, merge_method=merge_method, errmsg=errmsg, calibration_method=calibration_method)
;    if errmsg ne '' then return, retval
;
;;---Rotate from mlon to mlt.
;    mlt_image_var = mlon_image_to_mlt_image(mlon_image_var, output=mlt_image_var)
;    options, mlt_image_var, 'requested_time_range', time_range
;    return, mlt_image_var
;
;end