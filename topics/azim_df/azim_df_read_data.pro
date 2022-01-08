;+
; "Read" data from or calculated from the primitive data.
;-


; xxx_r_gsm, primitive.
pro azim_df_read_data_r_gsm, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    r_gsm_var = prefix+'r_gsm'
    if ~check_if_update(r_gsm_var, time_range) then return

    time_var = 'orbit_ut'
    cdf_load_var, r_gsm_var, range=time_range, time_var=time_var, filename=data_file
end

; xxx_b_sm, primitive.
pro azim_df_read_data_b_sm, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    b_sm_var = prefix+'b_sm'
    if ~check_if_update(b_sm_var, time_range) then return

    time_var = 'bfield_ut'
    cdf_load_var, b_sm_var, range=time_range, time_var=time_var, filename=data_file
end

; xxx_r_sm, calculated from xxx_r_gsm, save to file.
pro azim_df_read_data_r_sm, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    r_sm_var = prefix+'r_sm'
    if ~check_if_update(r_sm_var, time_range) then return

    r_gsm_var = prefix+'r_gsm'
    time_var = 'orbit_ut'
    if ~cdf_has_var(r_sm_var, filename=data_file) then begin
        times = cdf_read_var(time_var, filename=data_file)
        r_gsm = cdf_read_var(r_gsm_var, filename=data_file)
        r_sm = cotran(r_gsm, times, 'gsm2sm')
        cdf_save_var, r_sm_var, value=r_sm, filename=data_file
        setting = dictionary($
            'display_type', 'vector', $
            'short_name', 'R', $
            'unit', 'Re', $
            'coord', 'SM', $
            'coord_labels', constant('xyz'), $
            'depend_0', time_var)
        cdf_save_setting, setting, varname=r_sm_var, filename=data_file
    endif
    cdf_load_var, r_sm_var, range=time_range, time_var=time_var, filename=data_file
end

; xxx_mlt, calculated from xxx_r_gsm, save to file.
pro azim_df_read_data_mlt, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    mlt_var = prefix+'mlt'
    if ~check_if_update(mlt_var, time_range) then return

    r_gsm_var = prefix+'r_gsm'
    time_var = 'orbit_ut'
    if ~cdf_has_var(mlt_var, filename=data_file) then begin
        times = cdf_read_var(time_var, filename=data_file)
        r_gsm = cdf_read_var(r_gsm_var, filename=data_file)
        r_mag = cotran(r_gsm, times, 'gsm2mag')
        mlon = atan(r_mag[*,1],r_mag[*,0])*constant('deg')
        mlt = mlon2mlt(mlon, times)
        cdf_save_var, mlt_var, value=mlt, filename=data_file
        setting = dictionary($
            'display_type', 'scalar', $
            'short_name', 'MLT', $
            'unit', 'hr', $
            'range', [-12.,12], $
            'depend_0', time_var)
        cdf_save_setting, setting, varname=mlt_var, filename=data_file
    endif
    cdf_load_var, mlt_var, range=time_range, time_var=time_var, filename=data_file
end

; xxx_mlt_sm, calculated from xxx_r_sm, save to file.
pro azim_df_read_data_mlt_sm, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    mlt_sm_var = prefix+'mlt_sm'
    if ~check_if_update(mlt_sm_var, time_range) then return

    r_sm_var = prefix+'r_sm'
    time_var = 'orbit_ut'
    if ~cdf_has_var(mlt_sm_var, filename=data_file) then begin
        times = cdf_read_var(time_var, filename=data_file)
        r_sm = cdf_read_var(r_sm_var, filename=data_file)
        mlt_sm = azim_df_calc_pseudo_mlt(r_sm)
        cdf_save_var, mlt_sm_var, value=mlt_sm, filename=data_file
        setting = dictionary($
            'display_type', 'scalar', $
            'short_name', 'MLT!USM!N', $
            'unit', 'hr', $
            'range', [-12,12], $
            'depend_0', time_var)
        cdf_save_setting, setting, varname=mlt_sm_var, filename=data_file
    endif
    cdf_load_var, mlt_var, range=time_range, time_var=time_var, filename=data_file
end

; xxx_mlat, calculated from xxx_r_gsm, save to file.
pro azim_df_read_data_mlat, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    mlat_var = prefix+'mlat'
    if ~check_if_update(mlat_var, time_range) then return

    r_gsm_var = prefix+'r_gsm'
    time_var = 'orbit_ut'
    if ~cdf_has_var(mlat_var, filename=data_file) then begin
        times = cdf_read_var(time_var, filename=data_file)
        r_gsm = cdf_read_var(r_gsm_var, filename=data_file)
        r_mag = cotran(r_gsm, times, 'gsm2mag')
        mlat = asin(r_mag[*,2]/snorm(r_mag))*constant('deg')
        cdf_save_var, mlat_var, value=mlat, filename=data_file
        setting = dictionary($
            'display_type', 'scalar', $
            'short_name', 'MLat', $
            'unit', 'deg')
        cdf_save_setting, setting, varname=mlat_var, filename=data_file
    endif
    cdf_load_var, mlat_var, range=time_range, time_var=time_var, filename=data_file
end

; xxx_bmod_gsm, calculated from xxx_r_gsm, save to file.
pro azim_df_read_data_bmod_gsm, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    bmod_gsm_var = prefix+'bmod_gsm'
    if ~check_if_update(bmod_gsm_var, time_range) then return

    r_gsm_var = prefix+'r_gsm'
    time_var = 'orbit_ut'
    model = 't89'
    ndim = 3
    par = 2.
    if ~cdf_has_var(bmod_gsm_var, filename=data_file) then begin
        times = cdf_read_var(time_var, filename=data_file)
        r_gsm = cdf_read_var(r_gsm_var, filename=data_file)
        ntime = n_elements(times)
        bmod_gsm = fltarr(ntime, ndim)
        for ii=0, ntime-1 do begin
            tilt = geopack_recalc(times[ii])
            ; in-situ position
            rx = r_gsm[ii,0]
            ry = r_gsm[ii,1]
            rz = r_gsm[ii,2]
            ; in-situ B field.
            geopack_igrf_gsm, rx,ry,rz, bx,by,bz
            geopack_t89, par, rx,ry,rz, dbx,dby,dbz
            bmod_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endfor
        cdf_save_var, bmod_gsm_var, value=bmod_gsm, filename=data_file
        setting = dictionary($
            'display_type', 'vector', $
            'short_name', 'B!S!U'+strupcase(model)+'!N!S', $
            'unit', 'nT', $
            'model', strupcase(model), $
            'coord', 'GSM', $
            'coord_labels', constant('xyz'), $
            'par', par, $
            'depend_0', time_var)
        cdf_save_setting, setting, varname=bmod_gsm_var, filename=data_file
    endif
    cdf_load_var, bmod_gsm_var, range=time_range, time_var=time_var, filename=data_file
end

; xxx_bmod_tilt, calculated from xxx_bmod_gsm, do not save to file.
pro azim_df_read_data_bmod_tilt, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    bmod_tilt_var = prefix+'bmod_tilt'
    if ~check_if_update(bmod_tilt_var, time_range) then return

    bmod_gsm_var = prefix+'bmod_gsm'
    azim_df_read_data_bmod_gsm, time_range=time_range, filename=data_file, probe=probe
    get_data, bmod_gsm_var, times, bmod_gsm
    bmod_sm = cotran(bmod_gsm, times, 'gsm2sm')
    bmod_tilt = azim_df_calc_tilt(bmod_sm)
    store_data, bmod_tilt_var, times, bmod_tilt
    model = get_setting(bmod_gsm_var, 'model')
    add_setting, bmod_tilt_var, /smart, {$
        display_type: 'scalar', $
        short_name: tex2str('alpha')+'!D'+strupcase(model), $
        coord: 'SM', $
        model: strupcase(model), $
        unit: 'deg'}
end

; xxx_b_tilt, calculated from xxx_b_sm, do not save to file.
pro azim_df_read_data_b_tilt, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    b_tilt_var = prefix+'b_tilt'
    if ~check_if_update(b_tilt_var, time_range) then return

    b_sm_var = prefix+'b_sm'
    azim_df_read_data_b_sm, time_range=time_range, filename=data_file, probe=probe
    get_data, b_sm_var, times, b_sm
    b_tilt = azim_df_calc_tilt(b_sm)
    store_data, b_tilt_var, times, b_tilt
    add_setting, b_tilt_var, /smart, {$
        display_type: 'scalar', $
        short_name: tex2str('alpha')+'!Dmeas', $
        coord: 'SM', $
        unit: 'deg'}
end

; xxx_db_tilt, calculated from xxx_b_sm and xxx_r_gsm, do not save to file.
pro azim_df_read_data_db_tilt, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    db_tilt_var = prefix+'db_tilt'
    if ~check_if_update(db_tilt_var, time_range) then return

    b_tilt_var = prefix+'b_tilt'
    azim_df_read_data_b_tilt, time_range=time_range, filename=data_file, probe=probe
    bmod_tilt_var = prefix+'bmod_tilt'
    azim_df_read_data_bmod_tilt, time_range=time_range, filename=data_file, probe=probe

    get_data, b_tilt_var, times, b_tilt
    bmod_tilt = get_var_data(bmod_tilt_var, at=times)
    db_tilt = b_tilt-bmod_tilt
    store_data, db_tilt_var, times, db_tilt
    add_setting, db_tilt_var, /smart, {$
        display_type: 'scale', $
        short_name: tex2str('Delta')+tex2str('alpha'), $
        coord: 'SM', $
        unit: 'deg'}
end

; xxx_bmag, calculated from xxx_b_sm, do not save to file.
pro azim_df_read_data_bmag, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    bmag_var = prefix+'bmag'
    if ~check_if_update(bmag_var, time_range) then return

    b_sm_var = prefix+'b_sm'
    azim_df_read_data_b_sm, time_range=time_range, filename=data_file, probe=probe
    get_data, b_sm_var, times, b_sm
    bmag = snorm(b_sm)
    store_data, bmag_var, times, bmag
    add_setting, bmag_var, /smart, {$
        display_type: 'scalar', $
        short_name: '|B|', $
        unit: 'nT'}
end


; xxx_b_gsm, calculated from xxx_b_sm, do not save to file.
pro azim_df_read_data_b_gsm, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    b_gsm_var = prefix+'b_gsm'
    if ~check_if_update(b_gsm_var, time_range) then return

    b_sm_var = prefix+'b_sm'
    azim_df_read_data_b_sm, time_range=time_range, filename=data_file, probe=probe
    get_data, b_sm_var, times, b_sm
    b_gsm = cotran(b_sm, times, 'sm2gsm')
    store_data, b_gsm_var, times, b_gsm
    add_setting, b_gsm_var, /smart, {$
        display_type: 'vector', $
        short_name: 'B', $
        unit: 'nT', $
        coord: 'GSM', $
        coord_labels: constant('xyz') }
end


; xxx_theta, calculated from xxx_b_tilt, save to file.
pro azim_df_read_data_theta, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    theta_var = prefix+'theta'
    if ~check_if_update(theta_var, time_range) then return

    bmod_tilt_var = prefix+'bmod_tilt'
    b_tilt_var = prefix+'b_tilt'
    time_var = 'bfield_ut'
    if ~cdf_has_var(theta_var, filename=data_file) then begin
        times = cdf_read_var(time_var, filename=data_file)
        full_time_range = minmax(times)
        azim_df_read_data_b_tilt, time_range=full_time_range, filename=data_file, probe=probe
        azim_df_read_data_bmod_tilt, time_range=full_time_range, filename=data_file, probe=probe
        b_tilt = get_var_data(b_tilt_var, at=times)
        bmod_tilt = get_var_data(bmod_tilt_var, at=times)
        theta = b_tilt-bmod_tilt
        smooth_window = constant('secofhour')
        time_step = sdatarate(times)
        smooth_width = smooth_window/time_step
        theta -= smooth(theta, smooth_width, /edge_truncate, /nan)
        cdf_save_var, theta_var, value=theta, filename=data_file
        setting = dictionary($
            'display_type', 'scalar', $
            'short_name', tex2str('theta'), $
            'unit', 'deg', $
            'coord', 'SM', $
            'smooth_window', smooth_window, $
            'depend_0', time_var)
        cdf_save_setting, setting, varname=theta_var, filename=data_file
    endif
    cdf_load_var, theta_var, range=time_range, time_var=time_var, filename=data_file
end


; xxx_theta_median, calculated from xxx_theta, do not save to file.
pro azim_df_read_data_theta_median, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    theta_median_var = prefix+'theta_median'
    if ~check_if_update(theta_median_var, time_range) then return

    theta_var = prefix+'theta'
    azim_df_read_data_theta, time_range=time_range, filename=data_file, probe=probe
    get_data, theta_var, times, theta
    window_size = 120.  ; sec.
    time_step = times[1]-times[0]
    theta_median = sliding_boxcar(theta, window_size/time_step, type='median')
    store_data, theta_median_var, times, theta_mean
    add_setting, theta_median_var, /smart, {$
        display_type: 'scalar', $
        short_name: tex2str('theta')+' median', $
        window_size: window_size, $
        unit: 'deg'}

end


; xxx_dis, calculated from xxx_r_gsm, do not save to file.
pro azim_df_read_data_dis, time_range=time_range, filename=data_file, probe=probe
    prefix = probe+'_'
    dis_var = prefix+'dis'
    if ~check_if_update(dis_var, time_range) then return

    r_gsm_var = prefix+'r_gsm'
    azim_df_read_data_r_gsm, time_range=time_range, filename=data_file, probe=probe
    get_data, r_gsm_var, times, r_gsm
    store_data, dis_var, times, snorm(r_gsm)
    add_setting, dis_var, /smart, {$
        display_type: 'scalar', $
        short_name: 'R', $
        unit: 'Re'}
end


pro azim_df_read_data, var, probe=probe, time_range=time_range, project=project

;---Check inputs.
    if n_elements(var) eq 0 then begin
        errmsg = handle_error('No input var ...')
        return
    endif
    if n_elements(probe) eq 0 then begin
        errmsg = handle_error('No input probe ...')
        return
    endif
    if n_elements(time_range) ne 2 then begin
        errmsg = handle_error('Invalid time_range ...')
        return
    endif


;---Init search.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    search_settings = project.search_types

    foreach search, search_settings do begin
        search_time_range = search.time_range
        index = lazy_where(time_range, '[]', search_time_range, count=count)
        if count eq 0 then continue
        probes = search.probes
        index = where(probes eq probe, count)
        if count eq 0 then continue

        ; See if data file exists.
        data_file = join_path([project.data_dir,search.data_file_suffix])
        if file_test(data_file) eq 0 then azim_df_load_primitive_data, time_range=search_time_range, probe=probe, project=project, data_file=data_file

        ; Check for primitive data.
        prefix = probe+'_'
        vars = prefix+['r_gsm','b_sm']
        load = 0
        foreach tvar, vars do if ~cdf_has_var(tvar, filename=data_file) then load = 1
        if load eq 1 then azim_df_load_primitive_data, time_range=search_time_range, probe=probe, project=project, data_file=data_file

        ; Read the requested var.
        call_procedure, 'azim_df_read_data_'+var, time_range=time_range, filename=data_file, probe=probe
    endforeach

end

var = 'r_gsm'
probe = 'tha'
time_range = time_double(['2014-08-28','2014-08-29'])
azim_df_read_data, var, probe=probe, time_range=time_range, project=project
end
