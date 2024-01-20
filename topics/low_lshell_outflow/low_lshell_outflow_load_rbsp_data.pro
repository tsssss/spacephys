
function low_lshell_outflow_load_rbsp_data, input_time_range, $
    probe=probe, filename=data_file, bad_e_time_ranges=bad_e_time_ranges, _extra=ex

    update = 0

;---Check input.
    time_range = time_double(input_time_range)
    if n_elements(probe) eq 0 then message, 'No input probe ...'

    if n_elements(data_file) eq 0 then begin
        base = time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_rbsp'+probe+'_data_v01.cdf'
        data_file = join_path([googledir(),'works','low_lshell_outflow','data',base])
    endif

    if file_test(data_file) eq 0 then begin
        cdf_save_setting, 'probe', probe, filename=data_file
        cdf_save_setting, 'prefix', 'rbsp'+probe+'_', filename=data_file
        cdf_save_setting, 'time_range', time_range, filename=data_file
        cdf_save_setting, 'mission', 'rbsp', filename=data_file
    endif
    sc_info = cdf_read_setting(filename=data_file)
    sc_info['data_file'] = data_file

    time_range = sc_info['time_range']
    probe = sc_info['probe']
    prefix = sc_info['prefix']
    mission_probe = 'rbsp'+probe


;---Orbit related vars.
    time_step = 60.
    pad_time = 30.*60
    orbit_time_range = time_range+[-1,1]*pad_time
    orbit_time_var = lets_prep_time_var_in_memory(orbit_time_range, time_step, output=prefix+'orbit_time')
    sc_info['orbit_time_var'] = orbit_time_var
    r_gsm_var = lets_read_orbit(orbit_time_range, probe=mission_probe, time_var=orbit_time_var, save_to=data_file)
    print, 'Loading '+r_gsm_var+' ...'
    
    ; The fmlat/fmlon/mlt of the orbit.
    var_info = dictionary($
        'mlat', prefix+'mlat', $
        'mlon', prefix+'mlon', $
        'mlt', prefix+'mlt', $
        'lshell', prefix+'lshell', $
        'dis', prefix+'dis' )
    mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var, save_to=data_file, time_var=orbit_time_var, var_info=var_info, update=update)
    foreach var, mlat_vars.values() do print, 'Loading '+var+' ...'
    
    
    external_models = ['t89','t96','t01','t04s']
    internal_models = ['dipole','igrf']
    hemispheres = ['north','south']
    sc_info['external_models'] = external_models
    sc_info['internal_models'] = internal_models
    foreach external_model, external_models do begin
        foreach internal_model, internal_models do begin
            foreach hemisphere, hemispheres do begin
                suffix = '_'+internal_model+'_'+external_model+'_'+hemisphere

                fpt_var = lets_trace_to_ionosphere(orbit_var=r_gsm_var, coord='gsm', time_var=orbit_time_var, $
                    internal_model=internal_model, external_model=external_model, hemisphere=hemisphere, save_to=data_file)
                print, 'Loading '+fpt_var+' ...'
                
                ; The fmlat/fmlon/mlt of the footpoint.
                var_info = dictionary($
                    'mlat', prefix+'fmlat'+suffix, $
                    'mlon', prefix+'fmlon'+suffix, $
                    'mlt', prefix+'fmlt'+suffix )
                                    
                mlat_vars = lets_read_mlat_vars(orbit_var=fpt_var, save_to=data_file, time_var=orbit_time_var, var_info=var_info)
                foreach var, mlat_vars.values() do print, 'Loading '+var+' ...'
                
                ; The B at the footpoint.
                bf_var = lets_read_geopack_bfield(var_info=prefix+'bf_gsm'+suffix, $
                    orbit_var=fpt_var, time_var=orbit_time_var, update=update, $
                    internal_model=internal_model, external_model=external_model, save_to=data_file)
                print, 'Loading '+bf_var+' ...'
            endforeach
            
            
            ; The B at sc position.
            suffix = '_'+internal_model+'_'+external_model
            bmod_var = lets_read_geopack_bfield(var_info=prefix+'bmod_gsm'+suffix, $
                orbit_var=r_gsm_var, time_var=orbit_time_var, $
                internal_model=internal_model, external_model=external_model, save_to=data_file)
            print, 'Loading '+bmod_var+' ...'
        endforeach
    endforeach



;---Bfield related vars.
    pad_time = 30.*60
    field_time_range = time_range+[-1,1]*pad_time
    field_time_step = 1d/32
    field_time_var = lets_prep_time_var_in_memory(field_time_range, field_time_step, output=prefix+'field_time')
    b_gsm_var = lets_read_bfield(field_time_range, probe=mission_probe, $
        time_var=field_time_var, save_to=data_file, resolution='hires', update=update, remove_spin_tone=11d)
    print, 'Loading '+b_gsm_var+' ...'
    
    b0_window = 15.*60
    sc_info['b0_window'] = b0_window
    bmod_var = prefix+'bmod_gsm_igrf_t89'
    b_vars = lets_decompose_bfield(b_var=b_gsm_var, $
        b0_window=b0_window, bmod_var=bmod_var, time_var=field_time_var, $
        save_to=data_file, update=update)
    foreach var, (b_vars.values()).toarray() do print, 'Loading '+var+' ...'
    


;---Efield related vars.
    e0_mgse_var = lets_read_efield(field_time_range, probe=mission_probe, $
        time_var=field_time_var, save_to=data_file)
    print, 'Loading '+e0_mgse_var+' ...'

    ; Remove bad E.
    if probe eq 'a' then begin
        bad_e_trs = list()
        bad_e_trs.add, ['2015-03-17/03:11:30','2015-03-17/03:27:00']
        bad_e_trs.add, ['2015-03-17/08:16:42','2015-03-17/08:16:44']
        bad_e_trs.add, ['2015-03-17/08:38:35','2015-03-17/08:38:38']
        bad_e_trs.add, ['2015-03-17/08:39:07','2015-03-17/08:39:11']
        bad_e_trs.add, ['2015-03-17/08:51:05','2015-03-17/08:51:08']
        bad_e_trs.add, ['2015-03-17/08:51:16','2015-03-17/08:51:37']
        bad_e_trs.add, ['2015-03-17/09:06:00','2015-03-17/09:08:34']
        bad_e_trs.add, ['2015-03-17/09:13:09','2015-03-17/09:19:26']
        bad_e_trs.add, ['2015-03-17/09:20:10','2015-03-17/09:40:00']
        bad_e_trs.add, ['2015-03-17/11:11:48','2015-03-17/11:11:57']
        bad_e_trs.add, ['2015-03-17/11:13:11','2015-03-17/11:13:14']
        bad_e_trs.add, ['2015-03-17/11:34:20','2015-03-17/11:34:40']
        bad_e_trs.add, ['2015-03-17/12:07:59','2015-03-17/12:12:54']
        bad_e_trs.add, ['2015-03-17/12:54:35','2015-03-17/12:55:27']
        bad_e_trs.add, ['2015-03-17/12:56:07','2015-03-17/12:56:39']
        bad_e_trs.add, ['2015-03-17/13:03:47','2015-03-17/13:09:01']
        bad_e_trs.add, ['2015-03-17/13:09:32','2015-03-17/13:10:26']
        bad_e_trs.add, ['2015-03-17/13:12:38','2015-03-17/13:13:22']
        bad_e_trs.add, ['2015-03-17/13:17:14','2015-03-17/13:17:24']
        bad_e_trs.add, ['2015-03-17/13:19:02','2015-03-17/13:22:44']
        bad_e_trs.add, ['2015-03-17/13:26:27','2015-03-17/13:27:29']
        bad_e_trs.add, ['2015-03-17/16:40:15','2015-03-17/16:41:17']
        bad_e_trs.add, ['2015-03-17/16:42:09','2015-03-17/16:53:08']
        bad_e_trs.add, ['2015-03-17/17:36:05','2015-03-17/17:37:50']
        bad_e_trs.add, ['2015-03-17/17:42:55','2015-03-17/17:44:25']
        bad_e_trs.add, ['2015-03-17/17:45:22','2015-03-17/17:45:35']
        bad_e_trs.add, ['2015-03-17/17:46:00','2015-03-17/17:51:19']
        bad_e_trs.add, ['2015-03-17/18:00:54','2015-03-17/18:03:44']
        bad_e_trs.add, ['2015-03-17/18:09:10','2015-03-17/18:13:30']
        bad_e_trs.add, ['2015-03-17/18:21:30','2015-03-17/18:25:09']
        bad_e_trs.add, ['2015-03-17/18:26:44','2015-03-17/18:50:10']
        bad_e_trs.add, ['2015-03-17/19:57:00','2015-03-17/21:12:10']
        bad_e_trs.add, ['2015-03-17/21:56:40','2015-03-17/21:58:30']
        bad_e_trs.add, ['2015-03-17/22:00:39','2015-03-17/22:05:59']
        bad_e_trs.add, ['2015-03-17/22:12:29','2015-03-17/22:26:30']
        bad_e_trs.add, ['2015-03-18/02:36:02','2015-03-18/02:40:23']
        bad_e_trs.add, ['2015-03-18/05:05:00','2015-03-18/05:40:30']
        bad_e_trs.add, ['2015-03-18/06:06:39','2015-03-18/06:08:29']
        bad_e_trs.add, ['2015-03-18/06:08:59','2015-03-18/06:10:19']
        bad_e_trs.add, ['2015-03-18/06:32:59','2015-03-18/06:58:30']
        bad_e_trs.add, ['2015-03-18/11:20:15','2015-03-18/11:26:10']
        bad_e_trs.add, ['2015-03-18/11:45:15','2015-03-18/12:30:00']
    endif else begin
        bad_e_trs = list()
        bad_e_trs.add, ['2015-03-16/23:59:59','2015-03-17/00:03:04']
        bad_e_trs.add, ['2015-03-17/08:53:10','2015-03-17/08:55:30']
        bad_e_trs.add, ['2015-03-17/09:00:24','2015-03-17/09:01:04']
        bad_e_trs.add, ['2015-03-17/09:04:10','2015-03-17/09:07:29']
        bad_e_trs.add, ['2015-03-17/09:09:29','2015-03-17/09:12:27']
        bad_e_trs.add, ['2015-03-17/09:13:31','2015-03-17/09:14:43']
        bad_e_trs.add, ['2015-03-17/09:16:05','2015-03-17/09:19:01']
        bad_e_trs.add, ['2015-03-17/09:33:03','2015-03-17/09:33:39']
        bad_e_trs.add, ['2015-03-17/09:34:15','2015-03-17/09:35:03']
        bad_e_trs.add, ['2015-03-17/09:36:39','2015-03-17/09:37:50']
        bad_e_trs.add, ['2015-03-17/13:08:17','2015-03-17/13:08:21']
        bad_e_trs.add, ['2015-03-17/13:09:30','2015-03-17/13:09:59']
        bad_e_trs.add, ['2015-03-17/13:10:46','2015-03-17/13:11:28']
        bad_e_trs.add, ['2015-03-17/13:15:04','2015-03-17/13:15:27']
        bad_e_trs.add, ['2015-03-17/13:30:58','2015-03-17/13:32:07']
        bad_e_trs.add, ['2015-03-17/13:33:40','2015-03-17/13:33:58']
        bad_e_trs.add, ['2015-03-17/13:34:21','2015-03-17/13:34:43']
        bad_e_trs.add, ['2015-03-17/16:04:29','2015-03-17/16:04:35']
        bad_e_trs.add, ['2015-03-17/17:38:47','2015-03-17/17:38:49']
        bad_e_trs.add, ['2015-03-17/17:46:41','2015-03-17/17:48:32']
        bad_e_trs.add, ['2015-03-17/17:50:35','2015-03-17/17:56:45']
        bad_e_trs.add, ['2015-03-17/18:01:25','2015-03-17/18:02:35']
        bad_e_trs.add, ['2015-03-17/18:04:40','2015-03-17/18:06:35']
        bad_e_trs.add, ['2015-03-17/18:07:50','2015-03-17/18:09:25']
        bad_e_trs.add, ['2015-03-17/18:40:10','2015-03-17/18:41:45']
        bad_e_trs.add, ['2015-03-17/18:45:05','2015-03-17/18:58:05']
        bad_e_trs.add, ['2015-03-17/23:08:35','2015-03-17/23:20:35']
        bad_e_trs.add, ['2015-03-17/23:37:20','2015-03-17/23:41:20']
        bad_e_trs.add, ['2015-03-17/23:43:30','2015-03-17/23:49:05']
        bad_e_trs.add, ['2015-03-17/23:52:55','2015-03-17/23:59:05']
        bad_e_trs.add, ['2015-03-18/00:14:00','2015-03-18/00:37:30']
        bad_e_trs.add, ['2015-03-18/02:06:55','2015-03-18/02:35:29']
        bad_e_trs.add, ['2015-03-18/02:56:04','2015-03-18/03:02:59']
        bad_e_trs.add, ['2015-03-18/03:07:34','2015-03-18/03:09:04']
        bad_e_trs.add, ['2015-03-18/03:09:54','2015-03-18/03:11:09']
        bad_e_trs.add, ['2015-03-18/03:15:04','2015-03-18/03:19:34']
        bad_e_trs.add, ['2015-03-18/03:21:18','2015-03-18/03:22:08']
        bad_e_trs.add, ['2015-03-18/03:22:59','2015-03-18/03:30:39']
        bad_e_trs.add, ['2015-03-18/08:28:59','2015-03-18/08:43:44']
        bad_e_trs.add, ['2015-03-18/08:48:09','2015-03-18/08:49:49']
        bad_e_trs.add, ['2015-03-18/11:57:40','2015-03-18/12:00:57']
        bad_e_trs.add, ['2015-03-18/12:02:03','2015-03-18/12:08:27']
        bad_e_trs.add, ['2015-03-18/12:09:57','2015-03-18/12:10:45']
        bad_e_trs.add, ['2015-03-18/12:12:07','2015-03-18/12:12:37']
        bad_e_trs.add, ['2015-03-18/12:13:09','2015-03-18/12:13:11']
    endelse

    e_mgse = get_var_data(e0_mgse_var, times=times, limits=lim)
    foreach tr, bad_e_trs do begin
        index = where_pro(times,'[]', time_double(tr), count=count) 
        if count eq 0 then continue
        e_mgse[index,*] = !values.f_nan
    endforeach
    store_data, e0_mgse_var, times, e_mgse, limits=lim


    edot0_mgse_var = lets_calc_edotb0(e_var=e0_mgse_var, b_var=b_vars['b0'], $
        var_info=prefix+'edot0_rbsp_mgse', save_to=data_file, time_var=field_time_var, update=update)
    print, 'Loading '+edot0_mgse_var+' ...'
    
    edot0_mgse = get_var_data(edot0_mgse_var, times=times, settings=e_settings)
    b_angle = e_settings['b_angle']
    min_b_angle = 15d
    index = where(abs(b_angle) le min_b_angle, count)
    if count ne 0 then begin
        edot0_mgse[index,0] = !values.f_nan
        store_data, edot0_mgse_var, times, edot0_mgse
    endif

    q_gsm2fac_var = lets_define_fac(b_var=b_vars['b0'], r_var=r_gsm_var, save_to=data_file, time_var=orbit_time_var)
    edot0_fac_var = prefix+'edot0_fac'
    b1_fac_var = prefix+'b1_fac'
    foreach var, [edot0_mgse_var,b_vars['b1']] do begin
        in_coord = strlowcase(get_setting(var,'coord'))
        out_var = streplace(var,in_coord,'fac')
        var_fac = lets_cotran([in_coord,'fac'], input=var, output=out_var, q_var=q_gsm2fac_var, save_to=data_file, update=update)
    endforeach

    
    scale_info = dictionary($
        's0', 1d/4, $
        's1', 30d*60, $
        'dj', 1d/4 )
    pf_var = prefix+'pf_fac'
    pf_var = lets_calc_pflux(e_var=edot0_fac_var, b_var=b1_fac_var, pf_var=pf_var, $
        scale_info=scale_info, save_to=data_file, update=update)
    
    
    cmap_var = prefix+'cmap'
    if check_if_update(cmap_var, time_range) then begin
        bf_gsm = get_var_data(prefix+'bf_gsm_dipole_t89_north', times=times)
        cmap = snorm(bf_gsm)/snorm(get_var_data(prefix+'b_gsm', at=times))
        store_data, cmap_var, times, cmap, limits={requested_time_range:time_range}
    endif
    
    pf_map_var = prefix+'pf_fac_map'
    if check_if_update(pf_map_var, time_range) then begin
        pf_fac = get_var_data(pf_var, times=times, limits=lim)
        cmap = get_var_data(cmap_var, at=times)
        ndim = 3
        for ii=0,ndim-1 do pf_fac[*,ii] *= cmap
        store_data, pf_map_var, times, pf_fac, limits=lim
    endif

    return, sc_info
end


id = '2015_0317'
event_info = low_lshell_outflow_load_data(id)
end