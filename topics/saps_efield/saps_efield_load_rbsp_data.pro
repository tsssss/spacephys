function saps_efield_load_rbsp_data, input_time_range, $
    probe=probe, filename=data_file, bad_e_trs=bad_e_trs, _extra=ex

    update = 0

;---Check input.
    time_range = time_double(input_time_range)
    if n_elements(probe) eq 0 then message, 'No input probe ...'

    if n_elements(data_file) eq 0 then begin
        base = time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_rbsp'+probe+'_data_v01.cdf'
        data_file = join_path([googledir(),'works','saps_efield','data',base])
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
                
                mlat_vars = lets_read_mlat_vars(orbit_var=fpt_var, save_to=data_file, time_var=orbit_time_var, var_info=var_info, update=update)
                foreach var, mlat_vars.values() do print, 'Loading '+var+' ...'
                
                ; The B at the footpoint.
                bf_var = lets_read_geopack_bfield(var_info=prefix+'bf_gsm'+suffix, $
                    orbit_var=fpt_var, time_var=orbit_time_var, update=update, $
                    internal_model=internal_model, external_model=external_model, save_to=data_file)
                print, 'Loading '+bf_var+' ...'
            endforeach
        endforeach
    endforeach


    foreach external_model, external_models do begin
        foreach internal_model, internal_models do begin
            ; The B at sc position.
            suffix = '_'+internal_model+'_'+external_model
            bmod_var = prefix+'bmod_gsm'+suffix
            if tnames(bmod_var) ne '' then del_data, bmod_var
            bmod_var = lets_read_geopack_bfield(var_info=bmod_var, $
                orbit_var=r_gsm_var, time_var=orbit_time_var, $
                internal_model=internal_model, external_model=external_model, save_to=data_file, update=update)
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
    options, b_gsm_var, 'mission_probe', 'rbsp'+probe
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
    options, edot0_fac_var, 'yrange', [-1,1]*15

    if n_elements(bad_e_trs) eq 0 then bad_e_trs = []
    foreach var, [e0_mgse_var,edot0_mgse_var,edot0_fac_var] do begin
        e_vec = get_var_data(var, times=times, limits=lim)
        foreach tr, bad_e_trs do begin
            index = where_pro(times,'[]', time_double(tr), count=count)
            if count ne 0 then begin
                e_vec[index,*] = !values.f_nan
            endif
        endforeach
        store_data, var, times, e_vec, limits=lim
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


;---E field spec.
    spec_var = lets_read_efield_spec(time_range, probe=mission_probe, update=update, save_to=data_file)
    b_spec_var = lets_read_this(func='rbsp_read_wave_spec_hz', $
        time_range, probe=mission_probe, $
        save_to=data_file)


    ; Cyclotron freqs.
    fc_vars = list()
    foreach species, ['e','o','he','p'] do begin
        fc_vars.add, lets_read_this(func='rbsp_read_gyro_freq', $
            time_range, probe=mission_probe, species=species, $
            save_to=data_file)
    endforeach

    var = prefix+'fce_half'
    fce = get_var_data(prefix+'fce', times=times)
    store_data, var, times, fce*0.5
    fc_vars.add, var
    var = prefix+'flh'
    fcp = get_var_data(prefix+'fcp', times=times)
    store_data, var, times, fcp*43
    fc_vars.add, var
    fc_vars = fc_vars.toarray()
    fc_colors = get_color(n_elements(fc_vars))
    foreach var, fc_vars, ii do options, var, 'colors', fc_colors[ii]

    spec_combo_var = spec_var+'_combo'
    store_data, spec_combo_var, data=[spec_var,fc_vars]
    options, spec_combo_var, 'yrange', get_setting(spec_var,'yrange')
    options, spec_combo_var, 'labels', ''


;---Particle data.
    dens_var = lets_read_this(func='rbsp_read_density', $
        time_range, probe=mission_probe, id='emfisis', suffix='', $
        save_to=data_file)
    
    foreach species, ['e','p','o'] do begin
        en_combo_var = lets_read_this(func='rbsp_read_en_spec_combo', $
            time_range, probe=mission_probe, species=species, errmsg=errmsg, $
            save_to=data_file)
        en_var = lets_read_this(func='rbsp_read_en_spec', $
            time_range, probe=mission_probe, species=species, errmsg=errmsg, $
            save_to=data_file)
        if species ne 'e' then begin
            energy_range = [50,5e4]
        endif else begin
            energy_range = [200,5e4]
        endelse
        pa_var = lets_read_this(func='rbsp_read_pa_spec', $
            time_range, probe=mission_probe, species=species, errmsg=errmsg, energy_range=energy_range, $
            save_to=data_file)
        
        t_var = lets_read_this(func='rbsp_read_temperature', $
            time_range, probe=mission_probe, species=species, errmsg=errmsg, energy_range=energy_range, $
            save_to=data_file)
    endforeach

    p_var = lets_read_this(func='rbsp_read_thermal_pressure', $
        time_range, probe=mission_probe, errmsg=errmsg, $
        save_to=data_file)
    


    return, sc_info
end


time_range = time_double(['2013-05-01','2013-05-02'])
probe = 'b'
time_range = time_double(['2013-11-11','2013-11-12'])
probe = 'a'
; 2014-09-12, -a.
time_range = time_double(['2015-03-17','2015-03-18'])
probe = 'b'
time_range = time_double(['2016-05-08','2016-05-09'])
probe = 'b'
time_range = time_double(['2015-01-04','2015-01-05'])
probe = 'b'
del_data, '*'
; 2015-06-23, -b and -a.
; 2015-08-15, -a and -b.
print, saps_efield_load_rbsp_data(time_range, probe=probe)
end