function streamer_wave_load_data_themis, input_time_range, probe=probe, filename=data_file, version=version, projec_id=projec_id, event_id=event_id, _extra=ex

    update = 0

;---Check input.
    time_range = time_double(input_time_range)
    if n_elements(probe) eq 0 then message, 'No input probe ...'
    if n_elements(version) eq 0 then version = 'v01'
    if n_elements(projec_id) eq 0 then projec_id = 'streamer_wave'
    if n_elements(event_id) eq 0 then event_id = time_string(time_range[0],tformat='YYYY_MMDD_hh')
    mission = 'th'

    if n_elements(data_file) eq 0 then begin
        base = projec_id+'_'+event_id+'_'+mission+'_'+probe+'_data_'+version+'.cdf'
        data_file = join_path([googledir(),'works',projec_id,'data',base])
    endif

    if file_test(data_file) eq 0 then begin
        cdf_save_setting, 'probe', probe, filename=data_file
        cdf_save_setting, 'prefix', mission+probe+'_', filename=data_file
        cdf_save_setting, 'time_range', time_range, filename=data_file
        cdf_save_setting, 'mission', mission, filename=data_file
    endif
    sc_info = cdf_read_setting(filename=data_file)
    sc_info['data_file'] = data_file

    time_range = sc_info['time_range']
    probe = sc_info['probe']
    prefix = sc_info['prefix']
    mission_probe = mission+probe


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
    field_time_step = 3d
    field_time_var = lets_prep_time_var_in_memory(field_time_range, field_time_step, output=prefix+'field_time')
    b_gsm_var = lets_read_bfield(field_time_range, probe=mission_probe, $
        time_var=field_time_var, save_to=data_file, update=update)
    print, 'Loading '+b_gsm_var+' ...'
    
    b0_window = 15.*60
    sc_info['b0_window'] = b0_window
    bmod_var = prefix+'bmod_gsm_igrf_t89'
    b_vars = lets_decompose_bfield(b_var=b_gsm_var, $
        b0_window=b0_window, bmod_var=bmod_var, time_var=field_time_var, $
        save_to=data_file, update=update)
    foreach var, (b_vars.values()).toarray() do print, 'Loading '+var+' ...'
    



end

id = '2008_0328'
event_info = streamer_wave_load_data(id)
end