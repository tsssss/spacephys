
function saps_efield_load_dmsp_data, input_time_range, $
    probe=probe, filename=data_file, _extra=ex

    update = 0

;---Check input.
    time_range = time_double(input_time_range)
    if n_elements(probe) eq 0 then message, 'No input probe ...'

    if n_elements(data_file) eq 0 then begin
        base = time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_dmsp'+probe+'_data_v01.cdf'
        data_file = join_path([googledir(),'works','saps_efield','data',base])
    endif

    if file_test(data_file) eq 0 then begin
        cdf_save_setting, 'probe', probe, filename=data_file
        cdf_save_setting, 'prefix', 'dmsp'+probe+'_', filename=data_file
        cdf_save_setting, 'time_range', time_range, filename=data_file
        cdf_save_setting, 'mission', 'dmsp', filename=data_file
    endif
    sc_info = cdf_read_setting(filename=data_file)
    sc_info['data_file'] = data_file

    time_range = sc_info['time_range']
    probe = sc_info['probe']
    prefix = sc_info['prefix']
    mission_probe = 'dmsp'+probe

;---SSUSI.
    ssusi_id = 'energy'
    dmsp_mlt_image_var = lets_read_this(time_range, probe=mission_probe, $
        func='dmsp_read_mlt_image', save_to=data_file, update=update, id=ssusi_id)

;---Orbit related vars.
    mlat_vars = lets_read_this(time_range, probe=mission_probe, $
        func='dmsp_read_mlat_vars', save_to=data_file, update=update)
    
    ; The fmlat/fmlon/mlt of the orbit.
    foreach var, mlat_vars do print, 'Loading '+var+' ...'


;---B0 and dB in dmsp_xyz.
    b0_xyz_var = lets_read_this(time_range, probe=mission_probe, $
        func='dmsp_read_bfield_madrigal', read_b0=1, $
        save_to=data_file, update=update)
    db_xyz_var = lets_read_this(time_range, probe=mission_probe, $
        func='dmsp_read_bfield_madrigal', read_b0=0, $
        save_to=data_file, update=update)
    print, 'Loading '+b0_xyz_var+' ...'
    print, 'Loading '+db_xyz_var+' ...'

    b_vec = get_var_data(db_xyz_var, times=times)
    index = where(finite(snorm(b_vec)), count)
    b_vec = sinterpol(b_vec[index,*], times[index], times)
    store_data, db_xyz_var, times, b_vec

    b_vec = get_var_data(b0_xyz_var, times=times)
    bmag_var = prefix+'bmag'
    store_data, bmag_var, times, snorm(b_vec)
    add_setting, bmag_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', '|B|', $
        'unit', 'nT')
    options, bmag_var, 'ytitle', '|B| (nT)'

;---Velocity.
    v_xyz_var = lets_read_this(time_range, probe=mission_probe, $
        func='dmsp_read_ion_vel_madrigal', save_to=data_file, update=update)
    print, 'Loading '+v_xyz_var+' ...'

;---en spec.
    foreach species, ['e','p'] do begin
        var = lets_read_this(time_range, probe=mission_probe, $
            func='dmsp_read_en_spec', save_to=data_file, update=update)
        print, 'Loading '+var+' ...'
    endforeach
    
    return, sc_info

end