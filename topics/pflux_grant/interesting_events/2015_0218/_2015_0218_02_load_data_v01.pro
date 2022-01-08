function _2015_0218_02_load_orbit_ut, event_info, varname=var, $
    time_step=time_step, pad_time=pad_time

    time_range = event_info['time_range']
    data_file = event_info['data_file']

    if n_elements(var) eq 0 then var = 'orbit_ut'
    if ~cdf_has_var(var, filename=data_file) then begin
        time_step = event_info['orbit_time_step']
        data_time_range = event_info['orbit_time_range']
        times = make_bins(data_time_range, time_step)
        cdf_save_var, var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', data_time_range )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    return, cdf_read_var(var, filename=data_file)

end

pro _2015_0218_02_load_r_gsm, event_info, varname=var, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = 'orbit_ut'
    times = _2015_0218_02_load_orbit_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'r_gsm'
    if ~cdf_has_var(var, filename=data_file) then begin
        rbsp_read_orbit, time_range, probe=probe, coord='gsm'
        interp_time, var, times
        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end



pro _2015_0218_02_load_model_vars, event_info

    best_model = 't89'
    direction = -1  ; to northern hemisphere.
    igrf = 0
    refine = 1
    event_info['model_setting'] = dictionary($
        'model', best_model, $
        'igrf', igrf , $
        'refine', refine, $
        'direction', -1 )
    data_file = event_info['data_file']
    prefix = event_info['prefix']


    model_vars = [$
        best_model+'_par', $
        prefix+['c0map','fmlt','fmlon','fmlat','bmod_gsm','fpt_gsm']+'_'+best_model]
    load_data = 0
    foreach var, model_vars do begin
        if ~cdf_has_var(var, filename=data_file) then begin
            load_data = 1
            break
        endif
    endforeach
    if load_data then begin
        _2015_0218_02_load_r_gsm, event_info, varname=r_var, time_var=time_var
        read_geopack_info, r_var, model=best_model, $
            direction=direction, refine=refine, igrf=igrf
        times = _2015_0218_02_load_orbit_ut(event_info, varname=time_var)
        foreach var, model_vars do begin
            interp_time, var, times
            data = get_var_data(var, limits=limits)
            cdf_save_var, var, value=data, filename=data_file
            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=var
        endforeach
    endif

    foreach var, model_vars do begin
        if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    endforeach

end



function _2015_0218_02_load_field_ut, event_info, varname=var, $
    time_step=time_step, pad_time=pad_time

    time_range = event_info['time_range']
    data_file = event_info['data_file']

    if n_elements(var) eq 0 then var = 'field_ut'
    if ~cdf_has_var(var, filename=data_file) then begin
        time_step = event_info['field_time_step']
        data_time_range = event_info['field_time_range']
        times = make_bins(data_time_range, time_step)
        cdf_save_var, var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', data_time_range )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    return, cdf_read_var(var, filename=data_file)

end


pro _2015_0218_02_load_b_gsm, event_info, varname=var, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']


    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'b_gsm'
    if ~cdf_has_var(var, filename=data_file) then begin
        rbsp_read_bfield, time_range, probe=probe, coord='gsm', resolution='hires'
        interp_time, var, times
        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end


pro _2015_0218_02_load_b0_gsm, event_info, varname=var, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']


    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'b0_gsm'
    if ~cdf_has_var(var, filename=data_file) then begin
        b_gsm_var = prefix+'b_gsm'
        _2015_0218_02_load_b_gsm, event_info, varname=b_gsm_var, time_var=time_var
        model = event_info['model_setting'].model
        bmod_gsm_var = prefix+'bmod_gsm_'+model
        copy_data, bmod_gsm_var, var
        interp_time, var, to=b_gsm_var
        b0_gsm = get_var_data(var)
        b_gsm = get_var_data(b_gsm_var)
        b1_gsm = b_gsm-b0_gsm
        for ii=0,2 do b0_gsm[*,ii] += mean(b1_gsm[*,ii],/nan)
        get_data, b_gsm_var, limits=limits
        store_data, var, times, b0_gsm, limits=limits

        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end


pro _2015_0218_02_load_b1_gsm, event_info, varname=var, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']


    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'b1_gsm'
    if ~cdf_has_var(var, filename=data_file) then begin
        b_gsm_var = prefix+'b_gsm'
        _2015_0218_02_load_b_gsm, event_info, varname=b_gsm_var, time_var=time_var
        b0_gsm_var = prefix+'b0_gsm'
        _2015_0218_02_load_b0_gsm, event_info, varname=b0_gsm_var, time_var=time_var
        dif_data, b_gsm_var, b0_gsm_var, newname=var
        get_data, b_gsm_var, limits=limits
        store_data, var, limits=limits
        options, var, 'short_name', 'dB'

        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end


pro _2015_0218_02_load_cmap, event_info, varname=cmap_var, time_var=field_time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']
    model = event_info['model_setting'].model

    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'cmap'
    if ~cdf_has_var(var, filename=data_file) then begin
        _2015_0218_02_load_model_vars, event_info
        b0_gsm_var = prefix+'b0_gsm'
        _2015_0218_02_load_b0_gsm, event_info, varname=b0_gsm_var, time_var=time_var
        b0_gsm = get_var_data(b0_gsm_var)
        bmod_gsm_var = prefix+'bmod_gsm_'+model
        bmod_gsm = get_var_data(bmod_gsm_var, at=times, quadratic=1)
        cmap_var = prefix+'c0map_'+model
        c0map = get_var_data(cmap_var, at=times, quadratic=1)
        cmap = c0map*snorm(bmod_gsm)/snorm(b0_gsm)
        store_data, var, times, cmap

        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end


pro _2015_0218_02_load_e_mgse, event_info, varname=var, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']


    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'e_mgse'
    if ~cdf_has_var(var, filename=data_file) then begin
        rbsp_read_efield, time_range, probe=probe, resolution='hires', coord='mgse'
        interp_time, var, times

        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end


pro _2015_0218_02_load_edot0_mgse, event_info, varname=var, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']


    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'edot0_mgse'
    if ~cdf_has_var(var, filename=data_file) then begin
        e_mgse_var = prefix+'e_mgse'
        _2015_0218_02_load_e_mgse, event_info, varname=e_mgse_var, time_var=time_var
        b0_gsm_var = prefix+'b0_gsm'
        _2015_0218_02_load_b0_gsm, event_info, varname=b0_gsm_var, time_var=time_var
        e_mgse = get_var_data(e_mgse_var, limits=limits)
        b0_gsm = get_var_data(b0_gsm_var)
        b0_mgse = cotran(b0_gsm, times, 'gsm2mgse', probe=probe)
        e_mgse[*,0] = -(e_mgse[*,1]*b0_mgse[*,1]+e_mgse[*,2]*b0_mgse[*,2])/b0_mgse[*,0]
        store_data, var, times, e_mgse, limits=limits
        ; b_angle < -20 deg in the time_range, so no filter on b_angle.

        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end


pro _2015_0218_02_load_q_gsm2fac_var, event_info, varname=q_gsm2fac_var, time_var=field_time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']


    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'q_gsm2fac'
    if ~cdf_has_var(var, filename=data_file) then begin
        r_gsm_var = prefix+'r_gsm'
        _2015_0218_02_load_r_gsm, event_info, varname=e_mgse_var, time_var=time_var
        b0_gsm_var = prefix+'b0_gsm'
        _2015_0218_02_load_b0_gsm, event_info, varname=b0_gsm_var, time_var=time_var
        define_fac, b0_gsm_var, r_gsm_var, time_var=b0_gsm_var

        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end


pro _2015_0218_02_load_fac_vars, event_info, varnames=fac_vars, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']


    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(fac_vars) eq 0 then fac_vars = prefix+['b1','e','edot0']+'_fac'
    load_data = 0
    foreach var, fac_vars do begin
        if ~cdf_has_var(var, filename=data_file) then begin
            load_data = 1
            break
        endif
    endforeach
    if load_data then begin
        q_gsm2fac_var = prefix+'q_gsm2fac'
        _2015_0218_02_load_q_gsm2fac_var, event_info, varname=q_gsm2fac_var, time_var=time_var
        b1_gsm_var = prefix+'b1_gsm'
        _2015_0218_02_load_b1_gsm, event_info, varname=b1_gsm_var, time_var=time_var
        b1_fac_var = prefix+'b1_fac'
        to_fac, b1_gsm_var, to=b1_fac_var, q_var=q_gsm2fac_var
        e_mgse_var = prefix+'e_mgse'
        _2015_0218_02_load_e_mgse, event_info, varname=e_mgse_var, time_var=time_var
        edot0_mgse_var = prefix+'edot0_mgse'
        _2015_0218_02_load_edot0_mgse, event_info, varname=edot0_mgse_var, time_var=time_var
        foreach e_var, prefix+['e','edot0'] do begin
            e_mgse = get_var_data(e_var+'_mgse', limits=limits)
            e_gsm = cotran(e_mgse, times, 'mgse2gsm', probe=probe)
            e_gsm_var = e_var+'_gsm'
            store_data, e_gsm_var, times, e_gsm, limits=limits
            options, e_gsm_var, 'coord', 'MGSE'
            e_fac_var = e_var+'_fac'
            to_fac, e_gsm_var, to=e_fac_var, q_var=q_gsm2fac_var
        endforeach

        foreach var, fac_vars do begin
            data = get_var_data(var, limits=limits)
            cdf_save_var, var, value=data, filename=data_file
            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=var
        endforeach
    endif

    foreach var, fac_vars do begin
        if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    endforeach

end


pro _2015_0218_02_load_vsc, event_info, varname=var, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']


    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'vsc'

    if ~cdf_has_var(var, filename=data_file) then begin
        timespan, time_range[0], total(time_range*[-1,1]), /second
        rbsp_load_efw_waveform, probe=probe, type='calibrated', datatype='vsvy', /noclean
        vsc_var = prefix+'efw_vsvy'
        interp_time, vsc_var, times
        vsvy = get_var_data(vsc_var)
        vsc = total(vsvy[*,0:1],2)*0.5
        store_data, var, times, vsc
        add_setting, var, /smart, dictionary($
            'display_type', 'scalar', $
            'short_name', 'Vsc', $
            'unit', 'V', $
            'boom_pair', '12')

        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file


end


pro _2015_0218_02_load_efw_density, event_info, varname=var, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']


    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = prefix+'efw_density'

    if ~cdf_has_var(var, filename=data_file) then begin
        timespan, time_range[0], total(time_range*[-1,1]), /second
        vsc_var = prefix+'vsc'
        _2015_0218_02_load_vsc, event_info, varname=vsc_var, time_var=time_var

        rbsp_efw_density_fit_from_uh_line, vsc_var, probe, newname=var, $
            dmin=1, dmax=3000., setval=!values.f_nan
        add_setting, var, /smart, dictionary($
            'display_type', 'scalar', $
            'short_name', 'EFW n', $
            'unit', 'cm!U-3!N' )

        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end


pro _2015_0218_02_load_pflux_vars, event_info

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = 'field_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    filter = [0.25,1200]    ; sec.
    scale_info = {s0:0.25d, s1:1200d, dj:1d/8, ns:0d}    ; pflux calc.
    event_info['pflux_filter'] = filter
    event_info['pflux_scale_info'] = scale_info

    types = ['','dot0']
    pflux_vars = []
    foreach type, types do begin
        new_vars = prefix+'pf'+type+'_'+['fac','fac_mor_spec_1']
        pflux_vars = [pflux_vars, new_vars]
    endforeach
    load_data = 0
    foreach var, pflux_vars do begin
        if ~cdf_has_var(var, filename=data_file) then begin
            load_data = 1
            break
        endif
    endforeach
    if load_data then begin
        _2015_0218_02_load_fac_vars, event_info, time_var=time_var
        foreach type, types do begin
            e_var = prefix+'e'+type+'_fac'
            b_var = prefix+'b1_fac'
            pf_var = prefix+'pf'+type+'_fac'

            stplot_calc_pflux_mor, e_var, b_var, pf_var, $
                scaleinfo=scale_info
            get_data, pf_var, times
            nrec = n_elements(times)
            ndim = 3
            pffac = fltarr(nrec,ndim)
            for ii=0, ndim-1 do begin
                var = pf_var+'_mor_spec_'+string(ii+1,format='(I0)')
                get_data, var, times, dat, val
                index = where(val ge filter[0] and val le filter[1])
                pffac[*,ii] = total(dat[*,index], 2)
            endfor
            store_data, pf_var, times, pffac
            get_data, e_var, limits=lim
            add_setting, pf_var, /smart, dictionary($
                'short_name', 'S', $
                'unit', 'mW/m', $
                'coord', lim.coord, $
                'coord_labels', lim.coord_labels )


            data = get_var_data(pf_var, limits=limits)
            cdf_save_var, pf_var, value=data, filename=data_file
            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=pf_var

            pf_spec_var = pf_var+'_mor_spec_1'
            get_data, pf_spec_var, times, pfspec, periods
            pf_freq_var = prefix+'pf_spec_freq'
            if ~cdf_has_var(pf_freq_var, filename=data_file) then begin
                freqs = 1d/periods
                cdf_save_var, pf_freq_var, value=freqs, filename=data_file
                settings = dictionary($
                    'var_type', 'support_data', $
                    'unit', 'Hz' )
                cdf_save_setting, settings, filename=data_file, varname=var
            endif
            cdf_save_var, pf_spec_var, value=pfspec, filename=data_file
            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['depend_1'] = pf_freq_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=pf_spec_var
        endforeach
    endif

    foreach var, pflux_vars do begin
        if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    endforeach

end


pro _2015_0218_02_load_hope_velocity, event_info

    ; Use Cristian's data to overwrite rbsp_x_v_gsm and _fac.
    species = ['e','p','o']
    event_info['species'] = species
    data_file = event_info['data_file']
    probe = event_info['probe']
    prefix = event_info['prefix']
    time_range = event_info['time_range']

    root_dir = join_path([srootdir(),prefix+'201502_vel_eflux'])
    files = file_search(join_path([root_dir,'*.tplot']))
    foreach species_name, species do begin
        search_str = species_name
        if species_name eq 'p' then search_str = 'h'
        index = where(stregex(files, '_'+search_str+'_1.tplot') ne -1, count)
        if count eq 0 then stop
        file = files[index]
        tplot_restore, filename=file

        case species_name of
            'e': species_str = 'ele'
            'p': species_str = 'h'
            'o': species_str = 'o'
        endcase
        old_var = 'V_gsm_corr_'+species_str
        new_var = prefix+species_name+'_v_gsm'
        get_data, old_var, times, v_gsm
        get_data, new_var, common_times
        time_step = 22d
        smooth_window = 300d
        smooth_width = smooth_window/time_step
        for ii=0,2 do v_gsm[*,ii] = smooth(v_gsm[*,ii], smooth_width, /nan, /edge_mirror)
;        v_gsm = sinterpol(v_gsm, times, common_times)
;        store_data, new_var, common_times, v_gsm
        store_data, new_var, times, v_gsm
        options, new_var, 'coord', 'GSM'
        fac_var = prefix+species_name+'_v_fac'
        q_gsm2fac_var = prefix+'q_gsm2fac'
        to_fac, new_var, to=fac_var, q_var=q_gsm2fac_var
    endforeach

end



pro _2015_0218_02_load_hope_moments, event_info

    species = ['e','p','o']
    event_info['species'] = species
    data_file = event_info['data_file']
    probe = event_info['probe']
    prefix = event_info['prefix']
    time_range = event_info['time_range']

    hope_vars = [prefix+'ele_n']
    moments = ['en_spec','pa_spec','density','temp','v_gsm']
    foreach species_name, species do begin
        hope_vars = [hope_vars,prefix+species_name+'_'+moments]
    endforeach
    load_data = 0
    foreach var, hope_vars do begin
        if ~cdf_has_var(var, filename=data_file) then begin
            load_data = 1
            break
        endif
    endforeach

    if load_data then begin
        rbsp_read_en_spec, time_range, probe=probe
        rbsp_read_pa_spec, time_range, probe=probe

        foreach species_name, species do begin
            rbsp_read_hope_moments, time_range, probe=probe, species=species_name
            old_var = prefix+species_name+'_t_avg'
            new_var = prefix+species_name+'_temp'
            rename_var, old_var, to=new_var
            old_var = prefix+species_name+'_vbulk'
            new_var = prefix+species_name+'_v_gsm'
            rename_var, old_var, to=new_var
        endforeach

        foreach species_name, species do begin
            species_type = (species_name eq 'e')? 'ele': 'ion'
            time_var = species_type+'_ut'

            vars = prefix+species_name+'_en_spec'
            val_var = species_type+'_en_bins'
            foreach var, vars do begin
                get_data, var, times, data, vals, limits=limits
                if ~cdf_has_var(time_var, filename=data_file) then begin
                    cdf_save_var, time_var, value=times, filename=data_file
                endif
                if ~cdf_has_var(val_var, filename=data_file) then begin
                    cdf_save_var, val_var, value=vals, filename=data_file
                    settings = dictionary($
                        'unit', 'eV', $
                        'var_type', 'support_data' )
                    cdf_save_setting, settings, filename=data_file, varname=var
                endif

                cdf_save_var, var, value=data, filename=data_file
                settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
                settings['depend_0'] = time_var
                settings['depend_1'] = val_var
                settings['var_type'] = 'data'
                cdf_save_setting, settings, filename=data_file, varname=var
            endforeach

            vars = prefix+species_name+'_pa_spec'
            val_var = 'pa_bins'
            foreach var, vars do begin
                get_data, var, times, data, vals, limits=limits
                if ~cdf_has_var(val_var, filename=data_file) then begin
                    cdf_save_var, val_var, value=vals, filename=data_file
                    settings = dictionary($
                        'unit', 'deg', $
                        'var_type', 'support_data' )
                    cdf_save_setting, settings, filename=data_file, varname=var
                endif

                cdf_save_var, var, value=data, filename=data_file
                settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
                settings['depend_0'] = time_var
                settings['depend_1'] = val_var
                settings['var_type'] = 'data'
                cdf_save_setting, settings, filename=data_file, varname=var
            endforeach

            vars = prefix+species_name+'_'+['density','temp','v_gsm']
            foreach var, vars do begin
                data = get_var_data(var, limits=limits)
                cdf_save_var, var, value=data, filename=data_file
                settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
                settings['depend_0'] = time_var
                settings['var_type'] = 'data'
                cdf_save_setting, settings, filename=data_file, varname=var
            endforeach


            rbsp_read_density, time_range, probe=probe
            time_var = 'ele_ut'
            var = prefix+'ele_n'
            interp_time, var, to=prefix+'e_density'
            data = get_var_data(var)
            cdf_save_var, var, value=data, filename=data_file
            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=var
        endforeach
    endif

    foreach var, hope_vars do begin
        if check_if_update(var, time_range, dtime=60) then cdf_load_var, var, filename=data_file
    endforeach

end


pro _2015_0218_02_load_hope_info, event_info

    species = event_info['species']
    data_file = event_info['data_file']
    probe = event_info['probe']
    prefix = event_info['prefix']
    time_range = event_info['hope_info_time_range']

    ; O+ ratio.
    o_dens = get_var_data(prefix+'o_density', in=time_range)
    p_dens = get_var_data(prefix+'p_density', in=time_range)
    o_num_dens_ratio = o_dens/(o_dens+p_dens)
    o_ratio = mean(o_num_dens_ratio)
    event_info['o_num_dens_ratio'] = o_ratio

    ; Number density, cc.
    num_dens = get_var_data(prefix+'efw_density', in=time_range)
    avg_num_dens = mean(num_dens)
    event_info['num_dens'] = avg_num_dens

    ; avg mass.
    p_ratio = 1-o_ratio
    p_mass = 1d
    o_mass = 16d
    avg_mass = p_ratio*p_mass+o_ratio*o_mass

    ; Alfven speed.
    ; MHD: va = B/sqrt(n*mass), mass = SUM_s (n_s*mass_s)/n, n = SUM_s n_s
    va0 = 22.0d     ; km/s, for B in nT, n in cc, m in atomic mass.
    bmag = snorm(get_var_data(prefix+'b0_gsm', in=time_range))
    avg_bmag = mean(bmag)
    event_info['bmag'] = avg_bmag
    event_info['va'] = va0*avg_bmag/sqrt(avg_num_dens*avg_mass)
    event_info['va_p'] = va0*avg_bmag/sqrt(avg_num_dens*p_mass)
    event_info['va_o'] = va0*avg_bmag/sqrt(avg_num_dens*o_mass)


    ; Gyro freq, Omega_i.
    ; MHD: Omega = eB/mass.
    f_g0 = 1.5e-2   ; in Hz.
    event_info['omega'] = f_g0*avg_bmag*avg_mass
    event_info['omega_p'] = f_g0*avg_bmag*p_mass
    event_info['omega_o'] = f_g0*avg_bmag*o_mass

    ; Ion thermal vel, v_i.
    ; MHD: v = sqrt(P/rho), where P = total (P_s), rho = n*mass, P = n*kb*T.
    p_temp = mean(get_var_data(prefix+'p_temp', in=time_range))
    o_temp = mean(get_var_data(prefix+'o_temp', in=time_range))
    event_info['vi_p'] = 310*sqrt(p_temp*1e-3)
    event_info['vi_o'] = 310*sqrt(o_temp*1e-3)
    event_info['vi'] = 310*sqrt((p_temp*p_ratio+o_temp*o_ratio)*1e-3/avg_mass)


    ; Plasma bulk velocity, v_f
    ; MHD: U = (n_s*m_s*U_s)/rho    Eq (2.64)
    p_v_gsm = get_var_data(prefix+'p_v_gsm', in=time_range)
    o_v_gsm = get_var_data(prefix+'o_v_gsm', in=time_range)
    vf_p = mean(snorm(p_v_gsm))
    vf_o = mean(snorm(o_v_gsm))
    event_info['vf_p'] = vf_p
    event_info['vf_o'] = vf_o
    v_gsm = (p_ratio*p_mass*p_v_gsm+o_ratio*o_mass*o_v_gsm)/avg_mass
    event_info['vf'] = mean(snorm(v_gsm))

end


function _2015_0218_02_load_asf_ut, event_info, varname=var, $
    time_step=time_step, pad_time=pad_time

    time_range = event_info['time_range']
    data_file = event_info['data_file']

    if n_elements(var) eq 0 then var = 'asf_ut'
    if ~cdf_has_var(var, filename=data_file) then begin
        time_step = event_info['asf_time_step']
        data_time_range = event_info['asf_time_range']
        times = make_bins(data_time_range, time_step)
        cdf_save_var, var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', data_time_range )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    return, cdf_read_var(var, filename=data_file)

end

pro _2015_0218_02_load_asf_ewo, event_info, varname=var, time_var=time_var

    site = event_info['site']
    data_file = event_info['data_file']

    mlat_range = [61.5,64.5]
    event_info['asf_mlat_range'] = mlat_range

    if n_elements(time_var) eq 0 then time_var = 'asf_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = 'thg_asf_ewo'
    if ~cdf_has_var(var, filename=data_file) then begin
        themis_read_mltimg_ewo, time_range, mlat_range=mlat_range, sites=site
        val_var = 'mlt_bins'
        get_data, var, times, data, val, limits=limits
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['depend_1'] = val_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
        cdf_save_var, val_var, value=val, filename=data_file
        settings = dictionary($
            'var_type', 'support_data', $
            'unit', 'hr', $
            'short_name', 'MLT bins' )
        cdf_save_setting, settings, filename=data_file, varname=val_var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end

pro _2015_0218_02_load_asf_keo, event_info, varname=var, time_var=time_var

    site = event_info['site']
    data_file = event_info['data_file']

    mlt_range = [-2.5,-0.5]
    event_info['asf_mlt_range'] = mlt_range

    if n_elements(time_var) eq 0 then time_var = 'asf_ut'
    times = _2015_0218_02_load_field_ut(event_info, varname=time_var)
    time_range = minmax(times)

    if n_elements(var) eq 0 then var = 'thg_asf_keo'
    if ~cdf_has_var(var, filename=data_file) then begin
        themis_read_mltimg_keo, time_range, mlt_range=mlt_range, sites=site
        val_var = 'mlat_bins'
        get_data, var, times, data, val, limits=limits
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['depend_1'] = val_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
        cdf_save_var, val_var, value=val, filename=data_file
        settings = dictionary($
            'var_type', 'support_data', $
            'unit', 'deg', $
            'short_name', 'MLat bins' )
        cdf_save_setting, settings, filename=data_file, varname=val_var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

end



function _2015_0218_02_load_data_v01, filename=data_file

    if n_elements(data_file) eq 0 then begin
        data_file = join_path([googledir(),'works','works','pflux_grant','data','2015_0218_02_conjunction.cdf'])
    endif

    if file_test(data_file) eq 0 then begin
        cdf_save_setting, 'probe', 'a', filename=data_file
        cdf_save_setting, 'time_range', time_double(['2015-02-18/02:05','2015-02-18/02:15']), filename=data_file
    endif
    event_info = cdf_read_setting(filename=data_file)
    event_info['data_file'] = data_file
    event_info['prefix'] = 'rbsp'+event_info['probe']+'_'
    event_info['site'] = 'gbay'

    prefix = event_info['prefix']
    time_range = event_info['time_range']


    orbit_time_var = 'orbit_ut'
    time_step = 60.
    pad_time = 30.*60
    event_info['orbit_time_var'] = orbit_time_var
    event_info['orbit_time_step'] = time_step
    event_info['orbit_time_range'] = time_range+[-1,1]*pad_time
    times = _2015_0218_02_load_orbit_ut(event_info, $
        time_step=time_step, pad_time=pad_time)
    r_gsm_var = prefix+'r_gsm'
    _2015_0218_02_load_r_gsm, event_info, varname=r_gsm_var, time_var=orbit_time_var
    _2015_0218_02_load_model_vars, event_info


    field_time_var = 'field_ut'
    time_step = 1d/16
    pad_time = 30.*60
    event_info['field_time_var'] = field_time_var
    event_info['field_time_step'] = time_step
    event_info['field_time_range'] = time_range+[-1,1]*pad_time
    times = _2015_0218_02_load_field_ut(event_info, $
        time_step=time_step, pad_time=pad_time)
    b_gsm_var = prefix+'b_gsm'
    _2015_0218_02_load_b_gsm, event_info, varname=b_gsm_var, time_var=field_time_var
    b0_gsm_var = prefix+'b0_gsm'
    _2015_0218_02_load_b0_gsm, event_info, varname=b0_gsm_var, time_var=field_time_var
    cmap_var = prefix+'cmap'
    _2015_0218_02_load_cmap, event_info, varname=cmap_var, time_var=field_time_var
    b1_gsm_var = prefix+'b1_gsm'
    _2015_0218_02_load_b1_gsm, event_info, varname=b1_gsm_var, time_var=field_time_var
    e_mgse_var = prefix+'e_mgse'
    _2015_0218_02_load_e_mgse, event_info, varname=e_mgse_var, time_var=field_time_var
    edot0_mgse_var = prefix+'edot0_mgse'
    _2015_0218_02_load_edot0_mgse, event_info, varname=edot0_mgse_var, time_var=field_time_var

    q_gsm2fac_var = prefix+'q_gsm2fac'
    _2015_0218_02_load_q_gsm2fac_var, event_info, varname=q_gsm2fac_var, time_var=field_time_var
    fac_vars = prefix+['b1','e','edot0']+'_fac'
    _2015_0218_02_load_fac_vars, event_info, varnames=fac_vars, time_var=field_time_var

    _2015_0218_02_load_pflux_vars, event_info

    vsc_var = prefix+'vsc'
    _2015_0218_02_load_vsc, event_info, varname=vsc_var, time_var=field_time_var
    efw_density_var = prefix+'efw_density'
    _2015_0218_02_load_efw_density, event_info, varname=efw_density_var, time_var=field_time_var


    _2015_0218_02_load_hope_moments, event_info
;    event_info['hope_info_time_range'] = time_double(['2015-02-18/02:08','2015-02-18/02:11'])
    event_info['hope_info_time_range'] = event_info['time_range']
    _2015_0218_02_load_hope_velocity, event_info
    _2015_0218_02_load_hope_info, event_info


    asf_time_var = 'asf_ut'
    time_step = 3.
    pad_time = 0.
    event_info['asf_time_var'] = asf_time_var
    event_info['asf_time_step'] = time_step
    event_info['asf_time_range'] = time_range+[-1,1]*pad_time
    times = _2015_0218_02_load_asf_ut(event_info, $
        time_step=time_step, pad_time=pad_time)
    asf_ewo_var = 'thg_asf_ewo'
    _2015_0218_02_load_asf_ewo, event_info, varname=asf_ewo_var, time_var=asf_time_var
    asf_ewo_var = 'thg_asf_keo'
    _2015_0218_02_load_asf_keo, event_info, varname=asf_keo_var, time_var=asf_time_var


    return, event_info

end


tmp = _2015_0218_02_load_data()
end
