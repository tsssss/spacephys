function alfven_arc_load_rbsp_data_pflux_setting, event_info

    key = 'pflux_setting'
    if ~event_info.haskey(key) then begin
        filter = [0.25d,1800]    ; sec.
;        filter = [10d,1800]    ; sec.
        scale_info = {s0:min(filter), s1:max(filter), dj:1d/8, ns:0d}
        event_info[key] = dictionary($
            'filter', filter, $
            'scale_info', scale_info )
    endif
    return, event_info[key]

end


function alfven_arc_load_rbsp_data_field_ut, event_info, time_var=var, get_name=get_name

    if n_elements(var) eq 0 then var = 'field_ut'
    if keyword_set(get_name) then return, var

    data_file = event_info['data_file']
    if ~cdf_has_var(var, filename=data_file) then begin
        time_range = event_info['field_time_range']
        time_step = event_info['field_time_step']
        times = make_bins(time_range, time_step)
        cdf_save_var, var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', time_range )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    return, cdf_read_var(var, filename=data_file)

end

function alfven_arc_load_rbsp_data_orbit_ut, event_info, time_var=var, get_name=get_name

    if n_elements(var) eq 0 then var = 'orbit_ut'
    if keyword_set(get_name) then return, var

    data_file = event_info['data_file']
    if ~cdf_has_var(var, filename=data_file) then begin
        time_range = event_info['orbit_time_range']
        time_step = event_info['orbit_time_step']
        times = make_bins(time_range, time_step)
        cdf_save_var, var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', time_range )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    return, cdf_read_var(var, filename=data_file)

end



function alfven_arc_load_rbsp_data_hope_moments, event_info

    species = ['e','p','o']
    event_info['species'] = species
    data_file = event_info['data_file']
    probe = event_info['probe']
    prefix = event_info['prefix']
    time_range = event_info['time_range']
    secofday = constant('secofday')
    time_range = time_range[0]-(time_range[0] mod secofday)+[0,secofday]

    hope_vars = []
    coord = 'gsm'
    moments = ['en_spec'+['','_'+['anti','perp','para']],'n','t',$
        'density','temp',['vbulk','nflux','eflux','enthalpy']+'_'+coord]
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
        foreach species_name, species do begin
            en_spec_combo_vars = rbsp_read_en_spec_combo(time_range, probe=probe, species=species_name)
            en_spec_var = rbsp_read_en_spec(time_range, probe=probe, species=species_name)

        ;---Energy spec vars.
            vars = [en_spec_var,(en_spec_combo_vars.values()).toarray()]
            species_type = (species_name eq 'e')? 'ele': 'ion'
            time_var = species_type+'_ut'
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

        ;---Moment vars.
            mom_vars = rbsp_read_hope_moments(time_range, probe=probe, species=species_name, coord=coord)

            foreach var, mom_vars do begin
                get_data, var, times, data, limits=limits
                if ~cdf_has_var(time_var, filename=data_file) then begin
                    cdf_save_var, time_var, value=times, filename=data_file
                endif
                common_times = cdf_read_var(time_var, filename=data_file)
                if n_elements(common_times) ne n_elements(times) then begin
                    interp_time, var, common_times
                    data = get_var_data(var)
                endif
                cdf_save_var, var, value=data, filename=data_file
                settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
                settings['depend_0'] = time_var
                settings['var_type'] = 'data'
                cdf_save_setting, settings, filename=data_file, varname=var
            endforeach
        
        ;---HOPE L3 density and temperature.
            dens_var = rbsp_read_density_hope(time_range, probe=probe, species=species_name)
            temp_var = rbsp_read_temperature(time_range, probe=probe, species=species_name)
            foreach var, [dens_var,temp_var] do begin
                get_data, var, times, data, limits=limits
                if ~cdf_has_var(time_var, filename=data_file) then begin
                    cdf_save_var, time_var, value=times, filename=data_file
                endif
                cdf_save_var, var, value=data, filename=data_file
                settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
                settings['depend_0'] = time_var
                settings['var_type'] = 'data'
                cdf_save_setting, settings, filename=data_file, varname=var
            endforeach
        endforeach
    endif

    foreach var, hope_vars do begin
        if check_if_update(var, time_range, dtime=60) then cdf_load_var, var, filename=data_file
    endforeach

    return, hope_vars
end

function alfven_arc_load_rbsp_data_density, event_info

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']
    time_range = event_info['time_range']
    secofday = constant('secofday')
    time_range = time_range[0]-(time_range[0] mod secofday)+[0,secofday]

    density_emfisis_var = rbsp_read_density_emfisis(time_range, probe=probe, get_name=1)
    density_efw_var = rbsp_read_density_efw(time_range, probe=probe, get_name=1)
    density_hope_var = prefix+'density_hope'

    vars = [density_emfisis_var,density_efw_var,density_hope_var]
    labels = ['EMFISIS','EFW','HOPE']
    colors = constant('rgb')

    time_var = 'density_ut'
    if ~cdf_has_var(time_var, filename=data_file) then begin
        time_step = 6
        times = make_bins(time_range, time_step)
        
        cdf_save_var, time_var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', time_range )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif


    foreach var, vars, var_id do begin
        if ~cdf_has_var(var, filename=data_file) then begin
            routine = 'rbsp_read_density_'+strlowcase(labels[var_id])
            if strlowcase(labels[var_id]) eq 'hope' then begin
                var = call_function(routine, time_range, probe=probe, species='e')
                var = rename_var(var, output=density_hope_var)
            endif else begin
                var = call_function(routine, time_range, probe=probe)
            endelse
            if n_elements(times) eq 0 then times = cdf_read_var(time_var, filename=data_file)
            interp_time, var, times
            
            cdf_save_var, var, value=get_var_data(var, limits=limits), filename=data_file
            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=var
        endif
        if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    endforeach

    var = prefix+'density_combo'
    var = stplot_merge(vars, output=var)
    add_setting, var, smart=1, dictionary($
        'ylog', 1, $
        'display_type', 'stack', $
        'unit', 'cm!E-3!N', $
        'labels', 'N '+labels, $
        'colors', colors )
    return, var

end


function alfven_arc_load_rbsp_data_r_gsm, event_info, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_rbsp_data_orbit_ut(get_name=1)
    times = alfven_arc_load_rbsp_data_orbit_ut(event_info)
    time_range = event_info['orbit_time_range']

    coord = 'gsm'
    var = rbsp_read_orbit(time_range, probe=probe, coord=coord, get_name=1)
    
    if ~cdf_has_var(var, filename=data_file) then begin
        var = rbsp_read_orbit(time_range+[-1,1]*60, probe=probe, coord=coord)
        interp_time, var, times
        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    
    return, var

end

function alfven_arc_load_rbsp_data_load_model_setting, event_info

    key = 'model_setting'
    if ~event_info.haskey(key) then begin
        model_setting = dictionary($
            'external_model', 't89', $
            'models', ['t89','t96','t01','t04s'], $
            't89_par', 2, $
            'igrf', 0, $
            'refine', 1, $
            'direction', 'north' )
        internal_model = (model_setting['igrf'] eq 0)? 'dipole': 'igrf'
        model_setting['internal_model'] = internal_model
        event_info[key] = model_setting
    endif
    return, event_info[key]

end


function alfven_arc_load_rbsp_data_load_model_vars, event_info, time_var=time_var

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_rbsp_data_orbit_ut(get_name=1)
    data_file = event_info['data_file']
    prefix = event_info['prefix']

    model_setting = alfven_arc_load_rbsp_data_load_model_setting(event_info)
    models = model_setting['models']
    refine = model_setting['refine']
    all_model_vars = []
    

    foreach direction, ['south','north'] do begin
        if direction eq 'north' then begin
            north = 1
            south = 0
        endif else begin
            north = 0
            south = 1
        endelse
        
        
        model_vars = []
        vars = prefix+['fmlt','fmlon','fmlat','bf_gsm_'+models]
        foreach internal_model, ['igrf','dipole'] do begin
            load_data = 0
            suffix = '_'+internal_model+'_'+direction
;            if direction eq 'south' then suffix += '_'+direction
            
            the_vars = vars+suffix
            model_vars = [model_vars,the_vars]
            foreach var, the_vars do begin
                if ~cdf_has_var(var, filename=data_file) then begin
                    load_data = 1
                    break
                endif
            endforeach
            if load_data then begin
                r_var = alfven_arc_load_rbsp_data_r_gsm(event_info, time_var=time_var)
                igrf = (internal_model eq 'igrf')? 1: 0
                vinfo = geopack_trace_to_ionosphere(r_var, models=models, $
                    igrf=igrf, south=south, north=north, refine=refine, suffix=suffix)
                foreach var, the_vars do begin
                    if tnames(var) eq '' then stop
                    data = get_var_data(var, limits=limits)
                    cdf_save_var, var, value=data, filename=data_file
                    settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
                    settings['depend_0'] = time_var
                    settings['var_type'] = 'data'
                    cdf_save_setting, settings, filename=data_file, varname=var
                endforeach
            endif
        endforeach
        
        foreach var, model_vars do begin
            if check_if_update(var, time_range) then begin
                cdf_load_var, var, filename=data_file
                get_data, var, times, data
                ndim = n_elements(data)/n_elements(times)
                window = 600d
                time_step = sdatarate(times)
                width = window/time_step
                if ndim eq 1 then begin
                    data = smooth(data, width, edge_mirror=1)
                endif else begin
                    for ii=0,ndim-1 do begin
                        data[*,ii] = smooth(data[*,ii], width, edge_mirror=1)
                    endfor
                endelse
                store_data, var, times, data
            endif
        endforeach
        
        all_model_vars = [all_model_vars,model_vars]
    endforeach
    



    return, all_model_vars

end

function alfven_arc_load_rbsp_data_bmod_gsm, event_info, time_var=time_var

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_rbsp_data_orbit_ut(get_name=1)
    data_file = event_info['data_file']
    prefix = event_info['prefix']

    model_setting = alfven_arc_load_rbsp_data_load_model_setting(event_info)
    models = model_setting['models']
    t89_par = model_setting['t89_par']
    coord = 'gsm'

    bmod_vars = []
    vars = prefix+'bmod_gsm_'+models
    foreach internal_model, ['igrf','dipole'] do begin
        load_data = 0
        the_vars = vars+'_'+internal_model
        bmod_vars = [bmod_vars,the_vars]
        foreach var, the_vars do begin
            if ~cdf_has_var(var, filename=data_file) then begin
                load_data = 1
                break
            endif
        endforeach
        if load_data then begin
            r_var = alfven_arc_load_rbsp_data_r_gsm(event_info, time_var=time_var)
            igrf = (internal_model eq 'igrf')? 1: 0
            vinfo = geopack_read_bfield(r_var, models=models, igrf=igrf, suffix='_'+internal_model, t89_par=t89_par, coord=coord)
            
            foreach var, the_vars do begin
                if tnames(var) eq '' then stop
                data = get_var_data(var, limits=limits)
                cdf_save_var, var, value=data, filename=data_file
                settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
                settings['depend_0'] = time_var
                settings['var_type'] = 'data'
                cdf_save_setting, settings, filename=data_file, varname=var
            endforeach
        endif
    endforeach

    foreach var, bmod_vars do begin
        if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    endforeach

    return, bmod_vars

end


function alfven_arc_load_rbsp_data_b_gsm, event_info, time_var=field_time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_rbsp_data_field_ut(get_name=1)
    times = alfven_arc_load_rbsp_data_field_ut(event_info)
    time_range = event_info['field_time_range']

    coord = 'gsm'
    resolution = 'hires'
    var = rbsp_read_bfield(time_range, probe=probe, coord=coord, resolution=resolution, get_name=1)
    if ~cdf_has_var(var, filename=data_file) then begin
        var = rbsp_read_bfield(time_range, probe=probe, coord=coord, resolution=resolution)
        interp_time, var, times
        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif

    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

    return, var
end



function alfven_arc_load_rbsp_data_b0_gsm, event_info, time_var=time_var

    prefix = event_info['prefix']
    data_file = event_info['data_file']
    time_range = event_info['field_time_range']

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_rbsp_data_field_ut(get_name=1)

    var = prefix+'b0_gsm'
    if ~cdf_has_var(var, filename=data_file) then begin
        b_gsm_var = alfven_arc_load_rbsp_data_b_gsm(event_info, time_var=time_var)
        bmod_gsm_vars = alfven_arc_load_rbsp_data_bmod_gsm(event_info)

        model = 't89'
        igrf = 0
        internal_model = (igrf eq 1)? 'igrf': 'dipole'
        suffix = '_'+internal_model
        bmod_gsm_var = prefix+'bmod_gsm_'+model+suffix

        b_gsm = get_var_data(b_gsm_var, times=times)
        bmod_gsm = get_var_data(bmod_gsm_var, at=times)
        b1_gsm = b_gsm-bmod_gsm
        ndim = 3
        time_step = total(times[0:1]*[-1,1])
        window = event_info['b0_window']
        width = window/time_step
        for ii=0,ndim-1 do begin
            b1_gsm[*,ii] -= smooth(b1_gsm[*,ii], width, edge_mirror=1, nan=1)
        endfor

        b0_gsm = b_gsm-b1_gsm
        store_data, var, times, b0_gsm
        add_setting, var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', strupcase(model)+' B', $
            'unit', 'nT', $
            'coord', 'GSM', $
            'coord_labels', constant('xyz'), $
            'model', model, $
            'internal_model', internal_model )


        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

    return, var

end

function alfven_arc_load_rbsp_data_b1_gsm, event_info

    prefix = event_info['prefix']
    data_file = event_info['data_file']

    b_gsm_var = alfven_arc_load_rbsp_data_b_gsm(event_info)
    b0_gsm_var = alfven_arc_load_rbsp_data_b0_gsm(event_info)
    b_gsm = get_var_data(b_gsm_var, times=times)
    b0_gsm = get_var_data(b0_gsm_var, at=times)
    b1_gsm = b_gsm-b0_gsm

    b1_gsm_var = prefix+'b1_gsm'
    store_data, b1_gsm_var, times, b1_gsm
    add_setting, b1_gsm_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'dB', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', constant('xyz') )
    return, b1_gsm_var

end



function alfven_arc_load_rbsp_data_e_mgse_2015_0416, event_info, time_var=field_time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_rbsp_data_field_ut(get_name=1)
    times = alfven_arc_load_rbsp_data_field_ut(event_info)
    time_range = event_info['field_time_range']

    if probe eq 'a' then pair = '24' else pair = '12'
    ; rbspa_e_mgse_[v24,spinfit_phasef] and rbspb_e_mgse_[survey,spinfit_phasef]
    if probe eq 'a' then begin
        ; V1 is bad, V2 and V4 are good, V3 is mostly good around the time of the event (but bad before).
        vsvy_var = rbsp_read_vsc(time_range, probe=probe, id='spin_plane')
        interp_time, vsvy_var, times
        vsvy = get_var_data(vsvy_var)
        rbsp_efw_read_boom_flag, time_range, probe=probe
        flag_var = prefix+'boom_flag'
        interp_time, flag_var, times
        boom_flags = get_var_data(flag_var) ne 0
        for ii=0,3 do begin
            bad_index = where(boom_flags[*,ii] eq 0, count)
            if count ne 0 then vsvy[bad_index,ii] = !values.f_nan
        endfor

        ; Get vsc.
        vsc = (vsvy[*,2]+vsvy[*,3])*0.5
        spin_period = 11d   ; the rough number works fine, no need to get the accurate number
        dt = sdatarate(times)
        width = spin_period/dt
        vsc = smooth(vsc, width, nan=1, edge_zero=1)

        cp0 = rbsp_efw_get_cal_params(time_range[0])
        cp = (probe eq 'a')? cp0.a: cp0.b
        boom_length = cp.boom_length
        boom_shorting_factor = cp.boom_shorting_factor
        the_boom_length = boom_length*boom_shorting_factor
        ntime = n_elements(times)

        ; Try to use V2,V3,V4.
        ; calc e_uvw
        eu = (vsc-vsvy[*,1])/(the_boom_length[0]*0.5)*1000
        ev = (vsvy[*,2]-vsvy[*,3])/the_boom_length[1]*1000
        ew = fltarr(ntime)
        e_uvw = [[eu],[ev],[ew]]
        for ii=0,1 do begin
            data = e_uvw[*,ii]
            offset1 = smooth(data, width, nan=1, edge_zero=1)
            offset2 = smooth(offset1, width, nan=1, edge_zero=1)
            data -= offset2
            e_uvw[*,ii] = data
        endfor
        var = prefix+'e_uvw_v234'
        store_data, var, times, e_uvw
        add_setting, var, smart=1, dictionary('coord','uvw', 'coord_labels',constant('uvw'), 'colors',constant('rgb'))
        e_mgse_var = prefix+'e_mgse_v234'
        e_mgse = cotran(e_uvw, times, 'uvw2mgse', probe=probe)
        store_data, e_mgse_var, times, e_mgse
        add_setting, e_mgse_var, smart=1, id='efield', dictionary('coord', 'mgse')

        ; Try to use V2 and V4 only.
        ; calc e_uvw
        ev = (vsvy[*,2]-vsvy[*,3])/the_boom_length[1]*1000
        eu = shift(ev,width*0.25)
        ew = fltarr(ntime)
        e_uvw = [[eu],[ev],[ew]]
        for ii=0,1 do begin
            data = e_uvw[*,ii]
            offset1 = smooth(data, width, nan=1, edge_zero=1)
            offset2 = smooth(offset1, width, nan=1, edge_zero=1)
            data -= offset2
            e_uvw[*,ii] = data
        endfor
        var = prefix+'e_uvw_v24'
        store_data, var, times, e_uvw
        add_setting, var, smart=1, dictionary('coord','uvw', 'coord_labels',constant('uvw'), 'colors',constant('rgb'))
        e_mgse = cotran(e_uvw, times, 'uvw2mgse', probe=probe)
        e_mgse_var = prefix+'e_mgse_v24'
        store_data, e_mgse_var, times, e_mgse
        add_setting, e_mgse_var, smart=1, id='efield', dictionary('coord', 'mgse')


        ; Try spinfit.
        e_mgse_var = rbsp_read_efield_spinfit_phasef(time_range, probe=probe, id=pair)

        e_mgse_vars = prefix+'e_mgse_'+['spinfit_phasef','v'+pair]
    endif else begin
        ; Try spinfit
        e_mgse_var = rbsp_read_efield_spinfit_phasef(time_range, probe=probe, id=pair)

        ; Load survey resolutions efield.
        e_mgse_var = rbsp_read_efield_survey(time_range, probe=probe)
        e_mgse_var = rename_var(e_mgse_var,output=prefix+'e_mgse_survey')

        e_mgse_vars = prefix+'e_mgse_'+['spinfit_phasef','survey']
    endelse


    ; Calculate edot0.
    b0_gsm_var = alfven_arc_load_rbsp_data_b0_gsm(event_info)
    b0_mgse = cotran(get_var_data(b0_gsm_var), times, 'gsm2mgse', probe=probe)
    edot0_angle_var = prefix+'edot0_angle'
    edot0_angle = asin(b0_mgse[*,0]/snorm(b0_mgse))*constant('deg')
    store_data, edot0_angle_var, times, edot0_angle
    add_setting, edot0_angle_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'edot0 angle', $
        'unit', 'deg' )

    edot0_vars = e_mgse_vars
    foreach e_var, e_mgse_vars, var_id do begin
        edot0_var = streplace(e_var,'e_mgse','edot0_mgse')
        edot0_vars[var_id] = edot0_var
        interp_time, e_var, times
        e_mgse = get_var_data(e_var)
        edot0_mgse = e_mgse
        edot0_mgse[*,0] = -total(e_mgse[*,1:2]*b0_mgse[*,1:2],2)/b0_mgse[*,0]

        var = edot0_var
        store_data, var, times, edot0_mgse
        add_setting, var, smart=1, id='efield', dictionary('coord', 'mgse')
    endforeach

    e_mgse_vars = edot0_vars
    e_gsm_vars = e_mgse_vars
    foreach mgse_var, e_mgse_vars, var_id do begin
        gsm_var = streplace(mgse_var, 'mgse', 'gsm')
        e_gsm_vars[var_id] = gsm_var
        vec_mgse = get_var_data(mgse_var, times=times)
        vec_gsm = cotran(vec_mgse, times, 'mgse2gsm', probe=probe)
        store_data, gsm_var, times, vec_gsm
        add_setting, gsm_var, smart=1, id='efield', dictionary('coord', 'gsm')
    endforeach

    ; Calc Poynting flux.
    r_gsm_var = prefix+'r_gsm'
    define_fac, b0_gsm_var, r_gsm_var, time_var=b0_gsm_var

    b1_gsm_var = alfven_arc_load_rbsp_data_b1_gsm(event_info)
    b1_fac_var = streplace(b1_gsm_var, 'gsm', 'fac')
    to_fac, b1_gsm_var, to=b1_fac_var

    e_fac_vars = e_gsm_vars
    foreach e_gsm_var, e_gsm_vars, var_id do begin
        e_fac_var = streplace(e_gsm_var, 'gsm', 'fac')
        e_fac_vars[var_id] = e_fac_var
        to_fac, e_gsm_var, to=e_fac_var
    endforeach

    model_setting = event_info['model_setting']
    external_model = model_setting['external_model']
    internal_model = model_setting['internal_model']
    bf_var = prefix+'bf_gsm_'+external_model+'_'+internal_model+'_south'
    b0_var = prefix+'b0_gsm'
    b0_gsm = get_var_data(b0_var, times=times)
    bf_gsm = get_var_data(bf_var, at=times)
    cmap = snorm(bf_gsm)/snorm(b0_gsm)


    pflux_setting = alfven_arc_load_rbsp_data_pflux_setting(event_info)
    scale_info = pflux_setting['scale_info']

    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    event_info['fac_labels'] = fac_labels
    pf_fac_vars = e_fac_vars
    foreach e_fac_var, e_fac_vars, var_id do begin
        pf_fac_var = streplace(e_fac_var, 'edot0', 'pfdot0')
        pf_fac_vars[var_id] = pf_fac_var
        stplot_calc_pflux_mor, e_fac_var, b1_fac_var, pf_fac_var, scaleinfo=scale_info
        add_setting, pf_fac_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'S', $
            'unit', 'mW/m!E2!N', $
            'coord', 'FAC', $
            'coord_labels', fac_labels )

        ; Normalize to 100 km.
        pf_fac_map_var = pf_fac_var+'_map'
        pf_fac = get_var_data(pf_fac_var, times=times)
        pf_fac_map = pf_fac
        ndim = 3
        for ii=0,ndim-1 do pf_fac_map[*,ii] = cmap*pf_fac[*,ii]
        store_data, pf_fac_map_var, times, pf_fac_map
        add_setting, pf_fac_map_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'S', $
            'unit', 'mW/m!E2!N', $
            'model', model, $
            'internal_model', internal_model, $
            'coord', 'FAC', $
            'coord_labels', event_info['fac_labels'] )
        pf_fac_vars = [pf_fac_vars,pf_fac_map_var]

        ; pflux spec.
        spec_unit = tex2str('mu')+'W/m!E2!N'
        spec_ct = 66
        spec_zrange = [-1,1]*10
        spec_zstep = 10
        spec_zminor = 5
        spec_ztickv = make_bins(spec_zrange, spec_zstep, inner=1)
        spec_zticks = n_elements(spec_ztickv)-1

        pf_spec_var = pf_fac_var+'_mor_spec_1'
        get_data, pf_spec_var, times, pfspec, ps
        fs = 1d3/ps
        store_data, pf_spec_var, times, pfspec*1e3, fs
        add_setting, pf_spec_var, smart=1, dictionary($
            'display_type', 'spec', $
            'short_name', 'S!D||!N', $
            'unit', spec_unit, $
            'ytitle', 'Freq!C(mHz)', $
            'yrange', minmax(fs), $
            'ylog', 1, $
            'zlog', 0, $
            'zrange', spec_zrange, $
            'ztickv', spec_ztickv, $
            'zticks', spec_zticks, $
            'zminor', spec_zminor, $
            'zticklen', spec_zticklen, $
            'color_table', spec_ct )


        ; Calculate ExB velocity.
        b0_fac = b0_gsm
        b0_fac[*,0] = snorm(b0_gsm)
        b0_fac[*,1:2] = 0
        coef_fac = b0_fac
        coef_fac[*,0] /= (b0_fac[*,0])^2
        e_fac = get_var_data(e_fac_var, at=times)
        vexb_fac = vec_cross(e_fac,coef_fac)*1e3
        vexb_fac_var = streplace(e_fac_var, 'edot0','vexbdot0')
        store_data, vexb_fac_var, times, vexb_fac
        add_setting, vexb_fac_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'ExB V', $
            'unit', 'km/s', $
            'coord', 'FAC', $
            'coord_labels', event_info['fac_labels'] )
    endforeach

    return, pf_fac_vars
end



function alfven_arc_load_rbsp_data_plasma_param, event_info

    data_file = event_info['data_file']
    probe = event_info['probe']
    prefix = event_info['prefix']
    the_time_range = event_info['plasma_param_time_range']
    plasma_param = dictionary()

    ; Number density, cc.
    num_dens = get_var_data(prefix+'e_density', in=the_time_range)
    avg_num_dens = mean(num_dens,nan=1)
    plasma_param['e_num_dens'] = avg_num_dens

    ; O+ ratio.
    o_dens = get_var_data(prefix+'o_density', in=the_time_range)
    p_dens = get_var_data(prefix+'p_density', in=the_time_range)
    o_num_dens_ratio = o_dens/(o_dens+p_dens)
    o_ratio = mean(o_num_dens_ratio,nan=1)  ; n_o/n_p.
    plasma_param['o_num_dens_partial'] = mean(o_dens,nan=1)
    plasma_param['p_num_dens_partial'] = mean(p_dens,nan=1)
    plasma_param['o_num_dens_ratio'] = o_ratio
    plasma_param['r_o'] = 1/(1+1/o_ratio)   ; n_o/(n_o+n_p).

    ; Use HOPE density?
    plasma_param['num_dens'] = plasma_param['e_num_dens']
    avg_num_dens = plasma_param['num_dens']


    ; avg mass.
    p_ratio = 1-o_ratio
    p_mass = 1d
    o_mass = 16d
    avg_mass = p_ratio*p_mass+o_ratio*o_mass
    plasma_param['avg_ion_mass'] = avg_mass

    ; Alfven speed.
    ; MHD: va = B/sqrt(n*mass), mass = SUM_s (n_s*mass_s)/n, n = SUM_s n_s
    va0 = 22.0d     ; km/s, for B in nT, n in cc, m in atomic mass.
    va0 = 1e-9/sqrt(1e6*!dpi*4e-7*1.67e-27)*1e-3    ; km/s, B in nT, n in cc, m in m_p.
    bmag = snorm(get_var_data(prefix+'b0_gsm', in=the_time_range))
    avg_bmag = mean(bmag)
    plasma_param['bmag'] = avg_bmag
    plasma_param['va'] = va0*avg_bmag/sqrt(avg_num_dens*avg_mass)
    plasma_param['va_p'] = va0*avg_bmag/sqrt(avg_num_dens*p_mass)
    plasma_param['va_o'] = va0*avg_bmag/sqrt(avg_num_dens*o_mass)

    ; Gyro freq, Omega_i.
    ; MHD: Omega = eB/mass, f = eB/(2*pi*mass).
    f_g0 = 1.6e-19*1e-9/1.67e-27/2/!dpi   ; in Hz.
    plasma_param['omega'] = f_g0*avg_bmag/avg_mass*2*!dpi
    plasma_param['omega_p'] = f_g0*avg_bmag/p_mass*2*!dpi
    plasma_param['omega_o'] = f_g0*avg_bmag/o_mass*2*!dpi
    plasma_param['fg'] = f_g0*avg_bmag/avg_mass
    plasma_param['fg_p'] = f_g0*avg_bmag/p_mass
    plasma_param['fg_o'] = f_g0*avg_bmag/o_mass

    ; Ion thermal vel, v_i.
    ; MHD: v = sqrt(P/rho), where P = total (P_s), rho = n*mass, P = n*kb*T.
    p_temp = mean(get_var_data(prefix+'p_temp', in=the_time_range),nan=1)
    o_temp = mean(get_var_data(prefix+'o_temp', in=the_time_range),nan=1)
    e_temp = mean(get_var_data(prefix+'e_temp', in=the_time_range),nan=1)
    plasma_param['t_p'] = p_temp
    plasma_param['t_o'] = o_temp
    plasma_param['t_e'] = e_temp
    plasma_param['vi_p'] = 310*sqrt(p_temp*1e-3)
    plasma_param['vi_o'] = 310*sqrt(o_temp*1e-3)
    avg_i_temp = (p_temp*p_ratio+o_temp*o_ratio)
    plasma_param['vi'] = 310*sqrt(avg_i_temp*1e-3/avg_mass)
    plasma_param['t_i'] = avg_i_temp
    plasma_param['v_e'] = 310*sqrt(e_temp*1e-3*1836)
    
    
    ; Ion acoustic gyroradius.
    ; rho^2 = kb*T_e/(m_i*f_gi^2)
    ;       = (kb*T_i/m_i)*T_e/T_i/f_gi^2
    rho = plasma_param['vi']*sqrt(plasma_param['t_e']/plasma_param['t_i'])/plasma_param['fg']
    plasma_param['rho'] = rho

    ; Plasma bulk velocity, v_f
    ; MHD: U = (n_s*m_s*U_s)/rho    Eq (2.64)
    p_v_gsm = get_var_data(prefix+'p_vbulk_gsm', in=the_time_range, times=times)
    o_v_gsm = get_var_data(prefix+'o_vbulk_gsm', in=the_time_range)
    vf_p = mean(snorm(p_v_gsm))
    vf_o = mean(snorm(o_v_gsm))
    plasma_param['vf_p'] = vf_p
    plasma_param['vf_o'] = vf_o
    v_gsm = (p_ratio*p_mass*p_v_gsm+o_ratio*o_mass*o_v_gsm)/avg_mass
    store_data, prefix+'vf_gsm', times, v_gsm
    add_setting, prefix+'vf_gsm', smart=1, dictionary($
        'display_type', 'vector', $
        'unit', 'km/s', $
        'short_name', 'MHD V', $
        'coord', 'GSM', $
        'coord_labels', ['x','y','z'] )
    q_gsm2fac_var = prefix+'q_gsm2fac'
    vf_fac_var = prefix+'vf_fac'
    to_fac, prefix+'vf_gsm', to=vf_fac_var, q_var=q_gsm2fac_var
    var = prefix+'vf_fac'
    add_setting, vf_fac_var, smart=1, dictionary($
        'display_type', 'vector', $
        'unit', 'km/s', $
        'short_name', 'MHD V', $
        'coord', 'FAC', $
        'coord_labels', event_info['fac_labels'] )
    get_data, vf_fac_var, times, v_fac, limits=lim
    ndim = 3
    vf_fac = fltarr(ndim)
;    for ii=0,ndim-1 do vf_fac[ii] = abs(mean(v_fac[*,ii],nan=1))
    for ii=0,ndim-1 do vf_fac[ii] = mean(abs(v_fac[*,ii]),nan=1)
    plasma_param['vf_fac'] = vf_fac
    plasma_param['vf'] = vf_fac[1]
    
;    vf_exb_fac = fltarr(ndim)
;    vexb_fac = get_var_data(prefix+'vexbdot0_fac', in=the_time_range)
;    for ii=0,ndim-1 do vf_fac[ii] = mean(abs(vexb_fac[*,ii]),nan=1)
;    plasma_param['vf_exb_fac'] = vf_fac
;    plasma_param['vf_exb'] = vf_fac[1]

    ; Plasma beta.
    mu0 = !dpi*4e-7
    p_b = (plasma_param['bmag'])^2/(2*mu0)*1e-9   ; nPa.
    plasma_param['p_mag'] = p_b
    coef = 1.6e-19*1e9*1e6
    p_p = mean(get_var_data(prefix+'p_temp', in=the_time_range)*$
        get_var_data(prefix+'p_density', in=the_time_range))*coef
    p_o = mean(get_var_data(prefix+'o_temp', in=the_time_range)*$
        get_var_data(prefix+'o_density', in=the_time_range))*coef
    p_e = mean(get_var_data(prefix+'e_temp', in=the_time_range)*$
        get_var_data(prefix+'e_density', in=the_time_range))*coef
    plasma_param['p_p'] = p_p
    plasma_param['p_o'] = p_o
    plasma_param['p_e'] = p_e
    plasma_param['beta_p'] = p_p/p_b
    plasma_param['beta_o'] = p_o/p_b
    plasma_param['beta_e'] = p_e/p_b
    plasma_param['beta'] = (p_p+p_o+p_e)/p_b


    ; Ion acoustic gyroradius.
    tt = plasma_param['t_e']+plasma_param['t_p']+plasma_param['t_o']
    tt = plasma_param['t_e']
    mm = plasma_param['avg_ion_mass']
    c_s = sqrt(tt/mm)*sqrt(1.6e-19/1.67e-27) ; m/s.
    w_i = plasma_param['omega']/2/!dpi
    rho_i = c_s/w_i*1e-3    ; km.
    plasma_param['rho_i'] = rho_i

    ; Gyro radius.
    r_p = plasma_param['vi_p']*1.67e-27/1.6e-19/plasma_param['bmag']*1e3/1e-9*1e-3  ; km.
    r_o = plasma_param['vi_o']*16*1.67e-27/1.6e-19/plasma_param['bmag']*1e3/1e-9*1e-3  ; km.
    plasma_param['r_p'] = r_p
    plasma_param['r_o'] = r_o

;    the_time = event_info['snapshot_time']
;    r_gsm = get_var_data(prefix+'r_gsm', at=the_time)
;    plasma_param['r_gsm'] = r_gsm
;    plasma_param['r_sm'] = cotran(r_gsm, the_time, 'gsm2sm')

    event_info['plasma_param'] = plasma_param
    return, plasma_param

end



function alfven_arc_load_rbsp_data_2015_0416, input_time_range, probe=probe, filename=data_file, version=version, _extra=ex

    time_range = time_double(input_time_range)
    if n_elements(probe) eq 0 then message, 'No input probe ...'
    
    if n_elements(version) eq 0 then version = 'v01'
    if n_elements(data_file) eq 0 then begin
        base = time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_rbsp'+probe+'_data_'+version+'.cdf'
        data_file = join_path([googledir(),'works','pflux_grant','alfven_arc','data',base])
    endif
    

    if file_test(data_file) eq 0 then begin
        cdf_save_setting, 'probe', probe, filename=data_file
        cdf_save_setting, 'prefix', 'rbsp'+probe+'_', filename=data_file
        cdf_save_setting, 'time_range', time_range, filename=data_file
    endif
    event_info = cdf_read_setting(filename=data_file)
    event_info['data_file'] = data_file

    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']
    event_info['sc_color'] = sgcolor('dark_violet')
    
    
    
;---Orbit related vars.
    orbit_time_var = 'orbit_ut'
    time_step = 60.
    pad_time = 30.*60
    event_info['orbit_time_var'] = orbit_time_var
    event_info['orbit_time_step'] = time_step
    event_info['orbit_time_range'] = time_range+[-1,1]*pad_time
    times = alfven_arc_load_rbsp_data_orbit_ut(event_info, time_var=orbit_time_var)
    r_gsm_var = alfven_arc_load_rbsp_data_r_gsm(event_info, time_var=orbit_time_var)
    model_vars = alfven_arc_load_rbsp_data_load_model_vars(event_info, time_var=orbit_time_var)


;---HOPE vars.
    hope_vars = alfven_arc_load_rbsp_data_hope_moments(event_info)
    
    
;---Field related vars.
    field_time_var = 'field_ut'
    time_step = 1d/16
    pad_time = 30.*60
    event_info['field_time_var'] = field_time_var
    event_info['field_time_step'] = time_step
    event_info['field_time_range'] = time_range+[-1,1]*pad_time
    times = alfven_arc_load_rbsp_data_field_ut(event_info, time_var=field_time_var)
    b_gsm_var = alfven_arc_load_rbsp_data_b_gsm(event_info, time_var=field_time_var)
    density_var = alfven_arc_load_rbsp_data_density(event_info)

    ; Seperate B0 and B1.
    event_info['b0_window'] = 20.*60
    bmod_var = alfven_arc_load_rbsp_data_bmod_gsm(event_info)
    var = alfven_arc_load_rbsp_data_b0_gsm(event_info, time_var=field_time_var)
    b1_gsm_var = alfven_arc_load_rbsp_data_b1_gsm(event_info)
    
    ; Load E and Edot0 convert to FAC and Calculate Poynting flux.
    e_mgse_var = alfven_arc_load_rbsp_data_e_mgse_2015_0416(event_info, time_var=field_time_var, _extra=ex)
    
    ; Convert data to FAC.
;    fac_vars = alfven_arc_load_rbsp_data_fac_vars(event_info)
    
    ; Calc fields PS and E/B ratio.
    event_info['ebr_time_range'] = time_range
;    event_info['ebr_time_range'] = time_double(['2013-05-01/07:35','2013-05-01/07:42'])
;    ebr_vars = alfven_arc_load_rbsp_data_ebr_vars(event_info, fac_vars)

    ; Calculae pflux.
    ;pf_vars = alfven_arc_load_rbsp_data_pflux(event_info, fac_vars)


;---Plasma parameters.
    event_info['plasma_param_time_range'] = time_double(['2015-04-16/08:08','2015-04-17/08:12'])
    ;param_vars = alfven_arc_load_rbsp_data_plasma_param(event_info)



    return, event_info

end

event_info = alfven_arc_load_data('2015_0416_0800')
end