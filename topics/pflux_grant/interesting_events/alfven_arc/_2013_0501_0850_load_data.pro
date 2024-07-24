function _2013_0501_0850_load_data_dmsp_r_gsm, event_info, get_name=get_name, time_var=time_var

    probe = '18'
    prefix = 'dmsp'+probe+'_'
    var = prefix+'r_gsm'
    if keyword_set(get_name) then return, var
    time_var = var+'_ut'

    data_file = event_info['data_file']
    data_dir = file_dirname(data_file)
    if ~cdf_has_var(var, filename=data_file) then begin
        file = join_path([data_dir,'DMSPf'+probe+'_data_20130501_0855.tplot'])
        if file_test(file) eq 0 then message, 'Inconsistency ...'
        
        tplot_restore, filename=file
        get_data, prefix+'lon_tclip', times, glons
        get_data, prefix+'lat_tclip', times, glats
        get_data, prefix+'alt_tclip', times, diss
        
        deg = constant('deg')
        rad = constant('rad')
        re = constant('re')
        ntime = n_elements(times)
        ndim = 3
        r_geo = fltarr(ntime,ndim)
        glons *= rad
        glats *= rad
        diss = diss/re+1
        r_geo[*,0] = diss*cos(glats)*cos(glons)
        r_geo[*,1] = diss*cos(glats)*sin(glons)
        r_geo[*,2] = diss*sin(glats)
        r_gsm = cotran(r_geo, times, 'geo2gsm')
;        r_mag = cotran(r_geo, times, 'geo2mag')
;        mlons = atan(r_mag[*,1],r_mag[*,0])*deg
;        mlats = asin(r_mag[*,2]/diss)*deg
;        mlts = mlon2mlt(mlons, times)
        store_data, var, times, r_gsm
        add_setting, var, smart=1, dictionary($
            'display_type', 'vector', $
            'unit', 'Re', $
            'short_name', 'R', $
            'coord', 'GSM', $
            'coord_labels', constant('xyz') )
        
        time_step = 1d
        common_times = make_bins(times, time_step, inner=1)
        time_range = minmax(common_times)
        interp_time, var, common_times

        cdf_save_var, time_var, value=common_times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', time_range )
        cdf_save_setting, settings, filename=data_file, varname=time_var
        
        data = get_var_data(var, limits=limits)
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif
    if check_if_update(var) then cdf_load_var, var, filename=data_file

    return, var
    
    
end


function _2013_0501_0850_load_data_pflux_setting, event_info

    key = 'pflux_setting'
    if ~event_info.haskey(key) then begin
        filter = [0.25d,900]    ; sec.
        scale_info = {s0:min(filter), s1:max(filter), dj:1d/8, ns:0d}
        event_info[key] = dictionary($
            'filter', filter, $
            'scale_info', scale_info )
    endif
    return, event_info[key]

end


function _2013_0501_0850_load_data_field_ut, event_info, time_var=var, get_name=get_name

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

function _2013_0501_0850_load_data_orbit_ut, event_info, time_var=var, get_name=get_name

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



function _2013_0501_0850_load_data_hope_moments, event_info

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

function _2013_0501_0850_load_data_density, event_info

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

function _2013_0501_0850_load_data_e_mgse, event_info, time_var=field_time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = _2013_0501_0850_load_data_field_ut(get_name=1)
    times = _2013_0501_0850_load_data_field_ut(event_info)
    time_range = event_info['field_time_range']

    coord = 'mgse'
    resolution = 'survey'
    var = rbsp_read_efield(time_range, probe=probe, coord=coord, resolution=resolution, get_name=1)
    if ~cdf_has_var(var, filename=data_file) then begin
        var = rbsp_read_efield(time_range, probe=probe, coord=coord, resolution=resolution)
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

function _2013_0501_0850_load_data_b_gsm, event_info, time_var=field_time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = _2013_0501_0850_load_data_field_ut(get_name=1)
    times = _2013_0501_0850_load_data_field_ut(event_info)
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

function _2013_0501_0850_load_data_r_gsm, event_info, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = _2013_0501_0850_load_data_orbit_ut(get_name=1)
    times = _2013_0501_0850_load_data_orbit_ut(event_info)
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

function _2013_0501_0850_load_data_load_model_setting, event_info

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


function _2013_0501_0850_load_data_load_dmsp_model_vars, event_info

    data_file = event_info['data_file']
    prefix = event_info['prefix']

    model_setting = _2013_0501_0850_load_data_load_model_setting(event_info)
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
        prefix = 'dmsp'+['18']+'_'
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
                r_var = _2013_0501_0850_load_data_dmsp_r_gsm(event_info, time_var=time_var)
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
                window = 60d
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



function _2013_0501_0850_load_data_load_model_vars, event_info, time_var=time_var

    if n_elements(time_var) eq 0 then time_var = _2013_0501_0850_load_data_orbit_ut(get_name=1)
    data_file = event_info['data_file']
    prefix = event_info['prefix']

    model_setting = _2013_0501_0850_load_data_load_model_setting(event_info)
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
                r_var = _2013_0501_0850_load_data_r_gsm(event_info, time_var=time_var)
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

function _2013_0501_0850_load_data_bmod_gsm, event_info, time_var=time_var

    if n_elements(time_var) eq 0 then time_var = _2013_0501_0850_load_data_orbit_ut(get_name=1)
    data_file = event_info['data_file']
    prefix = event_info['prefix']

    model_setting = _2013_0501_0850_load_data_load_model_setting(event_info)
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
            r_var = _2013_0501_0850_load_data_r_gsm(event_info, time_var=time_var)
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


function _2013_0501_0850_load_data_b0_gsm, event_info, time_var=time_var

    prefix = event_info['prefix']
    data_file = event_info['data_file']
    time_range = event_info['field_time_range']

    if n_elements(time_var) eq 0 then time_var = _2013_0501_0850_load_data_field_ut(get_name=1)

    var = prefix+'b0_gsm'
    if ~cdf_has_var(var, filename=data_file) then begin
        b_gsm_var = _2013_0501_0850_load_data_b_gsm(event_info, time_var=time_var)
        bmod_gsm_vars = _2013_0501_0850_load_data_bmod_gsm(event_info)
        
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

function _2013_0501_0850_load_data_b1_gsm, event_info

    prefix = event_info['prefix']
    data_file = event_info['data_file']

    b_gsm_var = _2013_0501_0850_load_data_b_gsm(event_info)
    b0_gsm_var = _2013_0501_0850_load_data_b0_gsm(event_info)
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


function _2013_0501_0850_load_data_edot0_mgse, event_info, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']
    time_range = event_info['time_range']

    if n_elements(time_var) eq 0 then time_var = _2013_0501_0850_load_data_field_ut(get_name=1)

    var = prefix+'edot0_mgse'
    if ~cdf_has_var(var, filename=data_file) then begin
        e_mgse_var = _2013_0501_0850_load_data_e_mgse(event_info, time_var=time_var)
        b0_gsm_var = _2013_0501_0850_load_data_b0_gsm(event_info)
        common_times = _2013_0501_0850_load_data_field_ut(event_info)
        e_mgse = get_var_data(e_mgse_var, at=common_times)
        b0_gsm = get_var_data(b0_gsm_var, at=common_times)
        b0_mgse = cotran(b0_gsm, common_times, 'gsm2mgse', probe=probe)
        e_mgse[*,0] = -(e_mgse[*,1]*b0_mgse[*,1]+e_mgse[*,2]*b0_mgse[*,2])/b0_mgse[*,0]

        store_data, var, common_times, e_mgse
        add_setting, var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'Edot0', $
            'unit', 'mV/m', $
            'coord', 'mGSE', $
            'coord_labels', constant('xyz') )
            

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


function _2013_0501_0850_load_data_ebr_vars, event_info, fac_vars

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']
    time_range = event_info['ebr_time_range']

    pflux_setting = _2013_0501_0850_load_data_pflux_setting(event_info)
    scale_info = pflux_setting['scale_info']

    e_index = 2
    b_index = 1
    fac_labels = ['b','w','o']
    vars = list()
    foreach type, ['','dot0'] do begin
        ; |E| and |B1|
        e_var = prefix+'e'+type+'_fac'
        b_var = prefix+'b1_fac'
        ebr_var = prefix+'ebr_total'
        vars.add, stplot_calc_ebratio(time_range, $
            e_var=e_var, b_var=b_var, $
            scale_info=scale_info, output=ebr_var)
        vars.add, [e_var,b_var]+'_mor', extract=1
        
        ; E_out and B_west.
        e_var = stplot_index(e_var, e_index, output=prefix+'e'+type+fac_labels[e_index])
        b_var = stplot_index(b_var, b_index, output=prefix+'b'+fac_labels[b_index])
        ebr_var = prefix+'ebr_component'
        vars.add, stplot_calc_ebratio(time_range, $
            e_var=e_var, b_var=b_var, $
            scale_info=scale_info, output=ebr_var)
        vars.add, [e_var,b_var]+'_mor', extract=1
    endforeach

    vars = vars.toarray()
    return, vars

end


function _2013_0501_0850_load_data_fac_vars, event_info

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']
    time_range = event_info['time_range']

    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    event_info['fac_labels'] = fac_labels

    q_gsm2fac_var = prefix+'q_gsm2fac'
    b1_fac_var = prefix+'b1_fac'
    e_fac_var = prefix+'e_fac'
    edot0_fac_var = prefix+'edot0_fac'
    fac_vars = [b1_fac_var,e_fac_var,edot0_fac_var]
    save_vars = [q_gsm2fac_var,fac_vars]

    load_data = 0
    foreach var, save_vars do begin
        if cdf_has_var(var, filename=data_file) then continue
        load_data = 1
        break
    endforeach

    if load_data then begin
        time_var = _2013_0501_0850_load_data_field_ut(get_name=1)
        b0_gsm_var = _2013_0501_0850_load_data_b0_gsm(event_info)
        r_gsm_var = _2013_0501_0850_load_data_r_gsm(event_info)
        define_fac, b0_gsm_var, r_gsm_var, time_var=b0_gsm_var

        b1_gsm_var = _2013_0501_0850_load_data_b1_gsm(event_info)
        e_mgse_var = _2013_0501_0850_load_data_e_mgse(event_info)
        edot0_mgse_var = _2013_0501_0850_load_data_edot0_mgse(event_info)

        mgse_vars = [e_mgse_var,edot0_mgse_var]
        e_gsm_var = prefix+'e_gsm'
        edot0_gsm_var = prefix+'edot0_gsm'
        gsm_vars = [e_gsm_var,edot0_gsm_var]
        foreach mgse_var, mgse_vars, var_id do begin
            get_data, mgse_var, times, vec_mgse
            vec_gsm = cotran(vec_mgse, times, 'mgse2gsm', probe=probe)
            gsm_var = gsm_vars[var_id]
            store_data, gsm_var, times, vec_gsm
            add_setting, gsm_var, smart=1, dictionary($
                'display_type', 'vector', $
                'short_name', get_setting(mgse_var,'short_name'), $
                'unit', 'mV/m', $
                'coord', 'GSM', $
                'coord_labels', constant('xyz') )
        endforeach

        gsm_vars = [b1_gsm_var,e_gsm_var,edot0_gsm_var]
        foreach var, gsm_vars, var_id do begin
            to_fac, var, to=fac_vars[var_id]
            add_setting, fac_vars[var_id], smart=1, dictionary($
                'display_type', 'vector', $
                'short_name', get_setting(var,'short_name'), $
                'unit', get_setting(var,'unit'), $
                'coord', 'FAC', $
                'coord_labels', fac_labels)
        endforeach


        foreach var, save_vars do begin
            if tnames(var) eq '' then stop
            data = get_var_data(var, limits=limits)
            cdf_save_var, var, value=data, filename=data_file
            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=var
        endforeach
    endif

    
    foreach var, save_vars do begin
        if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    endforeach
    return, save_vars

end



function _2013_0501_0850_load_data_pflux, event_info, fac_vars

    prefix = event_info['prefix']

    pflux_setting = _2013_0501_0850_load_data_pflux_setting(event_info)
    scale_info = pflux_setting['scale_info']
    if n_elements(fac_vars) eq 0 then begin
        fac_vars = _2013_0501_0850_load_data_fac_vars(event_info)
    endif
    
    pf_fac_vars = []
    pf_spec_vars = []
    spec_unit = tex2str('mu')+'W/m!E2!N'
    spec_ct = 66
    spec_zrange = [-1,1]*10
    spec_zstep = 10
    spec_zminor = 5
    spec_ztickv = make_bins(spec_zrange, spec_zstep, inner=1)
    spec_zticks = n_elements(spec_ztickv)-1
    spec_zticklen = -0.3
    types = ['','dot0']
    
    model_setting = event_info['model_setting']
    model = model_setting['external_model']
    internal_model = model_setting['internal_model']
    bf_var = prefix+'bf_gsm_'+model+'_'+internal_model+'_north'
    b0_var = prefix+'b0_gsm'
    b0_gsm = get_var_data(b0_var, times=times)
    bf_gsm = get_var_data(bf_var, at=times)
    cmap = snorm(bf_gsm)/snorm(b0_gsm)
    ndim = 3
    
    foreach type, types do begin
        b1_fac_var = prefix+'b1_fac'
        e1_fac_var = prefix+'e'+type+'_fac'
        pf_fac_var = prefix+'pf'+type+'_fac'
        stplot_calc_pflux_mor, e1_fac_var, b1_fac_var, pf_fac_var, scaleinfo=scale_info
        add_setting, pf_fac_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'S', $
            'unit', 'mW/m!E2!N', $
            'coord', 'FAC', $
            'coord_labels', event_info['fac_labels'] )
        pf_fac_vars = [pf_fac_vars,pf_fac_var]

        ; Normalize to 100 km.
        pf_fac_map_var = pf_fac_var+'_map'
        pf_fac = get_var_data(pf_fac_var, times=times)
        pf_fac_map = pf_fac
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

        ; Spec.
        pf_spec_var = prefix+'pf'+type+'_fac_mor_spec_1'
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
        pf_spec_vars = [pf_spec_vars,pf_spec_var]
    endforeach
    
    ; Calculate ExB velocity.
    b0_gsm_var = prefix+'b0_gsm'
    b0_gsm = get_var_data(b0_gsm_var, times=times)
    b0_fac = b0_gsm
    b0_fac[*,0] = snorm(b0_gsm)
    b0_fac[*,1:2] = 0
    coef_fac = b0_fac
    coef_fac[*,0] /= (b0_fac[*,0])^2
    foreach type, ['','dot0'] do begin
        e_fac_var = prefix+'e'+type+'_fac'
        e_fac = get_var_data(e_fac_var, at=times)
        vexb_fac = vec_cross(e_fac,coef_fac)*1e3
        vexb_fac_var = prefix+'vexb'+type+'_fac'
        store_data, vexb_fac_var, times, vexb_fac
        add_setting, vexb_fac_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'ExB V', $
            'unit', 'km/s', $
            'coord', 'FAC', $
            'coord_labels', event_info['fac_labels'] )
    endforeach
    
    return, pf_spec_vars
    
end


function _2013_0501_0850_load_data_asi_ut, event_info, time_var=var, get_name=get_name

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


function _2013_0501_0850_load_data_asi, event_info, filename=data_file, time_var=time_var


    asi_setting = dictionary($
        'sites', ['fsmi','gako','atha','tpas','fsim'], $
        'min_elev', [0.5,0.5,0.5,0.5,10], $
        'best_site', 'atha', $
        'mlt_range', [-4d,2], $
        'mlat_range', [55d,70] )
    event_info['asi_setting'] = asi_setting
    prefix = event_info['prefix']
    data_file = event_info['data_file']
    time_range = event_info['asi_time_range']
    if n_elements(time_var) eq 0 then time_var = _2013_0501_0850_load_data_asi_ut(event_info, get_name=1)
    
    min_elev = asi_setting['min_elev']
    merge_method = 'merge_elev'
    calibration_method = 'simple'

    mlt_image_var = themis_asf_read_mlt_image(time_range, get_name=1)
    if ~cdf_has_var(mlt_image_var, filename=data_file) then begin
        sites = asi_setting['sites']
        mlt_image_var = themis_asf_read_mlt_image(time_range, sites=sites, min_elev=min_elev, merge_method=merge_method, calibration_method=calibration_method)
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
            min_elev=min_elev, merge_method=merge_method, calibration_method=calibration_method, $
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

;---Ewogram.
    

;---Keogram.
    fmlt_vars = prefix+'fmlt_'+['dipole','igrf']
;    mlt_range = []
;    foreach var, fmlt_vars do begin
;        mlt_range = [mlt_range,minmax(get_var_data(var,in=time_range))]
;    endforeach
;    mlt_range = minmax(mlt_range)
    mlt_range = []
    model_setting = _2013_0501_0850_load_data_load_model_setting(event_info)
    internal_model = model_setting['internal_model']
    mlt_var = prefix+'fmlt_'+internal_model+'_north'
    models = model_setting['models']
    external_model = model_setting['external_model']
    model_index = where(models eq external_model)
    mlt_range = reform((get_var_data(mlt_var,at=time_range))[*,model_index])
    mlt_range = round(mlt_range*10)*0.1
    mlat_range = [60,67]
    keo_var = themis_read_mlt_image_rect_keo(mlt_image_var=mlt_image_rect_var, $
        mlt_range=mlt_range, mlat_range=mlat_range)
    zlim, keo_var, 1,8e3, 0

    return, [mlt_image_var,mlt_image_rect_var,keo_var]

end



function _2013_0501_0850_load_data_plasma_param, event_info

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
    
    vf_exb_fac = fltarr(ndim)
    vexb_fac = get_var_data(prefix+'vexbdot0_fac', in=the_time_range)
    for ii=0,ndim-1 do vf_fac[ii] = mean(abs(vexb_fac[*,ii]),nan=1)
    plasma_param['vf_exb_fac'] = vf_fac
    plasma_param['vf_exb'] = vf_fac[1]

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

    the_time = event_info['snapshot_time']
    r_gsm = get_var_data(prefix+'r_gsm', at=the_time)
    plasma_param['r_gsm'] = r_gsm
    plasma_param['r_sm'] = cotran(r_gsm, the_time, 'gsm2sm')

    event_info['plasma_param'] = plasma_param
    return, plasma_param

end


function _2013_0501_0850_load_data_weygand_ut, event_info, time_var=var, get_name=get_name

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


function _2013_0501_0850_load_data_weygand, event_info, time_var=time_var

    data_file = event_info['data_file']
    time_range = event_info['weygand_time_range']
    if n_elements(time_var) eq 0 then time_var = _2013_0501_0850_load_data_weygand_ut(get_name=1)

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


function _2013_0501_0850_load_data, filename=data_file

    root_dir = join_path([googledir(),'works','pflux_grant'])
    root_dir = join_path([diskdir('Research'),'archives','works','2024_pflux_grant'])
    if n_elements(data_file) eq 0 then begin
        data_file = join_path([root_dir,'alfven_arc','data','2013_0501_0850_all_data_v02.cdf'])
    endif
    

    if file_test(data_file) eq 0 then begin
        probe = 'b'
        time_range = time_double(['2013-05-01/07:30','2013-05-01/09:10'])
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
    
    snapshot_time = time_double('2013-05-01/08:54:30')
    event_info['snapshot_time'] = snapshot_time
    
    
    

;---DMSP data.
    var = _2013_0501_0850_load_data_dmsp_r_gsm(event_info)
    var = _2013_0501_0850_load_data_load_dmsp_model_vars(event_info)
    r_var = 'dmsp18_r_gsm'
    vinfo = geopack_read_bfield(r_var, models='t89', igrf=0, suffix='_dipole', t89_par=2, coord='gsm')
    get_data, 'dmsp18_bf_gsm_t89_dipole_south', times, b1
    get_data, 'dmsp18_bmod_gsm_t89_dipole', times, b0
    cmap = smooth(snorm(b1)/snorm(b0),100, edge_mirror=1)
    store_data, 'dmsp18_cmap', times, cmap
    add_setting, 'dmsp18_cmap', smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'Cmap', $
        'unit', '#', $
        'text', 'Mapping coef to map flux to 100 km altitude' )
    
    
    
;---SECS.
    weygand_time_var = 'weygand_ut'
    weygand_time_step = 10d
    event_info['weygand_time_var'] = weygand_time_var
    event_info['weygand_time_step'] = weygand_time_step
    event_info['weygand_time_range'] = time_double(['2013-05-01/07:00','2013-05-01/09:30'])
    times = _2013_0501_0850_load_data_weygand_ut(event_info, time_var=weygand_time_var)
    j_vars = _2013_0501_0850_load_data_weygand(event_info, time_var=weygand_time_var)
    
    
;---HOPE vars.
    hope_vars = _2013_0501_0850_load_data_hope_moments(event_info)


;---Orbit related vars.
    orbit_time_var = 'orbit_ut'
    time_step = 60.
    pad_time = 30.*60
    event_info['orbit_time_var'] = orbit_time_var
    event_info['orbit_time_step'] = time_step
    event_info['orbit_time_range'] = time_range+[-1,1]*pad_time
    times = _2013_0501_0850_load_data_orbit_ut(event_info, time_var=orbit_time_var)
    r_gsm_var = _2013_0501_0850_load_data_r_gsm(event_info, time_var=orbit_time_var)
    model_vars = _2013_0501_0850_load_data_load_model_vars(event_info, time_var=orbit_time_var)

;---Field related vars.
    field_time_var = 'field_ut'
    time_step = 1d/16
    pad_time = 30.*60
    event_info['field_time_var'] = field_time_var
    event_info['field_time_step'] = time_step
    event_info['field_time_range'] = time_range+[-1,1]*pad_time
    times = _2013_0501_0850_load_data_field_ut(event_info, time_var=field_time_var)
    b_gsm_var = _2013_0501_0850_load_data_b_gsm(event_info, time_var=field_time_var)
    e_mgse_var = _2013_0501_0850_load_data_e_mgse(event_info, time_var=field_time_var)
    density_var = _2013_0501_0850_load_data_density(event_info)

    ; Seperate B0 and B1.
    event_info['b0_window'] = 15.*60
    bmod_var = _2013_0501_0850_load_data_bmod_gsm(event_info)
    var = _2013_0501_0850_load_data_b0_gsm(event_info, time_var=field_time_var)
    b1_gsm_var = _2013_0501_0850_load_data_b1_gsm(event_info)
    
    ; Calculate Edot0.
    edot0_mgse_var = _2013_0501_0850_load_data_edot0_mgse(event_info, time_var=field_time_var)
    
    ; Convert data to FAC.
    fac_vars = _2013_0501_0850_load_data_fac_vars(event_info)
    
    ; Calc fields PS and E/B ratio.
    event_info['ebr_time_range'] = time_range
;    event_info['ebr_time_range'] = time_double(['2013-05-01/07:35','2013-05-01/07:42'])
    ebr_vars = _2013_0501_0850_load_data_ebr_vars(event_info, fac_vars)

    ; Calculae pflux.
    pf_vars = _2013_0501_0850_load_data_pflux(event_info, fac_vars)



;---ASI vars.
    asi_time_var = 'asi_ut'
    asi_time_step = 3d
    event_info['asi_time_var'] = asi_time_var
    event_info['asi_time_step'] = asi_time_step
    event_info['asi_time_range'] = time_double(['2013-05-01/07:00','2013-05-01/09:30'])
    times = _2013_0501_0850_load_data_asi_ut(event_info, time_var=asi_time_var)
    asi_vars = _2013_0501_0850_load_data_asi(event_info, time_var=asi_time_var)


;---Plasma parameters.
    event_info['plasma_param_time_range'] = time_double(['2013-05-01/08:40','2013-05-01/08:50'])
;    param_vars = _2013_0501_0850_load_data_plasma_param(event_info)


    plot_dir = join_path([root_dir,'alfven_arc','plot'])
    event_info['plot_dir'] = plot_dir


    return, event_info

end



pinfo = _2013_0501_0850_load_data()
end