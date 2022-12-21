function _2013_0501_load_data_hope_moments, event_info

    species = ['e','p','o']
    event_info['species'] = species
    data_file = event_info['data_file']
    probe = event_info['probe']
    prefix = event_info['prefix']
    time_range = event_info['time_range']
    secofday = constant('secofday')
    time_range = time_range[0]-(time_range[0] mod secofday)+[0,secofday]

    hope_vars = []
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
        stop
    endif
end

function _2013_0501_load_data_density, event_info

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']
    time_range = event_info['time_range']
    secofday = constant('secofday')
    time_range = time_range[0]-(time_range[0] mod secofday)+[0,secofday]

    density_emfisis_var = rbsp_read_density_emfisis(time_range, probe=probe)
    density_efw_var = rbsp_read_density_efw(time_range, probe=probe)
    density_hope_var = rbsp_read_density_hope(time_range, probe=probe)

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


    foreach var, vars do begin
        if cdf_has_var(var, filename=data_file) then begin
            cdf_load_var, var, filename=data_file
        endif else begin
            if n_elements(times) eq 0 then times = cdf_read_var(time_var, filename=data_file)
            interp_time, var, times
            
            cdf_save_var, var, value=get_var_data(var, limits=limits), filename=data_file

            settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
            settings['depend_0'] = time_var
            settings['var_type'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=var
        endelse
    endforeach

    var = prefix+'density_combo'
    var = stplot_merge(vars, output=var)
    add_setting, var, smart=1, dictionary($
        'ylog', 1, $
        'display_type', 'stack', $
        'unit', 'cm!U-3!N', $
        'labels', 'N '+labels, $
        'colors', colors )
    return, var

end

function _2013_0501_load_data_e_mgse, event_info, time_var=field_time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = _2013_0501_load_data_field_ut(get_name=1)
    times = _2013_0501_load_data_field_ut(event_info)
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

function _2013_0501_load_data_b_gsm, event_info, time_var=field_time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = _2013_0501_load_data_field_ut(get_name=1)
    times = _2013_0501_load_data_field_ut(event_info)
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


function _2013_0501_load_data_field_ut, event_info, time_var=var, get_name=get_name

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


function _2013_0501_load_data_load_model_vars, event_info, time_var=time_var

    if n_elements(time_var) eq 0 then time_var = _2013_0501_load_data_orbit_ut(get_name=1)
    data_file = event_info['data_file']

    model_setting = dictionary($
        'model', 't01', $
        'igrf', 0, $
        'refine', 1, $
        'direction', -1 )
    models = ['t89','t96','t01','t04s']
    event_info['model_setting'] = model_setting
    event_info['models'] = models
    
    prefix = event_info['prefix']
    vars = prefix+['fmlt','fmlon','fmlat','bf_gsm_'+models]
    direction = model_setting['direction']

    model_vars = []
    foreach internal_model, ['igrf','dipole'] do begin
        load_data = 0
        the_vars = vars+'_'+internal_model
        model_vars = [model_vars,the_vars]
        foreach var, the_vars do begin
            if ~cdf_has_var(var, filename=data_file) then begin
                load_data = 1
                break
            endif
        endforeach
        if load_data then begin
            r_var = _2013_0501_load_data_r_gsm(event_info, time_var=time_var)
            igrf = (internal_model eq 'igrf')? 1: 0
            vinfo = geopack_trace_to_ionosphere(r_var, models=models, $
                igrf=igrf, direction=direction, suffix='_'+internal_model)
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
        if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file
    endforeach

    return, model_vars

end


function _2013_0501_load_data_r_gsm, event_info, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = _2013_0501_load_data_orbit_ut(get_name=1)
    times = _2013_0501_load_data_orbit_ut(event_info)
    time_range = event_info['orbit_time_range']

    coord = 'gsm'
    var = rbsp_read_orbit(time_range, probe=probe, coord=coord, get_name=1)
    if ~cdf_has_var(var, filename=data_file) then begin
        var = rbsp_read_orbit(time_range, probe=probe, coord=coord)
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


function _2013_0501_load_data_orbit_ut, event_info, time_var=var, get_name=get_name

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


function _2013_0501_load_data, filename=data_file

    if n_elements(data_file) eq 0 then begin
        data_file = join_path([googledir(),'works','pflux_grant','2013_0501_arc','data','2013_0501_all_data_v02.cdf'])
    endif

    if file_test(data_file) eq 0 then begin
        probe = 'b'
        time_range = time_double(['2013-05-01/07:25','2013-05-01/07:55'])
        cdf_save_setting, 'probe', probe, filename=data_file
        cdf_save_setting, 'prefix', 'rbsp'+probe+'_', filename=data_file
        cdf_save_setting, 'time_range', time_range, filename=data_file
    endif
    event_info = cdf_read_setting(filename=data_file)
    event_info['data_file'] = data_file

    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']

;---Orbit related vars.
    orbit_time_var = 'orbit_ut'
    time_step = 60.
    pad_time = 30.*60
    event_info['orbit_time_var'] = orbit_time_var
    event_info['orbit_time_step'] = time_step
    event_info['orbit_time_range'] = time_range+[-1,1]*pad_time
    times = _2013_0501_load_data_orbit_ut(event_info, time_var=orbit_time_var)
    r_gsm_var = _2013_0501_load_data_r_gsm(event_info, time_var=orbit_time_var)
    model_vars = _2013_0501_load_data_load_model_vars(event_info, time_var=orbit_time_var)

;---Field related vars.
    field_time_var = 'field_ut'
    time_step = 1d/16
    pad_time = 30.*60
    event_info['field_time_var'] = field_time_var
    event_info['field_time_step'] = time_step
    event_info['field_time_range'] = time_range+[-1,1]*pad_time
    times = _2013_0501_load_data_field_ut(event_info, time_var=field_time_var)
    b_gsm_var = _2013_0501_load_data_b_gsm(event_info, time_var=field_time_var)
    e_mgse_var = _2013_0501_load_data_e_mgse(event_info, time_var=field_time_var)
    density_var = _2013_0501_load_data_density(event_info)
    
;---HOPE vars.
    hope_vars = _2013_0501_load_data_hope_moments(event_info)
    


end




tmp = _2013_0501_load_data()
end
