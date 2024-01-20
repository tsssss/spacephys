
function event_study_geopack_read_model_vars, event_info, time_var=time_var, r_gsm_var=orbit_var, model_setting=model_setting, get_name=get_name, update=update

    mission = event_info['mission']
    probe = event_info['probe']
    prefix = event_info['prefix']
    data_file = event_info['data_file']
    

    if n_elements(model_setting) eq 0 then model_setting = event_study_get_model_setting(event_info)
    models = model_setting['models']

    all_model_vars = []
    foreach direction, ['south','north'] do begin
        basic_vars = prefix+['fmlt','fmlon','fmlat','bf_gsm_'+models]
        foreach internal_model, ['igrf','dipole'] do begin
            load_data = 0
            key = internal_model+'_'+direction
            suffix = '_'+key
            model_vars = basic_vars+suffix
            all_model_vars = [all_model_vars,model_vars]
            if keyword_set(get_name) then continue

            foreach var, model_vars do begin
                if ~cdf_has_var(var, filename=data_file) then begin
                    load_data = 1
                    break
                endif
            endforeach

            if keyword_set(update) then load_data = 1
            if load_data then begin
                if direction eq 'north' then begin
                    north = 1
                    south = 0
                endif else begin
                    north = 0
                    south = 1
                endelse
                refine = model_setting['refine']

                orbit_var = event_study_read_orbit(event_info, time_var=time_var, coord='gsm')
                time_range = get_setting(orbit_var, 'requested_time_range')
                igrf = (internal_model eq 'igrf')? 1: 0
                vinfo = geopack_trace_to_ionosphere(orbit_var, models=models, $
                    igrf=igrf, south=south, north=north, refine=refine, suffix=suffix)
                foreach var, model_vars do begin
                    if tnames(var) eq '' then stop
                    data = get_var_data(var, limits=limits)
                    cdf_save_var, var, value=data, filename=data_file
                    settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
                    settings['depend_0'] = time_var
                    settings['var_type'] = 'data'
                    settings['requested_time_range'] = time_range
                    cdf_save_setting, settings, filename=data_file, varname=var
                endforeach
            endif


            foreach var, model_vars do begin
                if keyword_set(update) then del_data, var
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
        endforeach
    endforeach

    return, all_model_vars

end