function alfven_arc_load_themis_data_orbit_ut, event_info, time_var=var, get_name=get_name

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


function alfven_arc_load_themis_data_r_gsm, event_info, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_themis_data_orbit_ut(get_name=1)
    times = alfven_arc_load_themis_data_orbit_ut(event_info)
    time_range = event_info['orbit_time_range']

    coord = 'gsm'
    var = themis_read_orbit(time_range, probe=probe, coord=coord, get_name=1)
    
    if ~cdf_has_var(var, filename=data_file) then begin
        var = themis_read_orbit(time_range+[-1,1]*60, probe=probe, coord=coord)
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


function alfven_arc_load_themis_data_load_model_setting, event_info

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



function alfven_arc_load_themis_data_load_model_vars, event_info, time_var=time_var

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_themis_data_orbit_ut(get_name=1)
    data_file = event_info['data_file']
    prefix = event_info['prefix']

    model_setting = alfven_arc_load_themis_data_load_model_setting(event_info)
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
                r_var = alfven_arc_load_themis_data_r_gsm(event_info, time_var=time_var)
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



function alfven_arc_load_themis_data_field_ut, event_info, time_var=var, get_name=get_name

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



function alfven_arc_load_themis_data_b_gsm, event_info, time_var=field_time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_themis_data_field_ut(get_name=1)
    times = alfven_arc_load_themis_data_field_ut(event_info)
    time_range = event_info['field_time_range']

    coord = 'gsm'
    resolution = 'fgl'
    var = themis_read_bfield(time_range, probe=probe, coord=coord, id=resolution, get_name=1)
    if ~cdf_has_var(var, filename=data_file) then begin
        var = themis_read_bfield(time_range, probe=probe, coord=coord, id=resolution)
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



function alfven_arc_load_themis_data_e_dsl, event_info, time_var=field_time_var, bad_e_time_ranges=bad_e_time_ranges

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_themis_data_field_ut(get_name=1)
    times = alfven_arc_load_themis_data_field_ut(event_info)
    time_range = event_info['field_time_range']

    coord = 'themis_dsl'
    resolution = 'survey'
    if n_elements(bad_e_time_ranges) eq 0 then bad_e_time_ranges = list()
    var = themis_read_efield(time_range, probe=probe, coord=coord, id=resolution, get_name=1)
    if ~cdf_has_var(var, filename=data_file) then begin
        var = themis_read_efield(time_range, probe=probe, coord=coord, id=resolution)
        interp_time, var, times
        data = get_var_data(var, limits=limits)
        foreach tr, bad_e_time_ranges do begin
            index = where_pro(times,'[]',time_double(tr), count=count)
            if count eq 0 then continue
            data[index,*] = !values.f_nan
        endforeach
        store_data, var, times, data
        if n_elements(bad_e_time_ranges) eq 0 then bad_time_ranges = [] else bad_time_ranges = bad_e_time_ranges.toarray()
        options, var, 'bad_time_ranges', bad_time_ranges
        cdf_save_var, var, value=data, filename=data_file
        settings = (isa(limits,'struct'))? dictionary(limits): dictionary()
        settings['depend_0'] = time_var
        settings['var_type'] = 'data'
        cdf_save_setting, settings, filename=data_file, varname=var
    endif

    if check_if_update(var, time_range) then cdf_load_var, var, filename=data_file

    return, var
end


function alfven_arc_load_data_load_model_setting, event_info

    key = 'model_setting'
    if ~event_info.haskey(key) then begin
        model_setting = dictionary($
            'external_model', 't89', $
            'models', ['t89','t96','t01','t04s'], $
            't89_par', 2, $
            'igrf', 0, $
            'refine', 1, $
            'direction', 'south' )
        internal_model = (model_setting['igrf'] eq 0)? 'dipole': 'igrf'
        model_setting['internal_model'] = internal_model
        event_info[key] = model_setting
    endif
    return, event_info[key]

end



function alfven_arc_load_themis_data_bmod_gsm, event_info, time_var=time_var

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_themis_data_orbit_ut(get_name=1)
    data_file = event_info['data_file']
    prefix = event_info['prefix']

    model_setting = alfven_arc_load_data_load_model_setting(event_info)
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
            r_var = alfven_arc_load_themis_data_r_gsm(event_info, time_var=time_var)
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



function alfven_arc_load_themis_data_b0_gsm, event_info, time_var=time_var

    prefix = event_info['prefix']
    data_file = event_info['data_file']
    time_range = event_info['field_time_range']

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_themis_data_field_ut(get_name=1)

    var = prefix+'b0_gsm'
    window = event_info['b0_window']
    if cdf_has_var(var, filename=data_file) then begin
        vatt = cdf_read_setting(var, filename=data_file)
        if window ne vatt['window'] then cdf_del_var, var, filename=data_file
    endif
    if ~cdf_has_var(var, filename=data_file) then begin
        b_gsm_var = alfven_arc_load_themis_data_b_gsm(event_info, time_var=time_var)
        bmod_gsm_vars = alfven_arc_load_themis_data_bmod_gsm(event_info)
        
        igrf = 0
        internal_model = (igrf eq 1)? 'igrf': 'dipole'
        external_model = 't89'
        suffix = '_'+internal_model
        bmod_gsm_var = prefix+'bmod_gsm_'+external_model+suffix

        b_gsm = get_var_data(b_gsm_var, times=times)
        bmod_gsm = get_var_data(bmod_gsm_var, at=times)
        b1_gsm = b_gsm-bmod_gsm
        ndim = 3
        time_step = total(times[0:1]*[-1,1])
        width = window/time_step
        for ii=0,ndim-1 do begin
            b1_gsm[*,ii] -= smooth(b1_gsm[*,ii], width, edge_mirror=1, nan=1)
        endfor
        
        b0_gsm = b_gsm-b1_gsm
        store_data, var, times, b0_gsm
        add_setting, var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', strupcase(external_model)+' B', $
            'unit', 'nT', $
            'coord', 'GSM', $
            'coord_labels', constant('xyz'), $
            'window', window, $
            'external_model', external_model, $
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

function alfven_arc_load_themis_data_b1_gsm, event_info

    prefix = event_info['prefix']
    data_file = event_info['data_file']

    ; Themis B background too small, so we use B-B_mod as B1.
    b_gsm_var = alfven_arc_load_themis_data_b_gsm(event_info)
;    b0_gsm_var = alfven_arc_load_themis_data_b0_gsm(event_info)
;    b_gsm = get_var_data(b_gsm_var, times=times)
;    b0_gsm = get_var_data(b0_gsm_var, at=times)
;    b1_gsm = b_gsm-b0_gsm
    igrf = 0
    internal_model = (igrf eq 1)? 'igrf': 'dipole'
    external_model = 't89'
    suffix = '_'+internal_model
    bmod_gsm_var = prefix+'bmod_gsm_'+external_model+suffix
    b_gsm = get_var_data(b_gsm_var, times=times)
    bmod_gsm = get_var_data(bmod_gsm_var, at=times)
    b1_gsm = b_gsm-bmod_gsm

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


function alfven_arc_load_themis_data_edot0_dsl, event_info, time_var=time_var

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']
    time_range = event_info['time_range']

    if n_elements(time_var) eq 0 then time_var = alfven_arc_load_themis_data_field_ut(get_name=1)
    var = prefix+'edot0_dsl'
    
    if ~cdf_has_var(var, filename=data_file) then begin
        e_dsl_var = alfven_arc_load_themis_data_e_dsl(event_info, time_var=time_var)
        b0_gsm_var = alfven_arc_load_themis_data_b0_gsm(event_info)
        common_times = alfven_arc_load_themis_data_field_ut(event_info)
        interp_time, e_dsl_var, common_times   ; deal with nan.
        e_dsl = get_var_data(e_dsl_var, limits=lim)
        b0_gsm = get_var_data(b0_gsm_var, at=common_times)
        b0_dsl = cotran_pro(b0_gsm, common_times, 'gsm2themis_dsl', probe=probe)
        e_dsl[*,2] = -(e_dsl[*,0]*b0_dsl[*,0]+e_dsl[*,1]*b0_dsl[*,1])/b0_dsl[*,2]

        store_data, var, common_times, e_dsl
        add_setting, var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'Edot0', $
            'unit', 'mV/m', $
            'coord', 'THEMIS_DSL', $
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



function alfven_arc_load_themis_data_fac_vars, event_info

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
        time_var = alfven_arc_load_themis_data_field_ut(get_name=1)
        b0_gsm_var = alfven_arc_load_themis_data_b0_gsm(event_info)
        r_gsm_var = alfven_arc_load_themis_data_r_gsm(event_info)
        define_fac, b0_gsm_var, r_gsm_var, time_var=b0_gsm_var

        b1_gsm_var = alfven_arc_load_themis_data_b1_gsm(event_info)
        e_dsl_var = alfven_arc_load_themis_data_e_dsl(event_info)
        edot0_dsl_var = alfven_arc_load_themis_data_edot0_dsl(event_info)

        dsl_vars = [e_dsl_var,edot0_dsl_var]
        e_gsm_var = prefix+'e_gsm'
        edot0_gsm_var = prefix+'edot0_gsm'
        gsm_vars = [e_gsm_var,edot0_gsm_var]
        foreach dsl_var, dsl_vars, var_id do begin
            get_data, dsl_var, times, vec_dsl
            vec_gsm = cotran_pro(vec_dsl, times, 'themis_dsl2gsm', probe=probe)
            gsm_var = gsm_vars[var_id]
            store_data, gsm_var, times, vec_gsm
            add_setting, gsm_var, smart=1, dictionary($
                'display_type', 'vector', $
                'short_name', get_setting(dsl_var,'short_name'), $
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


function alfven_arc_load_themis_data_pflux_setting, event_info

    key = 'pflux_setting'
    if ~event_info.haskey(key) then begin
        filter = [0.25d,600]    ; sec.
;        filter = [10d,1800]    ; sec.
        scale_info = {s0:min(filter), s1:max(filter), dj:1d/8, ns:0d}
        event_info[key] = dictionary($
            'filter', filter, $
            'scale_info', scale_info )
    endif
    return, event_info[key]

end


function alfven_arc_load_themis_data_pflux, event_info, fac_vars

    prefix = event_info['prefix']

    pflux_setting = alfven_arc_load_themis_data_pflux_setting(event_info)
    scale_info = pflux_setting['scale_info']
    if n_elements(fac_vars) eq 0 then begin
        fac_vars = alfven_arc_load_themis_data_fac_vars(event_info)
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
    external_model = model_setting['external_model']
    internal_model = model_setting['internal_model']
    bf_var = prefix+'bf_gsm_'+external_model+'_'+internal_model+'_north'
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



function alfven_arc_load_themis_data_particle_data, event_info

    species = ['e','p']
    event_info['species'] = species
    data_file = event_info['data_file']
    probe = event_info['probe']
    prefix = event_info['prefix']
    time_range = event_info['time_range']
    secofday = constant('secofday')
    time_range = time_range+[-1,1]*3600*2

    particle_vars = []
    coord = 'gsm'
;    moments = ['en_spec'+['','_'+['anti','perp','para']],'n','t',$
;        'density','temp',['vbulk','nflux','eflux','enthalpy']+'_'+coord]
    moments = ['en_spec'+['','_'+['anti','perp','para']], $
        'n','tavg','pavg',['vbulk','nflux','eflux','keflux','enthalpy','hflux']+'_'+coord]
    foreach species_name, species do begin
        particle_vars = [particle_vars,prefix+species_name+'_'+moments]
    endforeach

    load_data = 0
    foreach var, particle_vars do begin
        if ~cdf_has_var(var, filename=data_file) then begin
            load_data = 1
            break
        endif
    endforeach

    if load_data then begin
        foreach species_name, species do begin
            en_spec_combo = themis_read_en_spec_combo(time_range, probe=probe, species=species_name)
            en_spec_var = themis_read_en_spec(time_range, probe=probe, species=species_name)

        ;---Energy spec vars.
            vars = [en_spec_var,(en_spec_combo.values()).toarray()]
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
            mom_combo = themis_read_moment_combo(time_range, probe=probe, species=species_name, coord=coord)
            mom_vars = mom_combo.values()
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

        endforeach
        
    endif

    foreach var, particle_vars do begin
        if check_if_update(var, time_range, dtime=60) then cdf_load_var, var, filename=data_file
    endforeach

    return, particle_vars
end




function alfven_arc_load_themis_data_density, event_info

    prefix = event_info['prefix']
    probe = event_info['probe']
    data_file = event_info['data_file']
    time_range = event_info['time_range']
    secofday = constant('secofday')
    time_range = time_range[0]-(time_range[0] mod secofday)+[0,secofday]

    density_esa_var = themis_read_density_esa(time_range, probe=probe, get_name=1)
    density_efi_var = themis_read_density_efi(time_range, probe=probe, get_name=1)

    vars = [density_esa_var,density_efi_var]
    labels = ['ESA','EFI']
    colors = sgcolor(['red','green'])

    time_var = 'density_ut'
    if ~cdf_has_var(time_var, filename=data_file) then begin
        time_step = 3
        times = make_bins(time_range, time_step)
        
        cdf_save_var, time_var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', time_range )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif


    foreach var, vars, var_id do begin
        if ~cdf_has_var(var, filename=data_file) then begin
            routine = 'themis_read_density_'+strlowcase(labels[var_id])
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



function alfven_arc_load_themis_data, input_time_range, probe=probe, filename=data_file, _extra=ex

    time_range = time_double(input_time_range)
    if n_elements(probe) eq 0 then message, 'No input probe ...'
    

    if n_elements(data_file) eq 0 then begin
        base = time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_themis'+probe+'_data_v01.cdf'
        data_file = join_path([googledir(),'works','pflux_grant','alfven_arc','data',base])
    endif
    

    if file_test(data_file) eq 0 then begin
        cdf_save_setting, 'probe', probe, filename=data_file
        cdf_save_setting, 'prefix', 'th'+probe+'_', filename=data_file
        cdf_save_setting, 'time_range', time_range, filename=data_file
    endif
    event_info = cdf_read_setting(filename=data_file)
    event_info['data_file'] = data_file

    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']

;---Particle related vars.
    ; density_var = alfven_arc_load_themis_data_density(event_info)
    spec_vars = alfven_arc_load_themis_data_particle_data(event_info)


;---Orbit related vars.
    orbit_time_var = 'orbit_ut'
    time_step = 60.
    pad_time = 30.*60
    event_info['orbit_time_var'] = orbit_time_var
    event_info['orbit_time_step'] = time_step
    event_info['orbit_time_range'] = time_range+[-1,1]*pad_time
    times = alfven_arc_load_themis_data_orbit_ut(event_info, time_var=orbit_time_var)
    r_gsm_var = alfven_arc_load_themis_data_r_gsm(event_info, time_var=orbit_time_var)
    model_vars = alfven_arc_load_themis_data_load_model_vars(event_info, time_var=orbit_time_var)


;---Field related vars.
    field_time_var = 'field_ut'
    time_step = 1d/8
    pad_time = 30.*60
    event_info['field_time_var'] = field_time_var
    event_info['field_time_step'] = time_step
    event_info['field_time_range'] = time_range+[-1,1]*pad_time
    times = alfven_arc_load_themis_data_field_ut(event_info, time_var=field_time_var)
    b_gsm_var = alfven_arc_load_themis_data_b_gsm(event_info, time_var=field_time_var)
    e_dsl_var = alfven_arc_load_themis_data_e_dsl(event_info, time_var=field_time_var, _extra=ex)

    ; Seperate B0 and B1.
    event_info['b0_window'] = 10d
    bmod_var = alfven_arc_load_themis_data_bmod_gsm(event_info)
    b0_gsm_var = alfven_arc_load_themis_data_b0_gsm(event_info, time_var=field_time_var)
    b1_gsm_var = alfven_arc_load_themis_data_b1_gsm(event_info)
    
    ; Calculate Edot0.
    edot0_dsl_var = alfven_arc_load_themis_data_edot0_dsl(event_info, time_var=field_time_var)
    
    ; Convert data to FAC.
    fac_vars = alfven_arc_load_themis_data_fac_vars(event_info)

    ; Calculae pflux.
    pf_vars = alfven_arc_load_themis_data_pflux(event_info, fac_vars)


    return, event_info
end