function trace_to_ionosphere, var_info=var_info, hemisphere=hemisphere, stop_altitude=h0, $
    orbit_var=orbit_var, external_model=external_model, internal_model=internal_model, $
    get_name=get_name, $
    t89_use_kp=t89_use_kp, refine=refine, _extra=ex

    if n_elements(orbit_var) eq 0 then begin
        errmsg = 'No input orbit_var ...'
        return, retval
    endif
    mission_probe = get_var_setting(orbit_var, 'mission_probe')
    probe_info = resolve_probe(mission_probe)
    prefix = probe_info['prefix']
    coord_orig = 'gsm'
    if n_elements(hemisphere) eq 0 then hemisphere = 'north'
    if n_elements(h0) eq 0 then h0 = 100d   ; km.
    if n_elements(external_model) eq 0 then external_model = 't89'
    if n_elements(internal_model) eq 0 then internal_model = 'dipole'
    suffix = '_'+internal_model+'_'+external_model+'_'+hemisphere

    if n_elements(var_info) eq 0 then var_info = prefix+'f_'+coord_orig+suffix
    if keyword_set(get_name) then return, var_info


    get_data, orbit_var, times, r_coord
    coord_in = get_var_setting(orbit_var, 'coord')
    if strlowcase(coord_in) ne coord_orig then begin
        coord_msg = strlowcase([coord_in,coord_orig])
        probe = get_var_setting(orbit_var, 'probe')
        r_gsm = cotran(r_coord, times, coord_msg=coord_msg, probe=probe, _extra=ex)
    endif else begin
        r_gsm = temporary(r_coord)
    endelse
    if n_elements(time_range) eq 0 then time_range = minmax(times)
    if n_elements(time_range) eq 1 then time_range = minmax(times)
    if keyword_set(update) then is_success = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ndim = 3
    ntime = n_elements(times)
    f_gsm = fltarr(ntime,ndim)

;---Prepare external model parameters.
    has_external_model = 0
    avail_external_models = ['t89','t96','t01','t04s','ts04']
    index = where(strlowcase(external_model) eq avail_external_models, count)
    if count ne 0 then has_external_model = 1
    if has_external_model then begin
        t89_par = keyword_set(t89_use_kp)? !null: 2d
        par_time_range = time_range+[-1,1]*300
        par_var = geopack_read_par(par_time_range, model=external_model, t89_par=t89_par)
        pars = get_var_data(par_var, at=times)
    endif
    trace_dir = (hemisphere eq 'north')? -1: 1
    r0 = h0/constant('re')+1

    tmp = geopack_resolve_model(external_model)
    t89 = tmp.t89
    t96 = tmp.t96
    t01 = tmp.t01
    ts04 = tmp.ts04
    storm = tmp.storm
    igrf = (internal_model eq 'igrf')? 1: 0
    if ~keyword_set(refine) then refine = 1

    foreach time, times, time_id do begin
        ps = geopack_recalc(time)
        rx = r_gsm[time_id,0]
        ry = r_gsm[time_id,1]
        rz = r_gsm[time_id,2]

        geopack_trace, rx,ry,rz, trace_dir, reform(pars[time_id,*]), $
            fx,fy,fz, r0=r0, refine=refine, ionosphere=1, $
            t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm, igrf=igrf
        f_gsm[time_id,*] = [fx,fy,fz]
    endforeach

    settings = inherit_setting(get_var_setting(orbit_var),id='basic')
    settings['short_name'] = 'R'
    settings['coord'] = coord_orig
    settings['internal_model'] = internal_model
    settings['external_model'] = external_model
    settings['h0'] = h0
    settings['hemisphere'] = hemisphere
    return, save_data_to_memory(var_info, times, f_gsm, settings=settings)

end


function lets_trace_to_ionosphere, var_info=var_info, $
    hemisphere=hemisphere, stop_altitude=h0, $
    orbit_var=orbit_var, external_model=external_model, internal_model=internal_model, $
    t89_use_kp=t89_use_kp, $
    update=update, get_name=get_name, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var, $
    coord=coord, _extra=ex

    errmsg = ''
    retval = !null

    ; Need prefix
    if n_elements(orbit_var) eq 0 then begin
        errmsg = 'No input orbit_var ...'
        return, retval
    endif
    mission_probe = get_var_setting(orbit_var, 'mission_probe')
    probe_info = resolve_probe(mission_probe)
    prefix = probe_info['prefix']

    ; Get var_info.
    if n_elements(hemisphere) eq 0 then hemisphere = 'north'
    if n_elements(h0) eq 0 then h0 = 100d   ; km.
    if n_elements(external_model) eq 0 then external_model = 't89'
    if n_elements(internal_model) eq 0 then internal_model = 'dipole'
    igrf = internal_model eq 'igrf'
    if n_elements(suffix) eq 0 then suffix = '_'+internal_model+'_'+external_model+'_'+hemisphere
    default_coord = 'gsm'
    if n_elements(coord) eq 0 then coord = default_coord
    if n_elements(var_info) eq 0 then var_info = prefix+'f_'+coord+suffix
    if keyword_set(get_name) then return, var_info

    ; Check if update in memory.
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info

    ; Read var from routine.
    vec_default_var = trace_to_ionosphere(var_info=var_info, $
        orbit_var=orbit_var, $
        hemisphere=hemisphere, stop_altitude=h0, $
        external_model=external_model, internal_model=internal_model, t89_use_kp=t89_use_kp)
    
    ; Convert to the wanted coord.
    if default_coord ne coord then begin
        coord_msg = [default_coord,coord]
        vec_default_var = prefix+'f_'+default_coord+suffix
        vec_coord_var = var_info
        var_info = lets_cotran(coord_msg, input=vec_default_var, output=vec_coord_var)
    endif
    if n_elements(time_var) ne 0 then is_success = interp_var_to_time(var_info, time_var=time_var)
    
    ; Save to file.
    if n_elements(data_file) ne 0 then begin
        is_success = save_var_to_file(var_info, file=data_file, time_var=time_var)
    endif
    
    return, var_info

end

time_range = ['2015-03-17','2015-03-18']
probe = 'rbspb'
data_file = join_path([homedir(),'test.cdf'])
time_var = 'rbsp_orbit_time'
;if file_test(data_file) eq 1 then file_delete, data_file
foreach probe, ['rbspa','rbspb'] do begin
    orbit_var = lets_read_orbit(time_range, probe=probe, time_var=time_var, coord='gsm')
    fpt_var = lets_trace_to_ionosphere(orbit_var=orbit_var, coord='gsm', internal='igrf')
    prefix = get_var_prefix(orbit_var)
    fpt_mlat_vars = lets_read_mlat_vars(orbit_var=fpt_var, prefix=prefix+'f')
    fb_var = lets_read_geopack_bfield(orbit_var=fpt_var, external_model='n/a', var_info=prefix+'fb_gsm')
    b_var = lets_read_bfield(time_range, probe=probe)
    bmod_var = lets_read_geopack_bfield(orbit_var=orbit_var)
    b_vars = lets_decompose_bfield(b_var=b_var, bmod_var=bmod_var)
    cmap_var = lets_read_cmap(b_var=b_vars['b0'], fb_var=fb_var, update=1)
endforeach
end