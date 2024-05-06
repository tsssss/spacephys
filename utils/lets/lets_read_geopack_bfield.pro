

function read_geopack_bfield, orbit_var, external_model=external_model, internal_model=internal_model, $
    var_info=var_info, get_name=get_name, t89_use_kp=t89_use_kp, update=update, _extra=ex

    prefix = get_prefix(orbit_var)
    coord_orig = 'gsm'
    if n_elements(external_model) eq 0 then external_model = 't89'
    if n_elements(internal_model) eq 0 then internal_model = 'dipole'
    if n_elements(suffix) eq 0 then begin
        tmp = [internal_model,external_model]
        index = where(tmp ne 'n/a', count)
        if count eq 0 then begin
            suffix = ''
        endif else if count eq 1 then begin
            suffix = '_'+tmp[index[0]]
        endif else begin
            suffix = '_'+internal_model+'_'+external_model
        endelse
    endif
    if n_elements(var_info) eq 0 then var_info = prefix+'b_'+coord_orig+suffix
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

;---Internal and external model.
    b_internal = fltarr(ntime,ndim)
    b_external = fltarr(ntime,ndim)
    foreach time, times, time_id do begin
        ps = geopack_recalc(time)
        rx = r_gsm[time_id,0]
        ry = r_gsm[time_id,1]
        rz = r_gsm[time_id,2]

    ;---Internal model.
        if internal_model eq 'dipole' or internal_model eq 'dip' then begin
            geopack_dip, rx,ry,rz, bx,by,bz
        endif else if internal_model eq 'igrf' then begin
            geopack_igrf_gsm, rx,ry,rz, bx,by,bz
        endif
        b_internal[time_id,*] = [bx,by,bz]


    ;---External model.
        if ~has_external_model then continue
        if external_model eq 't04s' then begin
            routine = 'geopack_ts04'
        endif else routine = 'geopack_'+external_model
        par = reform(pars[time_id,*])
        call_procedure, routine, par, rx,ry,rz, dbx,dby,dbz
        b_external[time_id,*] = [dbx,dby,dbz]
    endforeach

    b_model = b_internal+b_external
    settings = get_var_setting(orbit_var)
    settings['short_name'] = strupcase(internal_model)+'+'+strupcase(external_model)+' B'
    settings['unit'] = 'nT'
    settings['coord'] = coord_orig
    settings['internal_model'] = internal_model
    settings['external_model'] = external_model
    return, save_data_to_memory(var_info, times, b_model, settings=settings)

end


function lets_read_geopack_bfield, var_info=var_info, $
    orbit_var=orbit_var, external_model=external_model, internal_model=internal_model, $
    t89_use_kp=t89_use_kp, $
    update=update, get_name=get_name, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var, $
    coord=coord, _extra=ex

    errmsg = ''
    retval = ''

    ; Need prefix
    if n_elements(orbit_var) eq 0 then begin
        errmsg = 'No input orbit_var ...'
        return, retval
    endif
    mission_probe = get_var_setting(orbit_var, 'mission_probe')
    probe_info = resolve_probe(mission_probe)
    prefix = probe_info['prefix']


    ; Get var_info.
    if n_elements(external_model) eq 0 then external_model = 't89'
    if n_elements(internal_model) eq 0 then internal_model = 'dipole'
    igrf = internal_model eq 'igrf'
    if n_elements(suffix) eq 0 then begin
        tmp = [internal_model,external_model]
        index = where(tmp ne 'n/a', count)
        if count eq 0 then begin
            suffix = ''
        endif else if count eq 1 then begin
            suffix = '_'+tmp[index[0]]
        endif else begin
            suffix = '_'+internal_model+'_'+external_model
        endelse
    endif
    default_coord = 'gsm'
    if n_elements(coord) eq 0 then coord = default_coord
    if n_elements(var_info) eq 0 then var_info = prefix+'b_'+coord+suffix
    if keyword_set(get_name) then return, var_info

    ; Check if update in memory.
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info

    ; Read var from routine.
    vec_default_var = read_geopack_bfield(orbit_var, var_info=var_info, $
        external_model=external_model, internal_model=internal_model, t89_use_kp=t89_use_kp, update=update)
    
    ; Convert to the wanted coord.
    if var_info ne vec_default_var then begin
        coord_msg = [default_coord,coord]
        var_info = lets_cotran(coord_msg, input=vec_default_var, output=var_info)
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
    print, lets_read_orbit(time_range, probe=probe, time_var=time_var, coord='gsm')
    print, lets_read_geopack_bfield(time_range, probe=probe, save_to=data_file, time_var=time_var, coord='gsm', internal='igrf')
    print, lets_read_bfield(time_range, probe=probe)
    stop
endforeach
end
