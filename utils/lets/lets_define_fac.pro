;+
; b_var=.
; r_var=.
; fac_coord=.
; fac_labels=.
;-
function lets_define_fac, b_var=b_var, r_var=r_var, $
    fac_coord=fac_coord, fac_labels=fac_labels, $
    update=update, get_name=get_name, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var, $
    _extra=ex

    errmsg = ''
    retval = !null

    ; Need prefix and coord.
    if n_elements(b_var) eq 0 then begin
        errmsg = 'No input b_var ...'
        return, retval
    endif

    prefix = get_prefix(b_var)
    mission_probe = get_var_setting(r_var, 'mission_probe')
    probe_info = resolve_probe(mission_probe)
    time_range = get_var_setting(r_var, 'requested_time_range')

    coord = strlowcase(get_setting(b_var,'coord'))
    r_coord = strlowcase(get_setting(r_var,'coord'))
    if coord ne r_coord then begin
        errmsg = 'b_var and r_var should be in the same coord ...'
        return, retval
    endif
    if n_elements(suffix) eq 0 then suffix = ''
    ; need to use lets_read or lets_read_this to ensure have mission in setting.
    mission = get_var_setting(b_var, 'mission', exist)
    if n_elements(fac_coord) eq 0 then fac_coord = mission+'_fac'
    if n_elements(fac_labels) ne 3 then fac_labels = ['b','w','o']
    var_info = prefix+'q_'+coord+'2'+fac_coord+suffix
    if keyword_set(get_name) then return, var_info

    ; Check if update in memory.
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info


    ; Read from routine.
    if n_elements(time_var) eq 0 then time_var = r_var
    get_data, time_var, times

    bvec = get_var_data(b_var, at=times)
    rvec = get_var_data(r_var, at=times)

    rhat = sunitvec(rvec)
    bhat = sunitvec(bvec)
    what = sunitvec(vec_cross(rhat, bhat))
    ohat = vec_cross(bhat, what)

    ntime = n_elements(times)
    ndim = 3
    m_xxx2fac = fltarr(ntime,ndim,ndim)
    m_xxx2fac[*,0,*] = bhat
    m_xxx2fac[*,1,*] = what
    m_xxx2fac[*,2,*] = ohat
    q_xxx2fac = mtoq(m_xxx2fac)    
    store_data, var_info, times, q_xxx2fac

    coord_labels = get_var_setting(b_var, 'coord_labels')
    if n_elements(coord_labels) ne 3 then coord_labels = constant('xyz')
    add_setting, var_info, {$
        display_type: 'quaternion', $
        in_coord: coord, $
        in_coord_labels: coord_labels, $
        out_coord: fac_coord, $
        out_coord_labels: fac_labels}

    ; Save to file.
    if n_elements(data_file) ne 0 then begin
        is_success = save_var_to_file(var_info, file=data_file, time_var=time_var)
    endif
    
    return, var_info

end