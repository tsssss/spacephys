;+
; Read RBSP dis in Re.
;-

function ml_rbsp_read_dis, input_time_range, probe=probe, errmsg=errmsg, get_name=get_name

    errmsg = ''
    retval = ''

;---Check input probe.
    probes = rbsp_probes()
    index = where(probes eq probe, count)
    if count eq 0 then begin
        errmsg = 'probe is unknown ...'
        return, retval
    endif
    prefix = 'rbsp'+probe+'_'

;---Return if only name is needed.
    dis_var = prefix+'dis'
    if keyword_set(get_name) then return, dis_var

;---Check input time range.
    if n_elements(input_time_range) ne 2 then begin
        errmsg = 'Invalid input_time_range ...'
        return, retval
    endif
    time_range = time_double(input_time_range)

;---Read r and calculate |r|.
    r_var = ml_rbsp_read_pos(input_time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval
    get_data, r_var, times, r_vec
    store_data, dis_var, times, snorm(r_vec)
    add_setting, dis_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'unit', 'Re', $
        'short_name', '|R|' )
    return, dis_var

end