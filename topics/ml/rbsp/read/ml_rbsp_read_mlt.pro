;+
; Read RBSP MLT in hr.
;-

function ml_rbsp_read_mlt, input_time_range, probe=probe, errmsg=errmsg, get_name=get_name

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
    mlt_var = prefix+'mlt'
    if keyword_set(get_name) then return, mlt_var

;---Check input time range.
    if n_elements(input_time_range) ne 2 then begin
        errmsg = 'Invalid input_time_range ...'
        return, retval
    endif
    time_range = time_double(input_time_range)

;---Read files.
    files = ml_rbsp_load_orbit_var(time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval

;---Read vars.
    var_list = list()
    var_list.add, dictionary('in_vars', mlt_var)
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    add_setting, mlt_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'unit', 'h', $
        'short_name', 'MLT' )
    return, mlt_var

end
