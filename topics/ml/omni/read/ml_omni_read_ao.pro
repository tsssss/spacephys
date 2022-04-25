;+
; Read AO = (AU+AL)*0.5.
;-

function ml_omni_read_ao, input_time_range, get_name=get_name

    errmsg = ''
    retval = ''

;---Return if only name is needed.
    prefix = 'omni_'
    var = prefix+'ao'
    if keyword_set(get_name) then return, var

;---Check input time range.
    if n_elements(input_time_range) ne 2 then begin
        errmsg = 'Invalid input_time_range ...'
        return, retval
    endif
    time_range = time_double(input_time_range)

    au_var = ml_omni_read_param_var(input_time_range, var='au')
    al_var = ml_omni_read_param_var(input_time_range, var='al')
    get_data, au_var, times, au
    get_data, al_var, times, al
    store_data, var, times, (au+al*0.5
    add_setting, var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'AO', $
        'unit', 'nT' )
    return, var

end