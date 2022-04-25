;+
; Read given parameters for a given time range
; var=. Can be 'sw_b_gse','sw_v_gse','sw_n','sw_p','sw_t','symh','symd','asyh','asyd','au','al','sw_timeshift'.
;-

function ml_omni_read_param_var, input_time_range, var=input_var, $
    get_name=get_name, errmsg=errmsg

    errmsg = ''
    retval = ''

;---Return if only name is needed.
    prefix = 'omni_'
    var = prefix+input_var
    if keyword_set(get_name) then return, var

;---Check input time range.
    if n_elements(input_time_range) ne 2 then begin
        errmsg = 'Invalid input_time_range ...'
        return, retval
    endif
    time_range = time_double(input_time_range)


;---Read files.
    files = ml_omni_load_param(time_range, errmsg=errmsg)
    if errmsg ne '' then return, retval

;---Read vars.
    var_list = list()
    var_list.add, dictionary('in_vars', var)
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval
    vatt = cdf_read_setting(var, filename=files[0])

    is_vec_var = input_var eq 'sw_v_gse' or input_var eq 'sw_b_gse'
    display_type = (is_vec_var)? 'vector': 'scalar'
    short_name = strsplit(input_var,'_',extract=1)
    if is_vec_var then short_name = short_name[0:-2]
    short_name = strjoin(strupcase(short_name),' ')

    settings = dictionary($
        'display_type', display_type, $
        'unit', vatt['UNITS'], $
        'short_name', short_name )
    if is_vec_var then begin
        settings['coord'] = 'GSE'
        settings['coord_labels'] = constant('xyz')
    endif
    add_setting, var, smart=1, settings
    return, var


end

time_range = ['2012-10-01','2019-10-05']
var = ml_omni_read_param_var(time_range, var='pc_index')
end