;+
; Need to have orbit_time_range, orbit_time_step, data_file
;-
function event_study_get_time, event_info, id=id, time_var=time_var, get_name=get_name, update=update, suffix=suffix

    if n_elements(id) eq 0 then message, 'Inconsistency ...'
    if n_elements(suffix) eq 0 then suffix = ''
    if n_elements(time_var) eq 0 then time_var = event_info['prefix']+id+'_time'+suffix
    if keyword_set(get_name) then return, time_var

    data_file = event_info['data_file']
    if keyword_set(update) then cdf_del_var, time_var, filename=data_file
    if ~cdf_has_var(time_var, filename=data_file) then begin
        time_range = event_info[id+'_time_range']
        time_step = event_info[id+'_time_step']
        times = make_bins(time_range, time_step)
        cdf_save_var, time_var, value=times, filename=data_file
        settings = dictionary($
            'time_step', time_step, $
            'data_time_range', time_range )
        cdf_save_setting, settings, filename=data_file, varname=time_var
    endif
    return, cdf_read_var(time_var, filename=data_file)

end