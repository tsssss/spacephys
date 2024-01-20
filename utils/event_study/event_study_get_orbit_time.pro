;+
; Need to have orbit_time_range, orbit_time_step, data_file
;-
function event_study_get_orbit_time, event_info, time_var=time_var, get_name=get_name, update=update, suffix=suffix

    id = 'orbit'
    pad_time = 1800d
    key = id+'_time_range'
    if ~event_info.haskey(key) then event_info[key] = event_info['time_range']+[-1,1]*pad_time
    key = id+'_time_step'
    if ~event_info.haskey(key) then event_info[key] = 60d

    return, event_study_get_time(event_info, id=id, time_var=time_var, get_name=get_name, update=update, suffix=suffix)

end