;+
; Adopted from azim_df_filter_vertex.
;-

function azim_dp_is_good_ramp, ramp

    retval = 0
    if n_elements(ramp) eq 0 then return, retval

;---Settings.
    height_range = [0.,80]      ; deg.
    width_range = [1.,25]*60.   ; sec.
    min_scaled_height = 8.      ; deg.
    max_abs_scaled_theta = 80.  ; deg.
    scale_width = 4.            ; hour.
    section_time_range = [0.,10]*60.    ; sec.
    section_min_theta = 0.              ; deg.
    section_min_ratio = 0.4             ; #.
    section_total_ratio = 0.8           ; #.


;---Filter.
    ; width.
    index = lazy_where(ramp.width, '()', width_range, count=count)
    if count eq 0 then return, retval

    ; height.
    index = lazy_where(ramp.height, '()', height_range, count=count)
    if count eq 0 then return, retval

    ; scaled_height.
    scaled_height = azim_df_scale_theta(ramp.height, ramp.mlt, width=scale_width)
    scaling_factor = scaled_height/ramp.height

    scaled_height = float(round(scaled_height))
    if scaled_height le min_scaled_height then return, retval

    scaled_theta_range = ramp.theta_range*scaling_factor
    if max(abs(scaled_theta_range)) ge max_abs_scaled_theta then return, retval

    ; shape.
    probe = ramp.probe
    the_time_range = ramp.time+section_time_range
    azim_dp_read_theta, the_time_range, probe=probe

    prefix = probe+'_'
    theta_var = prefix+'theta'
    theta = get_var_data(theta_var, times=times)
    index = where(theta ge section_min_theta, count)
    if count eq 0 then return, retval
    time_ranges = times[time_to_range(index,time_step=1)]
    durations = time_ranges[*,1]-time_ranges[*,0]
    section_duration = total(section_time_range*[-1,1])
    ; total duration.
    total_duration = total(durations)
    section_ratio = total_duration/section_duration
    section_ratio = round(section_ratio*10)/10.
    if section_ratio ge section_total_ratio then return, 1
    
    ; max duration.
    max_duration = max(durations, section_index)
    max_duration_time_range = reform(time_ranges[section_index,*])
    section_ratio = max_duration/section_duration
    section_ratio = round(section_ratio*10)/10.
    if section_ratio lt section_min_ratio then return, retval


;---Done.
    return, 1


end
