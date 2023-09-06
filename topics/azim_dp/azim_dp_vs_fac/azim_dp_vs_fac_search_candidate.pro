;+
; Search candidates when dp vs fac is likely to occur.
;-


    time_range = time_double(['2008-01-01','2018-01-01'])
    
    if check_if_update('ae',time_range) then omni_read_index, time_range

    
    ae = get_var_data('ae', times=times)
    min_ae = 1000.
    index = where(ae ge min_ae, count)
    if count eq 0 then stop
    substorm_times = times[time_to_range(index,time_step=1)]
    substorm_durations = substorm_times[*,1]-substorm_times[*,0]
    min_duration = 120. ; sec.
    index = where(substorm_durations ge min_duration, nsubstorm)
    if nsubstorm eq 0 then stop
    
    substorm_times = substorm_times[index,*]    
    ntime = n_elements(times)
    time_step = times[1]-times[0]
    flags = fltarr(ntime)
    pad_time = 3600.    ; sec.
    for substorm_id=0,nsubstorm-1 do begin
        the_time_range = reform(substorm_times[substorm_id,*])+[-1,1]*pad_time
        index = where_pro(times, '[]', the_time_range, count=count)
        if count eq 0 then continue
        flags[index] = 1
    endfor
    index = where(flags eq 1)
    candidate_times = times[time_to_range(index,time_step=1)]
    ncandidate = n_elements(candidate_times)*0.5
    candidate_list = list()
    for candidate_id=0,ncandidate-1 do begin
        candidate_list.add, reform(candidate_times[candidate_id,*])
    endfor
    
    
    event_list = list()
    out_file = join_path([homedir(),'azim_dp_vs_fac_substorm_list.txt'])
    ftouch, out_file
    foreach the_time_range, candidate_list do begin
        lprmsg, time_string(the_time_range)
        roi_list = azim_dp_vs_fac_search_roi(the_time_range, reset=1)
        mlt_list = azim_dp_vs_fac_search_mlt(roi_list, reset=1)
        if n_elements(mlt_list) eq 0 then continue
        event_list.add, mlt_list, /extract
        foreach tmp, mlt_list do azim_dp_candidate_write, tmp, filename=out_file
    endforeach


end