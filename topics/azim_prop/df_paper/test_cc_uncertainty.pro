;+
; The confidence interval of the lag-c.c curve is simply 1/sqrt(N).
; Thus, if a width is assigned to the proper lag, we can calculate the
; ratio of the c.c corresponds to the width and the maximum c.c.
; Then the ratio can be used to derive the percentage of confidence interval.
;-

function xcorr_z2c, z
    return, imsl_erf(z/sqrt(2))
end

function xcorr_c2z, c
    return, imsl_erf(c,/inverse)*sqrt(2)
end

pro test_cc_uncertainty

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    test = 1

;---Settings.
    all_time_lag_error = []
    event_ids = (events.keys()).sort()
    foreach event_id, event_ids do begin
        lprmsg, 'Processing '+event_id+' ...'
        if keyword_set(test_id) then if event_id ne test_id then continue

        event = events[event_id]
        probes = event.probes
        nprobe = n_elements(probes)
        time_lag_error = dblarr(nprobe)
        foreach probe, probes, ii do time_lag_error[ii] = event[probe].time_lag_error
        index = where(time_lag_error ne 0)
        time_lag_error = time_lag_error[index]
        all_time_lag_error = [all_time_lag_error,time_lag_error]
    endforeach

    print, mean(all_time_lag_error)
    stop


end
