;+
; Search coherent event, adopted from azim_df_filter_coherent_df.
;-
function azim_dp_search_coherent_event, coherent_candidate_list, min_dp_count=min_dp_count

    coherent_event_list = list()
    if n_elements(coherent_candidate_list) eq 0 then return, coherent_event_list


    max_angle_diff = 15.    ; deg.
    max_vmag_ratio = 0.3    ; #.

    foreach event, coherent_candidate_list do begin
        triad_list = event.triad_list
        ntriad = triad_list.length

    ;---Check vmag.
        vmags = fltarr(ntriad)
        foreach triad, triad_list, ii do vmags[ii] = triad.vmag_ramp
        vmag_mean = mean(vmags)
        vmag_stddev = stddev(vmags)
        vmag_ratio = vmag_stddev/vmag_mean
        vmag_ratio = round(vmag_ratio*10)/10.
        if vmag_ratio gt max_vmag_ratio then continue

    ;---Check vhat.
        vhat_angles = fltarr(ntriad)
        foreach triad, triad_list, ii do vhat_angles[ii] = atan(triad.vhat_ramp[1],triad.vhat_ramp[0])
        vhat_angles *= constant('deg')
        vhat_angle_stddev = stddev(vhat_angles)
        vhat_angle_stddev = round(vhat_angle_stddev)
        if vhat_angle_stddev gt max_angle_diff then continue

        coherent_event_list.add, event
    endforeach

    return, coherent_event_list

end