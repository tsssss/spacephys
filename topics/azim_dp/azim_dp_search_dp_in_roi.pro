;+
; Search DP in roi_list, c.f. azim_dp_search_roi.
;-

function azim_dp_search_dp_in_roi, roi_list, min_dp_count=min_dp_count

    event_list = list()
    if n_elements(roi_list) eq 0 then return, event_list
    if n_elements(min_dp_count) eq 0 then min_dp_count = 4

    foreach roi, roi_list do begin
        dp_list = list()
        foreach roi_time_range, roi.time_range_list, roi_id do begin
            probes = roi.probe_list[roi_id]
            foreach probe, probes do begin
                the_dp_list = azim_dp_read_dp(roi_time_range, probe=probe)
                if n_elements(the_dp_list) eq 0 then continue
                foreach dp, the_dp_list do begin
                    dp['probe'] = probe
                    dp_list.add, dp
                endforeach
            endforeach
        endforeach
        
        if dp_list.length lt min_dp_count then continue
        probe_list = list()
        foreach dp, dp_list do begin
            probe = dp.probe
            if probe_list.where(probe) ne !null then continue
            probe_list.add, probe
        endforeach
        if probe_list.length lt min_dp_count then continue
        
        ramp_times = dblarr(dp_list.length)
        foreach dp, dp_list, dp_id do ramp_times[dp_id] = dp.time
        time_range = minmax(ramp_times)
        
        event_list.add, dictionary($
            'time_range', time_range, $
            'probes', probe_list.toarray(), $
            'dp_list', dp_list, $
            'edge_list', list(), $
            'triad_list', list() )
    endforeach
    
    return, event_list

end
