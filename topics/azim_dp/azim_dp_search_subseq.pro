;+
; Search for subsequence, adopted from azim_df_search_subgroup.
;-

function azim_dp_search_subseq, event_list, min_dp_count=min_dp_count


    subseq_list = list()
    if n_elements(event_list) eq 0 then return, subseq_list
    if n_elements(min_dp_count) eq 0 then min_dp_count = 4

    foreach event, event_list do begin
        dp_list = event.dp_list
        if dp_list.length lt min_dp_count then continue

        ; Sort in ramp_time.
        ramp_times = dblarr(dp_list.length)
        foreach dp, dp_list, dp_id do ramp_times[dp_id] = dp.time
        index = sort(ramp_times)
        dp_list = dp_list[index]

        time_range = minmax(ramp_times)
        if keyword_set(test_time) then if product(time_range-test_time) ge 0 then continue else stop

    ;---Find consecutive different probes as a subsequence.
        the_subseq_list = list()
        probe_list = list()
        the_dp_list = list()
        foreach df, dp_list, df_id do begin
            probe_pos = probe_list.where(df.probe)
            if probe_pos eq !null then begin
                ; Add the current probe to list if it doesn't exist.
                probe_list.add, df.probe
                the_dp_list.add, df
            endif else begin
                probe_pos = probe_pos[0]
                ; The current probe exists in the current list.
                if probe_list.length lt min_dp_count then begin
                    ; If the current list is not long enough then toss and start over.
                    for iter=0, probe_pos do begin
                        probe_list.remove, 0
                        the_dp_list.remove, 0
                    end
                endif else begin
                    ; The current list is long enough, pop it.
                    current_list = list()   ; must copy manually
                    for iter=0, the_dp_list.length-1 do current_list.add, the_dp_list[iter]
                    the_subseq_list.add, current_list
                    ; Remove probes before the current probe.
                    for iter=0, probe_pos do begin
                        probe_list.remove, 0
                        the_dp_list.remove, 0
                    endfor
                endelse
                ; Add the new probe to the list.
                probe_list.add, df.probe
                the_dp_list.add, df
            endelse
            if df_id eq dp_list.length-1 and probe_list.length ge min_dp_count then begin
                ; Pop the current list.
                the_subseq_list.add, the_dp_list
                probe_list = list()
                the_dp_list = list()
            endif
        endforeach


    ;---Coerce to subseq_list.
        if the_subseq_list.length eq 0 then continue
        foreach dp_list, the_subseq_list do begin
            ramp_times = dblarr(dp_list.length)
            foreach dp, dp_list, dp_id do ramp_times[dp_id] = dp.time
            time_range = minmax(ramp_times)
            probes = strarr(dp_list.length)
            foreach dp, dp_list, dp_id do probes[dp_id] = dp.probe

            subseq_list.add, dictionary($
                'time_range', time_range, $
                'probes', probes, $
                'dp_list', dp_list, $
                'edge_list', list(), $
                'triad_list', list() )
        endforeach
    endforeach

    return, subseq_list

end
