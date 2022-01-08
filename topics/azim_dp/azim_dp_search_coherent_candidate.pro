;+
; Search coherent candidate.
;-

function azim_dp_search_coherent_candidate, subseq_list, min_dp_count=min_dp_count

    coherent_candidate_list = list()
    if n_elements(subseq_list) eq 0 then return, coherent_candidate_list
    if n_elements(min_dp_count) eq 0 then min_dp_count = 4

    dtime_max = 30.*60      ; sec.
    dtime_mean_max = 8.*60  ; sec.
    nedge_vertex = 2        ; #.
    ntriad_vertex = 3       ; #.
    min_mean_xcor = 0.5     ; #.
    triad_angle_range = [15.,165]   ; deg.
    max_angle_diff = 30.    ; deg.
    max_vmag_ratio = 0.2    ; #.
    min_triad_count = 2     ; #.



    foreach event, subseq_list do begin
    ;---Check coherency.
        dp_list = event.dp_list
        ndp = dp_list.length
        ramp_times = dblarr(ndp)
        ramp_widths = fltarr(ndp)
        foreach dp, dp_list, dp_id do ramp_times[dp_id] = dp.time
        foreach dp, dp_list, dp_id do ramp_widths[dp_id] = dp.width
        probes = event.probes
        nprobe = n_elements(probes)
        probe_list = list(probes, /extract)

        dtimes = ramp_times[1:-1]-ramp_times[0:-2]
        if max(dtimes) ge dtime_max then continue
        if mean(dtimes) ge dtime_mean_max then continue

        ; xcor.
        edge_combos = list()
        for probe_id=0,nprobe-2 do edge_combos.add, probes[probe_id:probe_id+1]
        edge_list = list()
        common_data_rate = 10.
        time_range = event.time_range+[-1,1]*3600
        common_times = make_bins(time_range, common_data_rate)
        foreach probe, probes do begin
            prefix = probe+'_'
            azim_dp_read_theta, time_range, probe=probe
            interp_time, prefix+'theta', common_times
        endforeach

        foreach edge_probes, edge_combos do begin
            edge_probes = edge_probes[sort(edge_probes)]

            ; Get info about the vertex.
            probe_index = intarr(nedge_vertex)
            foreach probe, edge_probes, ii do probe_index[ii] = where(probes eq probe)
            vertex_ramp_times = ramp_times[probe_index]
            vertex_ramp_widths = ramp_widths[probe_index,*]

            ; Sort by obs_time.
            index = sort(vertex_ramp_times)
            edge_probes = edge_probes[index]
            vertex_ramp_times = vertex_ramp_times[index]
            vertex_ramp_widths = vertex_ramp_widths[index]
            ; Must be positive or 0.
            time_lag_ramp = total(vertex_ramp_times*[-1,1])

        ;---xcor.
            ; This setting works.
            the_width = min(vertex_ramp_widths)
            the_time_range = the_width*[-1,1]*2.5
            nlag = floor(the_width*2/common_data_rate)
            section_times = make_bins(the_time_range, common_data_rate)
            nsection_time = n_elements(section_times)

            ; Read data.
            nrec = n_elements(section_times)
            data = fltarr(nrec,nedge_vertex)
            for probe_id=0,nedge_vertex-1 do begin
                the_times = vertex_ramp_times[probe_id]+section_times
                the_probe = edge_probes[probe_id]
                the_var = the_probe+'_theta'
                data[*,probe_id]= get_var_data(the_var, at=the_times)
            endfor

            ; Do xcor.
            lags = findgen(nlag)-round(nlag/2)
            xcors = c_correlate(data[*,0], data[*,1], lags)
            xcor_max = max(xcors, index)
            xcor_time_lag = lags[index]*common_data_rate
            ; xcor error.
            xx = total(data,2)/nedge_vertex
            del_x = xx-mean(xx)
            dx_dt = deriv(xx)/common_data_rate
            xcor_err = sqrt(1./(nrec-1)*(1-xcor_max)/xcor_max*2*mean(del_x^2)/mean(dx_dt^2))

            ; Save results.
            time_lag_ramp = total(vertex_ramp_times*[-1,1])
            time_lag_xcor = time_lag_ramp+xcor_time_lag
            time_lag_diff = abs(time_lag_ramp-time_lag_xcor)
            time_lag_ratio = time_lag_diff/max(abs([time_lag_ramp,time_lag_xcor]))

            edge_info = dictionary($
                'probes', edge_probes, $
                'time_lag_ramp', time_lag_ramp, $
                'time_lag_xcor', time_lag_xcor, $
                'xcor_max', xcor_max, $
                'xcor_err', xcor_err, $
                'time_lag_diff', time_lag_diff, $
                'time_lag_ratio', time_lag_ratio)
            edge_list.add, edge_info
        endforeach

        xcor_maxs = list()
        foreach edge_info, edge_list do xcor_maxs.add, edge_info.xcor_max
        xcor_maxs = xcor_maxs.toarray()
        mean_xcor = mean(xcor_maxs)
        mean_xcor = round(mean_xcor*10)/10.
        if mean_xcor le min_mean_xcor then continue
        event['edge_list'] = edge_list


    ;---Check triad.
        ndim = 3
        ramp_r_sms = fltarr(nprobe,ndim)
        foreach dp, dp_list, dp_id do ramp_r_sms[dp_id,*] = dp.r_sm
        xcor_times = dblarr(nprobe)
        xcor_times[0] = dp_list[0].time
        edge_list = event.edge_list
        foreach edge, edge_list, edge_id do xcor_times[edge_id+1] = xcor_times[edge_id]+edge.time_lag_xcor

        ntriad_vertex = 3
        triad_combos = choose_from(probes, ntriad_vertex)
        triad_list = list()
        foreach triad_probes, triad_combos do begin
            sorted_probes = triad_probes[sort(triad_probes)]
            key = strjoin(sorted_probes,'_')
            lprmsg, 'triad: '+key, log_file

            ; Get info about the vertex.
            probe_index = []
            foreach probe, triad_probes do probe_index = [probe_index, probe_list.where(probe)]
            the_ramp_times = ramp_times[probe_index]
            the_xcor_times = xcor_times[probe_index]
            the_r_sms = ramp_r_sms[probe_index,*]

            ; Sort by obs_time.
            index = sort(the_ramp_times)
            triad_probes = triad_probes[index]
            the_ramp_times = the_ramp_times[index]
            the_xcor_times = the_xcor_times[index]
            the_r_sms = the_r_sms[index,*]
            center_r_sm = reform(total(the_r_sms,1)/ntriad_vertex)

            triad = dictionary($
                'probes', triad_probes, $
                'times', the_ramp_times, $
                'r_sms', the_r_sms, $
                'center_r_sm', center_r_sm )


            ; Angles.
            triad_angles = triangle_angles(reform(transpose(the_r_sms), [1,ntriad_vertex,ndim]))
            triad_angles = reform(triad_angles)
            index = lazy_where(triad_angles, '()', triad_angle_range, count=count)
            if count ne ndim then begin
                lprmsg, 'bad geometry, skip ...', log_file
                continue
            endif
            triad['angles'] = triad_angles

            ; Solve for v_2d.
            timing = azim_df_solve_triad_timing(the_xcor_times, the_r_sms)
            if n_elements(timing) eq 0 then begin
                lprmsg, 'no solution for the_xcor_times, skip ...', log_file
                continue
            endif
            triad['vhat_xcor'] = timing.vhat
            triad['vmag_xcor'] = timing.vmag
            triad['omega_xcor'] = timing.omega

            timing = azim_df_solve_triad_timing(the_ramp_times, the_r_sms)
            if n_elements(timing) eq 0 then begin
                lprmsg, 'no solution for the_ramp_times, skip ...', log_file
                continue
            endif
            triad['vhat_ramp'] = timing.vhat
            triad['vmag_ramp'] = timing.vmag
            triad['omega_ramp'] = timing.omega

            ; Filter by vmag and vhat diff.
            vmag_xcor = triad.vmag_xcor
            vmag_ramp = triad.vmag_ramp
            vmag_ratio = abs(vmag_xcor-vmag_ramp)/vmag_ramp
            vmag_ratio = round(vmag_ratio*10)/10.
            if vmag_ratio gt max_vmag_ratio then begin
                lprmsg, 'vmag diff too large, skip ...', log_file
                continue
            endif

            vhat_xcor = triad.vhat_xcor
            vhat_ramp = triad.vhat_ramp
            vhat_angle = sang(vhat_xcor,vhat_ramp, /deg)
            vhat_angle = round(vhat_angle)
            if vhat_angle gt max_angle_diff then begin
                lprmsg, 'vhat diff too large, skip ...', log_file
                continue
            endif

            triad_list.add, triad
        endforeach

        lprmsg, 'Found '+string(triad_list.length,format='(I0)')+' triads with good timing ...', log_file
        if triad_list.length lt min_triad_count then continue
        event.triad_list = triad_list

        coherent_candidate_list.add, event
    endforeach

    return, coherent_candidate_list
end
