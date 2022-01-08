;+
; Search subgroups containing large DF from different sc.
;-

function azim_df_search_subgroup, search_setting, project=project, $
    reset=reset, test_time=test_time

    retval = !null
    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()

;---Settings for selecting subgroup.
    the_key = 'search_subgroup'
    if ~search_setting.haskey(the_key) then message, 'No settings for '+the_key+' ...'
    search_info = search_setting[the_key]


;---Check if file exists, to avoid search again.
    file_suffix = search_info.file_suffix
    out_file = join_path([project.data_dir,file_suffix])
    log_file = join_path([project.data_dir,file_suffix+'.log'])
    if keyword_set(reset) then begin
        lprmsg, 'Resetting '+the_key+'...'
        file_delete, out_file, /allow_nonexistent
        file_delete, log_file, /allow_nonexistent
    endif
    if file_test(out_file) eq 1 then begin
        if file_lines(out_file) eq 0 then file_delete, out_file
    endif
    if file_test(log_file) eq 1 then begin
        if file_lines(log_file) eq 0 then file_delete, log_file
    endif


    if file_test(out_file) eq 0 or keyword_set(test_time) then begin
        event_id = 0

        probe_min_count = project.probe_min_count
        triad_angle_range = project.triad_angle_range
        triad_min_count = project.triad_min_count
        dtime_max = 30.*60
        dtime_mean_max = 8.*60

        if file_test(out_file) eq 0 then ftouch, out_file
        if file_test(log_file) eq 0 then ftouch, log_file
        if keyword_set(test_time) then log_file = -1    ; output to console.

        lprmsg, 'Search subgroups ...', log_file
        lprmsg, tab+'Min # of consecutive probes: '+string(probe_min_count,format='(I0)'), log_file
        lprmsg, tab+'Triad angle range (deg): ['+strjoin(string(triad_angle_range,format='(I0)'),',')+']', log_file
        lprmsg, tab+'Triad min count (#): '+string(triad_min_count,format='(I0)'), log_file
        lprmsg, tab+'Max edge dtime (sec): '+string(dtime_max,format='(I0)'), log_file
        lprmsg, tab+'Max edge dtime mean (sec): '+string(dtime_mean_max,format='(I0)'), log_file


        candidates = azim_df_search_vertex(search_setting, project=project)
        foreach candidate, candidates do begin
            time_range = candidate.time_range
            if keyword_set(test_time) then if product(time_range-test_time) ge 0 then continue else stop
            lprmsg, '', log_file
            lprmsg, 'Processing candidate '+string(candidate.id,format='(I0)')+' ...', log_file
            lprmsg, strjoin(time_string(candidate.time_range),' to ')+' ...', log_file

            df_group = candidate.df_list
            if df_group.length lt probe_min_count then continue

        ;---Find consecutive different probes as a DF group.
            df_group_list = list()
            probe_list = list()
            df_list = list()
            foreach df, df_group, df_id do begin
                probe_pos = probe_list.where(df.probe)
                if probe_pos eq !null then begin
                    ; Add the current probe to list if it doesn't exist.
                    probe_list.add, df.probe
                    df_list.add, df
                endif else begin
                    probe_pos = probe_pos[0]
                    ; The current probe exists in the current list.
                    if probe_list.length lt probe_min_count then begin
                        ; If the current list is not long enough then toss and start over.
                        for iter=0, probe_pos do begin
                            probe_list.remove, 0
                            df_list.remove, 0
                        end
                    endif else begin
                        ; The current list is long enough, pop it.
                        current_list = list()   ; must copy manually
                        for iter=0, df_list.length-1 do current_list.add, df_list[iter]
                        df_group_list.add, current_list
                        ; Remove probes before the current probe.
                        for iter=0, probe_pos do begin
                            probe_list.remove, 0
                            df_list.remove, 0
                        endfor
                    endelse
                    ; Add the new probe to the list.
                    probe_list.add, df.probe
                    df_list.add, df
                endelse
                msg = tab+'current list: '
                foreach tmp, df_list do msg += tmp.probe+' '
                lprmsg, msg, log_file
                if df_id eq df_group.length-1 and probe_list.length ge probe_min_count then begin
                    ; Pop the current list.
                    df_group_list.add, df_list
                    probe_list = list()
                    df_list = list()
                endif
            endforeach

            foreach df_group, df_group_list do begin
                if df_group.length lt probe_min_count then message, 'Inconsistency, stop ...'

                ndf = df_group.length
                ndim = 3
                df_obs_times = dblarr(ndf)
                df_r_sms = fltarr(ndf,ndim)
                probes = strarr(df_group.length)
                foreach df, df_group, ii do begin
                    df_obs_times[ii] =  df.obs_time
                    probes[ii] = df.probe
                    df_r_sms[ii,*] = df.obs_r_sm
                endforeach
                probe_list = list(probes, /extract)

                ; Check dtime.
                dtimes = df_obs_times[1:ndf-1]-df_obs_times[0:ndf-2]
                lprmsg, tab+'edge dtimes (sec): '+strjoin(string(dtimes,format='(I0)'),','), log_file
                if max(dtimes) ge dtime_max then begin
                    lprmsg, 'dtime too large, skip ...', log_file
                    continue
                endif
                mean_dtime = mean(dtimes)
                lprmsg, tab+'edge mean dtime (sec): '+string(dtimes,format='(I0)'), log_file
                if mean_dtime ge dtime_mean_max then begin
                    lprmsg, 'dtime mean too large, skip ...', log_file
                    continue
                endif

                ; Check triad.
                ntriad_vertex = 3
                triad_combos = choose_from(probes, ntriad_vertex)
                triad_list = list()
                re = constant('re')
                km_per_s_2_deg_per_min = 1d/re*60*constant('deg')
                foreach triad_probes, triad_combos do begin
                    sorted_probes = triad_probes[sort(triad_probes)]
                    key = strjoin(sorted_probes,'_')

                    ; Get info about the vertex.
                    probe_index = []
                    foreach probe, triad_probes do probe_index = [probe_index, probe_list.where(probe)]
                    vertex_obs_times = df_obs_times[probe_index]
                    vertex_r_sms = df_r_sms[probe_index,*]

                    ; Sort by obs_time.
                    index = sort(vertex_obs_times)
                    triad_probes = triad_probes[index]
                    vertex_obs_times = vertex_obs_times[index]
                    vertex_r_sms = vertex_r_sms[index,*]
                    center_r_sm = reform(total(vertex_r_sms,1)/ntriad_vertex)

                    ; Angles.
                    triad_angles = triangle_angles(reform(transpose(vertex_r_sms), [1,ntriad_vertex,ndim]))
                    triad_angles = reform(triad_angles)
                    lprmsg, tab+'triad '+key+' angles (deg): '+strjoin(string(triad_angles,format='(I0)'),','), log_file
                    index = lazy_where(triad_angles, '()', triad_angle_range, count=count)
                    if count ne ntriad_vertex then begin
                        lprmsg, 'Not all angles in range, skip ...', log_file
                        continue
                    endif

                    triad_info = dictionary($
                        'probes', triad_probes, $
                        'times', vertex_obs_times, $
                        'r_sms', vertex_r_sms, $
                        'center_r_sm', center_r_sm, $
                        'angles', triad_angles)

                    ; 3-sc timing based on obs_time.
                    obs_times = triad_info.times
                    obs_r_sms = triad_info.r_sms[*,0:1]

                    ; Solve for v_2d using obs_times.
                    rr = transpose(obs_r_sms[1:*,*]-(obs_r_sms[0,*] ## [1,1]))
                    tt = obs_times[1:*]-obs_times[0]
                    index = where(snorm(transpose(rr)) eq 0, count)
                    if count ne 0 then begin
                        lprmsg, 'cannot solve triad, dr = 0, skip ...', log_file
                        continue
                    endif
                    index = where(tt eq 0, count)
                    if count ne 0 then begin
                        lprmsg, 'cannot solve triad, dt = 0, skip ...', log_file
                        continue
                    endif
                    vv = la_linear_equation(rr,tt)
                    vhat = sunitvec(vv)
                    rr_normal = dblarr(ndim-1)
                    for ii=0, ndim-2 do rr_normal[ii] = sdot(rr[*,ii],vhat)
                    fit_result = linfit(tt, rr_normal)
                    vmag = fit_result[1]*re
                    triad_info['vhat_obs_time'] = vhat
                    triad_info['vmag_obs_time'] = vmag
                    center_rxy = [center_r_sm[0:1],0]
                    omega = vec_cross([vhat,0],sunitvec(center_rxy))*vmag/snorm(center_rxy)*km_per_s_2_deg_per_min
                    triad_info['omega_obs_time'] = omega[2]     ; negative is eastward.

                    triad_info['vhat_xcor'] = [0.,0]
                    triad_info['vmag_xcor'] = 0.
                    triad_info['omega_xcor'] = 0.

                    triad_list.add, triad_info
                endforeach
                ntriad = triad_list.length
                lprmsg, tab+'Found '+string(ntriad,format='(I0)')+' triads with good geometry ...', log_file
                if ntriad lt triad_min_count then begin
                    lprmsg, 'Not enough # of triads, skip ...', log_file
                    continue
                endif


                event_id += 1
                subgroup = dictionary($
                    'id', event_id, $
                    'time_range', minmax(df_obs_times), $
                    'search_name', candidate.search_name, $
                    'region', candidate.region, $
                    'probes', probes, $
                    'df_list', df_group, $
                    'triad_list', triad_list, $
                    'edge_list', list())

                azim_df_subgroup_write_file, subgroup, filename=out_file
            endforeach
        endforeach
    endif

    if keyword_set(test_time) then stop
    return, azim_df_subgroup_read_file(out_file)
end
