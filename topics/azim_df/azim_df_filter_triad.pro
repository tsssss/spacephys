;+
; Check if a given candidate pass criteria for triads.
; Adopted from azim_df_filter_df_group.
;-

function azim_df_filter_triad, candidate, project=project, $
    settings=settings, $
    log_file=log_file, $
    ; Output for diagnostic purpose.
    vmag_ratios=vmag_ratios, $
    vhat_angles=vhat_angles

    tab = constant('4space')
    if n_elements(log_file) eq 0 then log_file = -1
    if n_elements(settings) eq 0 then begin
        settings = dictionary($
            'triad_angle_range', [15.,165], $
            'max_angle_diff', 30., $
            'max_vmag_ratio', 0.2, $
            'min_triad_count', 2)
    endif

    triad_angle_range = settings.triad_angle_range
    max_angle_diff = settings.max_angle_diff
    max_vmag_ratio = settings.max_vmag_ratio
    min_triad_count = settings.min_triad_count
    if n_elements(candidate) eq 0 then begin
        lprmsg, 'Settings for filtering triad ...', log_file
        lprmsg, 'Triad angle range (deg): ['+strjoin(string(triad_angle_range,format='(I0)'),',')+']', log_file
        lprmsg, 'Max angle diff (deg):'+string(max_angle_diff,format='(I0)'), log_file
        lprmsg, 'Max vmag diff ratio (#): '+string(max_vmag_ratio,format='(F3.1)'), log_file
        lprmsg, '', log_file
        return, settings
    endif

    retval = dictionary()
    if n_elements(project) eq 0 then project = azim_df_load_project()
    lprmsg, '', log_file
    msg = 'Processing candidate: '+strjoin(time_string(candidate.time_range),' to ')
    lprmsg, msg, log_file


;---Gather info.
    df_list = candidate.df_list
    probe_list = list()
    foreach df, df_list do begin
        probe = df.probe
        if probe_list.where(probe) eq !null then probe_list.add, probe
    endforeach
    nprobe = probe_list.length
    probes = probe_list.toarray()
    ndim = 3
    df_r_sms = fltarr(nprobe,ndim)
    foreach df, df_list, ii do df_r_sms[ii,*] = df.obs_r_sm
    df_obs_times = dblarr(nprobe)
    foreach df, df_list, ii do df_obs_times[ii] = df.obs_time
    df_xcor_times = dblarr(nprobe)
    df_xcor_times[0] = df_list[0].obs_time
    edge_list = candidate.edge_list
    if edge_list.length ne probe_list.length-1 then begin
        lprmsg, 'Inconsistent probe and edge, skip ...', log_file
        return, retval
    endif
    foreach edge, edge_list, ii do df_xcor_times[ii+1] = df_xcor_times[ii]+edge.time_lag_xcor


;---Check triad geometry, timing.
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
        obs_times = df_obs_times[probe_index]
        xcor_times = df_xcor_times[probe_index]
        obs_r_sms = df_r_sms[probe_index,*]

        ; Sort by obs_time.
        index = sort(obs_times)
        triad_probes = triad_probes[index]
        obs_times = obs_times[index]
        xcor_times = xcor_times[index]
        obs_r_sms = obs_r_sms[index,*]
        center_r_sm = reform(total(obs_r_sms,1)/ntriad_vertex)

        triad = dictionary($
            'probes', triad_probes, $
            'times', obs_times, $
            'r_sms', obs_r_sms, $
            'center_r_sm', center_r_sm )


        ; Angles.
        triad_angles = triangle_angles(reform(transpose(obs_r_sms), [1,ntriad_vertex,ndim]))
        triad_angles = reform(triad_angles)
        index = lazy_where(triad_angles, '()', triad_angle_range, count=count)
        lprmsg, tab+'angles (deg): '+strjoin(string(triad_angles,format='(I0)'),' '), log_file
        if count ne ndim then begin
            lprmsg, 'bad geometry, skip ...', log_file
            continue
        endif
        triad['angles'] = triad_angles

        ; Solve for v_2d.
        timing = azim_df_solve_triad_timing(xcor_times, obs_r_sms)
        if n_elements(timing) eq 0 then begin
            lprmsg, 'no solution for xcor_times, skip ...', log_file
            continue
        endif
        triad['vhat_xcor'] = timing.vhat
        triad['vmag_xcor'] = timing.vmag
        triad['omega_xcor'] = timing.omega

        timing = azim_df_solve_triad_timing(obs_times, obs_r_sms)
        if n_elements(timing) eq 0 then begin
            lprmsg, 'no solution for obs_times, skip ...', log_file
            continue
        endif
        triad['vhat_obs_time'] = timing.vhat
        triad['vmag_obs_time'] = timing.vmag
        triad['omega_obs_time'] = timing.omega

        ; Filter by vmag and vhat diff.
        vmag_xcor = triad.vmag_xcor
        vmag_obs_time = triad.vmag_obs_time
        vmag_ratio = abs(vmag_xcor-vmag_obs_time)/vmag_obs_time
        vmag_ratio = round(vmag_ratio*10)/10.
        lprmsg, tab+'vmag_ratio (#): '+string(vmag_ratio,format='(F4.1)'), log_file
        if vmag_ratio gt max_vmag_ratio then begin
            lprmsg, 'vmag diff too large, skip ...', log_file
            continue
        endif

        vhat_xcor = triad.vhat_xcor
        vhat_obs_time = triad.vhat_obs_time
        vhat_angle = sang(vhat_xcor,vhat_obs_time, /deg)
        vhat_angle = round(vhat_angle)
        lprmsg, tab+'vhat angle (deg): '+string(vhat_angle,format='(I0)'), log_file
        if vhat_angle gt max_angle_diff then begin
            lprmsg, 'vhat diff too large, skip ...', log_file
            continue
        endif

        triad_list.add, triad
    endforeach


    lprmsg, tab+'Found '+string(triad_list.length,format='(I0)')+' triads with good timing ...', log_file
    if triad_list.length lt min_triad_count then begin
        lprmsg, tab+'No enough triads with good timing, skip ...', log_file
        return, retval
    endif


;---Update probes, df_list, etc.
    candidate.triad_list = triad_list


    msg = 'Pass ...'
    lprmsg, msg, log_file
    return, candidate
end
