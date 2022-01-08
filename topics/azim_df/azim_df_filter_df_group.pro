;+
; Check if a given candidate is a DF group.
;-

function azim_df_filter_df_group, candidate, project=project, $
    log_file=log_file, $
    dtime_min_mean=dtime_min_mean, $
    dtime_range=dtime_range, $
    min_abs_mlt=min_abs_mlt, $
    min_mlt_count=min_mlt_count, $
    triad_angle_range=triad_angle_range, $
    triad_dtime_range=triad_dtime_range, $
    min_triad_count=min_triad_count, $
    test_time=test_time, $
    _extra=extra

    retval = dictionary()
    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()


;---Settings.
    dtime_range = [0.,30]*60
    dtime_min_mean = 8.*60
    min_abs_mlt = 4.
    min_mlt_count = 1
    triad_angle_range = [15.,165]
    triad_dtime_range = [1.,1e9]*60
    min_triad_count = 2


;---Load data. Do not need to load data in this step.
    ;if ~keyword_set(skip_load_data) then azim_df_load_basic_data, project=project, scale_width=scale_width


;---Check dtimes.
    dtime_max = dtime_range[1]
    dtime_min = 0.
    df_list = candidate.df_list
    obs_times = dblarr(df_list.length)
    foreach df, df_list, df_id do obs_times[df_id] = df.obs_time
    dtimes = obs_times[1:df_id-1]-obs_times[0:df_id-2]
    msg = tab+'dtimes (min): '+strjoin(string(dtimes,format='(F6.1)'),', ')
    lprmsg, msg, log_file
    ; maximum dtime.
    index = where(dtimes ge dtime_max, count)
    if count ne 0 then begin
        msg = tab+'dtime too large, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif
    ; mean dtime.
    index = where(dtimes gt dtime_min, count)
    msg = tab+'# >'+string(dtime_min,format='(I0)')+' sec: '+string(count,format='(I0)')
    mean_dtime = mean(dtimes[index])
    msg = tab+'dtimes_mean (min): '+string(mean_dtime,format='(F6.1)')
    if mean_dtime ge dtime_min_mean then begin
        msg = tab+'dtime mean too large, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif


;---Check distribution.
    mlts = intarr(df_list.length)
    foreach df, df_list, df_id do mlts[df_id] = round(df.obs_mlt)
    msg = tab+'mlts (hr): '+strjoin(string(mlts,format='(I0)'),', ')
    lprmsg, msg, log_file
    ; nightside.
    index = where(abs(mlts) le min_abs_mlt, count)
    if count lt min_mlt_count then begin
        msg = tab+'no DF in nightside, skip ...'
        lprmsg, msg, log_file
        return, retval
    endif


;---Gather info.
    probe_list = list()
    foreach df, df_list do begin
        probe = df.probe
        if probe_list.where(probe) eq !null then probe_list.add, probe
    endforeach
    nprobe = probe_list.length
    ndim = 3
    df_r_sms = fltarr(nprobe,ndim)
    foreach df, df_list, ii do df_r_sms[ii,*] = df.obs_r_sm
    df_obs_times = dblarr(nprobe)
    foreach df, df_list, ii do df_obs_times[ii] = df.obs_time
    df_widths = fltarr(nprobe)
    foreach df, df_list, ii do df_widths[ii] = df.width


;---Check triad geometry.
    ntriad_vertex = 3
    triad_combos = choose_from(probe_list.toarray(), ntriad_vertex)
    triad_list = dictionary()
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

        triad_list[key] = triad_info
    endforeach
    lprmsg, tab+'Found '+string(triad_list.count(),format='(I0)')+' triads ...', log_file


;---Filter by geometry.
    foreach key, triad_list.keys() do begin
        triad_info = triad_list[key]
        lprmsg, tab+'triad: '+key, log_file

        angles = triad_info.angles
        lprmsg, tab+tab+'angles (deg): '+strjoin(string(angles,format='(F5.1)'),','), log_file
        index = lazy_where(angles, '()', triad_angle_range, count=count)
        if count lt ntriad_vertex then begin
            lprmsg, tab+'Not all angles in range, skip ...', log_file
            triad_list.remove, key
            continue
        endif

        dtimes = triad_info.times
        dtimes = abs(dtimes-shift(dtimes,1))
        lprmsg, tab+tab+'dtimes (sec): '+strjoin(string(dtimes,format='(I0)'),','), log_file
        index = lazy_where(dtimes, '[]', triad_dtime_range, count=count)
        if count lt ntriad_vertex then begin
            lprmsg, tab+'Not all dtimes in range, skip ...', log_file
            triad_list.remove, key
            continue
        endif
    endforeach
    lprmsg, tab+'Found '+string(triad_list.count(),format='(I0)')+' triads with good geometry ...', log_file
    if triad_list.count() lt min_triad_count then begin
        lprmsg, tab+'No enough triads with good geometry, skip ...', log_file
        return, retval
    endif else candidate['triad_list'] = triad_list


;---Update probes, df_list.
    probe_list = list()
    foreach triad, triad_list do probe_list.add, triad.probes, /extract
    probes = sort_uniq(probe_list.toarray())
    probe_list = list(probes, /extract)
    flags = bytarr(df_list.length)
    foreach df, df_list, df_id do if probe_list.where(df.probe) eq !null then flags[df_id] = 1
    index = where(flags eq 0, ndf)
    if ndf eq 0 then message, 'Inconsistency ...'
    df_list = df_list[index]

    candidate.probes = probes
    candidate.df_list = df_list

    msg = 'Pass ...'
    lprmsg, msg, log_file
    return, candidate
end


search_step = 'subgroup'
candidates = list()
candidates.add, azim_df_search_post_midn_events(search_step=search_step), /extract
candidates.add, azim_df_search_pre_midn_events(search_step=search_step), /extract
events = list()
selected_index = list()
rejected_index = list()
log_file = -1
lprmsg, 'Start from '+string(candidates.length,format='(I0)')+' ...', log_file
foreach candidate, candidates, candidate_id do begin
    lprmsg, 'Processing '+string(candidate_id,format='(I0)')+' ...', log_file
    event = azim_df_filter_df_group(candidate, project=project, skip_load_data=1, log_file=log_file)
    if n_elements(event) eq 0 then begin
        rejected_index.add, candidate_id
    endif else begin
        events.add, event
        selected_index.add, candidate_id
    endelse
endforeach
selected_index = selected_index.toarray()
rejected_index = rejected_index.toarray()



subgroup_df_strs = list()
foreach candidate, candidates[selected_index] do foreach df, candidate.df_list do subgroup_df_strs.add, df.probe+'_'+time_string(df.obs_time)
uniq_subgroup_df_strs = sort_uniq(subgroup_df_strs.toarray())
lprmsg, '# of DFs: '+string(n_elements(uniq_subgroup_df_strs),format='(I0)')
stop

root_dir = join_path([project.plot_dir,'diagnostic_plot'])
in_dir = join_path([root_dir,'all_subgroups'])
plot_files = file_search(join_path([in_dir,'*.pdf']))
plot_base_files = fgetbase(plot_files)
if n_elements(plot_files) ne candidates.length then stop

foreach the_dir, ['selected','rejected'] do begin
    the_index = (the_dir eq 'selected')? selected_index: rejected_index
    dirname = join_path([root_dir,'filter_df_group',the_dir])
    if file_test(dirname,/directory) eq 1 then file_delete, dirname, /recursive
    file_mkdir, dirname
    foreach candidate, candidates[the_index], id do begin
        file_suffix = 'azim_df_event_'+strjoin(time_string(candidate.time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
        index = where(plot_base_files eq file_suffix, count)
        if count eq 0 then stop
        out_file = join_path([dirname,file_suffix])
        in_file = join_path([in_dir,file_suffix])
        file_copy, in_file, out_file
    endforeach
endforeach


end
