;+
; Resolve overlapping subgroups.
;-

function azim_df_resolve_overlap_subgroup_analyze_subgroup, candidate, project=project

    ; Get the r2 for linear fit on UT-MLT.
    df_list = candidate.df_list
    ndf = df_list.length
    obs_times = dblarr(ndf)
    foreach df, df_list, ii do obs_times[ii] = df.obs_time

    obs_mlts = fltarr(ndf)
    foreach df, df_list, ii do obs_mlts[ii] = df.obs_mlt

    obs_rxys = fltarr(ndf)
    foreach df, df_list, ii do obs_rxys[ii] = df.obs_rxy


    info = dictionary()

    key = 'mlt_fit'
    xxs = obs_times
    yys = obs_mlts
    fit_result = linfit(xxs, yys, yfit=yfit)
    ss_res = total((yys-yfit)^2)
    ss_tot = total((yys-mean(yys))^2)
    r_square = 1-ss_res/ss_tot
    info[key] = dictionary($
        'slope', fit_result[1], $
        'r_square', r_square)

    key = 'rxy_fit'
    xxs = obs_times
    yys = obs_rxys
    fit_result = linfit(xxs, yys, yfit=yfit)
    ss_res = total((yys-yfit)^2)
    ss_tot = total((yys-mean(yys))^2)
    r_square = 1-ss_res/ss_tot
    info[key] = dictionary($
        'slope', fit_result[1], $
        'r_square', r_square)

    return, info

end

function azim_df_resolve_overlap_subgroup, candidate_list, project=project, $
    settings=settings, log_file=log_file

    property_list = list()
    foreach candidate, candidate_list do begin
        property_list.add, azim_df_resolve_overlap_subgroup_analyze_subgroup(candidate, project=project)
    endforeach


    coefs = list()
    foreach info, property_list do begin
        mlt_r_square = info.mlt_fit.r_square
        rxy_r_square = info.rxy_fit.r_square
        coef = max([mlt_r_square,rxy_r_square])
        coefs.add, coef
    endforeach
    coefs = coefs.toarray()

    max_coef = max(coefs, index)

    return, candidate_list[index[0]]

end



;;---Load all df_group.
;    project = azim_df_load_project()
;    search_step = 'df_group'
;    candidates = list()
;    candidates.add, azim_df_search_post_midn_events(search_step=search_step), /extract
;    candidates.add, azim_df_search_pre_midn_events(search_step=search_step), /extract
;;stop
;;    dirname = 'df_group'
;;    azim_df_subgroup_gen_diagnostic_plot, candidates, dirname=dirname, project=project
;
;;---Sort by time.
;    times = list()
;    foreach candidate, candidates do times.add, candidate.time_range[0]
;    index = sort(times.toarray())
;    candidates = candidates[index]
;    foreach candidate, candidates, id do candidate.id = id


    project = azim_df_load_project()
    search_step = 'df_group'
    candidates = azim_df_search_all_events(project=project, search_step=search_step)


;---Get overlap_list
    overlap_list = list()
    current_list = list()
    isolated_list = list()
    event_list = list()
    foreach the_group, candidates, id do begin
        str_id = string(id,format='(I0)')
        lprmsg, 'current candidate: '+str_id
        if current_list.length ne 0 then begin
            pre_group = current_list[-1]
            pre_time = pre_group.time_range
            the_time = the_group.time_range
            if min(the_time) ge max(pre_time) then begin
                lprmsg, 'no overlap with previous candidate, pop current_list'
                if current_list.length gt 1 then begin
                    lprmsg, 'add current_list to overlap_list'
                    overlap_list.add, current_list
                    event_list.add, current_list
                endif else begin
                    lprmsg, 'add current_list to isolated_list'
                    isolated_list.add, current_list
                    event_list.add, current_list
                endelse
                lprmsg, 'reset current_list'
                current_list = list()
            endif
        endif
        lprmsg, 'add '+str_id+' to current_list'
        current_list.add, the_group
        if the_group.id eq candidates.length-1 then begin
            if current_list.length gt 1 then begin
                lprmsg, 'add current_list to overlap_list'
                overlap_list.add, current_list
                event_list.add, current_list
            endif else begin
                lprmsg, 'add current_list to isolated_list'
                isolated_list.add, current_list
                event_list.add, current_list
            endelse
        endif
    endforeach
    stop

;---Select test_list.
    root_dir = join_path([project.plot_dir,'diagnostic_plot'])
    foreach candidate_list, overlap_list, group_id do begin
        ; Copy file to groups.
        in_dir = join_path([root_dir,'azim_df_orig_df_groups'])
        out_dir = join_path([root_dir,'resolve_orig_df_groups',string(group_id+1,format='(I0)')])
        if file_test(out_dir) eq 0 then file_mkdir, out_dir
        foreach candidate, candidate_list do begin
            file_suffix = 'azim_df_event_'+strjoin(time_string(candidate.time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
            in_file = join_path([in_dir,file_suffix])
            out_file = join_path([out_dir,file_suffix])
            file_copy, in_file, out_file
        endforeach

        ; Copy selected file to dir.
        selected_event = azim_df_resolve_overlap_subgroup(candidate_list, project=project)
        file_suffix = 'azim_df_event_'+strjoin(time_string(selected_event.time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
        out_dir = join_path([root_dir,'resolve_orig_df_groups','selected'])
        if file_test(out_dir) eq 0 then file_mkdir, out_dir
        in_file = join_path([in_dir,file_suffix])
        out_file = join_path([out_dir,file_suffix])
        file_copy, in_file, out_file
    endforeach
end
