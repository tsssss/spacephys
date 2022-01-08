;+
; Find subgroups for all regions.
;-

function azim_df_find_subgroup, project=project, reset=reset, log_file=log_file

    retval = list()
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(log_file) eq 0 then log_file = -1
    file_suffix = project.name+'_find_subgroup.txt'
    out_file = join_path([project.data_dir,file_suffix])
    if keyword_set(reset) then file_delete, out_file, /allow_nonexistent
    if file_test(out_file) eq 1 then begin
        if file_lines(out_file) eq 0 then file_delete, out_file
    endif

    if file_test(out_file) eq 0 then begin
        ftouch, out_file

        regions = list()
        regions.add, dictionary($
            'name', 'pre_midn', $
            'mlt_range', [-9,0])
        regions.add, dictionary($
            'name', 'post_midn', $
            'mlt_range', [0,9])
        regions.add, dictionary($
            'name', 'around_midn', $
            'mlt_range', [-3,3])


        events = list()
        foreach region, regions do begin
            search_setting = dictionary()
            search_steps = 'search_'+['roi','vertex','subgroup']
            foreach search_step, search_steps do begin
                search_info = dictionary($
                    'file_suffix', 'azim_df_search_'+region.name+'_'+search_step+'.txt')
                if search_step eq 'search_roi' then begin
                    search_info['name'] = region.name
                    search_info['mlt_range'] = region.mlt_range
                endif
                search_setting[search_step] = search_info
            endforeach

            ; Candidates within region.
            candidates = azim_df_search_subgroup(search_setting, project=project)

            ; Resolve overlapping subgroups.
            candidate_list = azim_df_sort_overlap_candidate(candidates)
            foreach group, candidate_list, group_id do begin
                if group.length eq 1 then begin
                    group = group[0]
                endif else begin
                    group = azim_df_resolve_overlap_subgroup(group, project=project)
                endelse
                events.add, group
            endforeach
        endforeach


    ;---Sort by time.
        nevent = events.length
        times = dblarr(nevent)
        foreach event, events, ii do times[ii] = event.time_range[0]
        index = sort(times)
        events = events[index]


    ;---Resolve overlapping over regions.
        subgroups = list()
        event_list = azim_df_sort_overlap_candidate(events)
        foreach group, event_list do begin
            if group.length eq 1 then begin
                group = group[0]
            endif else begin
                ; This is a sanity check.
                ;if group.length eq 2 then continue
                ;foreach tmp, group do print, time_string(tmp.time_range)
                ngroup = group.length
                if ngroup ne 2 then message, 'Inconsistency ...'
                region_names = strarr(ngroup)
                for ii=0,ngroup-1 do region_names[ii] = group[ii].region
                ;lprmsg, region_names

                index = where(region_names eq 'around_midn', complement=index2)
                around_midn_group = group[index[0]]
                other_group = group[index2[0]]

                dtimes = around_midn_group.time_range-other_group.time_range
                if dtimes[0] le 0 and dtimes[1] ge 0 then begin
                    ; other_group encloses around_midn_group
                    group = other_group
                endif else if dtimes[0] ge 0 and dtimes[1] le 0 then begin
                    ; around_midn_group encloses other_group
                    group = around_midn_group
                endif else begin
                    ; merge the two.
                    df_list = list()
                    foreach df, around_midn_group.df_list do df_list.add, df
                    foreach df, other_group.df_list do df_list.add, df
                    df_str = strarr(df_list.length)
                    foreach df, df_list, ii do df_str[ii] = df.probe+'_'+time_string(df.time_range[0])
                    index = uniq(df_str, sort(df_str))
                    df_list = df_list[index]
                    group = other_group
                    group.df_list = df_list

                    ndf = df_list.length
                    probes = strarr(ndf)
                    foreach df, df_list, ii do probes[ii] = df.probe
                    group.probes = probes

                    obs_times = dblarr(ndf)
                    foreach df, df_list, ii do obs_times[ii] = df.obs_time
                    time_range = minmax(obs_times)
                    group.time_range = time_range
                    group.edge_list = list()
                    group.triad_list = list()
                endelse
            endelse
            subgroups.add, group
        endforeach


    ;---Filter by edge and triad.
        edge_flags = bytarr(subgroups.length)+1
        foreach group, subgroups, ii do begin
            group = azim_df_filter_edge(group, project=project, log_file=log_file)
            if n_elements(group) eq 0 then begin
                edge_flags[ii] = 0
                continue
            endif
            group = azim_df_filter_triad(group, project=project, log_file=log_file)
            if n_elements(group) eq 0 then begin
                edge_flags[ii] = 0
                continue
            endif
        endforeach
        index = where(edge_flags eq 1, count)
        if count eq 0 then return, retval
        subgroups = subgroups[index]


    ;---Output.
        foreach subgroup, subgroups, ii do subgroups[ii].id = ii+1
        foreach subgroup, subgroups do azim_df_subgroup_write_file, subgroup, filename=out_file
    endif

    subgroups = azim_df_subgroup_read_file(out_file)
    ndf = 0
    foreach subgroup, subgroups do ndf += subgroup.df_list.length
    lprmsg, '# of subgroups: '+string(subgroups.length,format='(I0)')
    lprmsg, '# of DFs: '+string(ndf,format='(I0)')

    return, subgroups

end

reset = 0
;log_file = join_path([homedir(),'azim_df_find_subgroup.log'])
;file_delete, log_file, /allow_nonexistent
;ftouch, log_file
candidates = azim_df_find_subgroup(project=project, reset=reset, log_file=log_file)

stop


;---Add test times.
test_event_times = time_double([$
    '2007-11-20/17:18:10', $
    '2008-01-09/11:27:45', $
    '2008-01-19/12:02:55', $
    '2008-02-29/08:26:50', $
    '2014-08-28/10:10:40', $
    '2014-12-26/01:05:25', $
    '2016-10-13/12:22:35', $
    '2016-12-11/09:46:35', $
    '2017-03-28/03:00:40'])
index = list()
foreach test_event_time, test_event_times do begin
    lprmsg, 'Processing event: '+time_string(test_event_time)
    foreach candidate, candidates, id do begin
        if test_event_time eq candidate.time_range[0] then begin
            index.add, id
            lprmsg, 'Found the test event ...'
            break
        endif
    endforeach
endforeach



;dirname = 'azim_df_find_subgroup'
;azim_df_subgroup_gen_diagnostic_plot, events, dirname=dirname, project=project

end