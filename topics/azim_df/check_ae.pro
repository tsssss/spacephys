;+
; Check the AE level of events, candidates, etc.
;-


omni_read_index, time_double(['2007-03-01','2017-10-01'])

;---All DF groups.
;    events = azim_df_find_dfgroup()
;    the_key = 'max_ae'
;    foreach event, events do begin
;        if event.haskey(the_key) then continue
;        omni_read_index, event.time_range
;        event[the_key] = max(get_var_data('ae', in=event.time_range))
;    endforeach
;
;    nevent = events.length
;    aes = fltarr(nevent)
;    foreach event, events, ii do aes[ii] = event.max_ae


;---Filter into subgroups.
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


    ncandidate = list()
    ndfs = list()
    all_df_list = list()
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
        df_list = list()
        candidates = azim_df_search_subgroup(search_setting, project=project)
        foreach candidate, candidates do df_list.add, candidate.df_list, /extract

        ndf = df_list.length
        df_strs = strarr(ndf)
        foreach df, df_list, ii do df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
        uniq_df_strs = sort_uniq(df_strs)
        ndf = n_elements(uniq_df_strs)
        lprmsg, region.name+' uniq DF: '+string(ndf,format='(I0)')
        lprmsg, 'subgroups: '+string(candidates.length,format='(I0)')

        ndfs.add, ndf
        all_df_list.add, df_list, /extract
    endforeach


    ndf = all_df_list.length
    df_strs = strarr(ndf)
    foreach df, all_df_list, ii do df_strs[ii] = df.probe+'_'+strjoin(time_string(df.time_range),'_')
    uniq_df_strs = sort_uniq(df_strs)
    ndf = n_elements(uniq_df_strs)
    lprmsg, 'All uniq DF: '+string(ndf,format='(I0)')

    aes = fltarr(ndf)
    foreach uniq_df_str, uniq_df_strs, ii do begin
        infos = strsplit(uniq_df_str,'_',/extract)
        time_range = time_double(infos[1:2])
        aes[ii] = max(get_var_data('ae', in=time_range))
    endforeach
    print, mean(aes)
    index = where(aes ge 300, count)
    print, count
    stop








;---All DF.
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


    ndfs = list()
    all_df_list = list()
    foreach region, regions do begin
        search_setting = dictionary()
        search_steps = 'search_'+['roi','vertex']
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
        df_list = list()
        candidates = azim_df_search_vertex(search_setting, project=project)
        foreach candidate, candidates do df_list.add, candidate.df_list, /extract

        ndf = df_list.length
        df_strs = strarr(ndf)
        foreach df, df_list, ii do df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
        uniq_df_strs = sort_uniq(df_strs)
        ndf = n_elements(uniq_df_strs)
        lprmsg, region.name+' uniq DF: '+string(ndf,format='(I0)')

        ndfs.add, ndf
        all_df_list.add, df_list, /extract
    endforeach


    ndf = all_df_list.length
    df_strs = strarr(ndf)
    foreach df, all_df_list, ii do df_strs[ii] = df.probe+'_'+strjoin(time_string(df.time_range),'_')
    uniq_df_strs = sort_uniq(df_strs)
    ndf = n_elements(uniq_df_strs)
    lprmsg, 'All uniq DF: '+string(ndf,format='(I0)')

    aes = fltarr(ndf)
    foreach uniq_df_str, uniq_df_strs, ii do begin
        infos = strsplit(uniq_df_str,'_',/extract)
        time_range = time_double(infos[1:2])
        aes[ii] = max(get_var_data('ae', in=time_range))
    endforeach
    print, mean(aes)
    index = where(aes ge 300, count)
    print, count
    stop
end
