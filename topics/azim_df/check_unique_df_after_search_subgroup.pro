;+
; Check how many unique vertex are found in each search and in total.
;-


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
    foreach df, all_df_list, ii do df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
    uniq_df_strs = sort_uniq(df_strs)
    ndf = n_elements(uniq_df_strs)
    lprmsg, 'All uniq DF: '+string(ndf,format='(I0)')
end
