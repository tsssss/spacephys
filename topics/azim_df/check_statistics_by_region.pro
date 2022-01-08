;+
; Print statistics of the regions: pre-, post-, around-midn.
;-

tab = constant('4space')
regions = ['post','pre']+'_midn'


;coherent_dfs = list()
;foreach region, regions do begin
;    routine = 'azim_df_search_'+region+'_events'
;    candidates = call_function(routine, search_step='coherent_df')
;    dfs = list()
;    foreach candidate, candidates do begin
;        dfs.add, candidate.df_list, /extract
;    endforeach
;    df_strs = strarr(dfs.length)
;    foreach df, dfs, ii do df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
;    uniq_df_strs = sort_uniq(df_strs)
;    lprmsg, region
;    lprmsg, '# of coherent DFs: '+string(n_elements(uniq_df_strs),format='(I0)')
;    coherent_dfs.add, dfs, /extract
;endforeach
;coherent_df_strs = strarr(coherent_dfs.length)
;foreach df, coherent_dfs, ii do coherent_df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
;coherent_df_strs = sort_uniq(coherent_df_strs)
;lprmsg, '# of coherent DFs: '+string(n_elements(coherent_df_strs),format='(I0)')
;stop



;    df_dfs = list()
;    foreach region, regions do begin
;        routine = 'azim_df_search_'+region+'_events'
;        candidates = call_function(routine, search_step='df')
;        dfs = list()
;        foreach candidate, candidates do begin
;            dfs.add, candidate.df_list, /extract
;        endforeach
;        df_strs = strarr(dfs.length)
;        df_mlts = fltarr(dfs.length)
;        df_rxys = fltarr(dfs.length)
;        foreach df, dfs, ii do begin
;            df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
;            df_mlts[ii] = df.obs_mlt
;            df_rxys[ii] = df.obs_rxy
;        endforeach
;        uniq_df_strs = sort_uniq(df_strs)
;        lprmsg, region
;        lprmsg, '# of DF group DFs: '+string(n_elements(uniq_df_strs),format='(I0)')
;        lprmsg, 'MLT range (hr): ['+strjoin(string(minmax(df_mlts),format='(F5.1)'),',')+']'
;        lprmsg, 'Rxy range (Re): ['+strjoin(string(minmax(df_rxys),format='(F5.1)'),',')+']'
;        df_dfs.add, dfs, /extract
;    endforeach
;    df_df_strs = strarr(df_dfs.length)
;    foreach df, df_dfs, ii do df_df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
;    uniq_df_df_strs = sort_uniq(df_df_strs)
;    lprmsg, '# of DF group DFs: '+string(n_elements(uniq_df_df_strs),format='(I0)')
;    stop



    large_dfs = list()
    foreach region, regions do begin
        lprmsg, region
        routine = 'azim_df_search_'+region+'_events'
    
        msg = '# of candidate periods: '
        search_step = 'triad'
        candidates = call_function(routine, search_step=search_step)
        msg += string(candidates.length,format='(I0)')
        lprmsg, tab+msg
    
        msg = 'Duration of candidate periods (hr): '
        duration = 0
        foreach candidate, candidates do duration += total(candidate.time_range*[-1,1])
        msg += sgnum2str(duration/3600.,ndec=1)
        lprmsg, tab+msg
    
    
        msg = '# of large DFs: '
        search_step = 'large_df'
        candidates = call_function(routine, search_step=search_step)
        nlarge_df = 0
        foreach candidate, candidates do begin
            nlarge_df += candidate.df_list.length
            large_dfs.add, candidate.df_list, /extract
        endforeach
        msg += string(nlarge_df,format='(I0)')
        lprmsg, tab+msg
    endforeach
    
    large_df_strs = strarr(large_dfs.length)
    foreach df, large_dfs, ii do large_df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
    uniq_large_df_strs = sort_uniq(large_df_strs)
    lprmsg, '# of large DFs: '+string(n_elements(uniq_large_df_strs),format='(I0)')
stop


subgroup_dfs = list()
foreach region, regions do begin
    routine = 'azim_df_search_'+region+'_events'
    candidates = call_function(routine, search_step='subgroup')
    dfs = list()
    foreach candidate, candidates do begin
        dfs.add, candidate.df_list, /extract
    endforeach
    df_strs = strarr(dfs.length)
    foreach df, dfs, ii do df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
    uniq_df_strs = sort_uniq(df_strs)
    lprmsg, region
    lprmsg, '# of subgroup DFs: '+string(n_elements(uniq_df_strs),format='(I0)')
    subgroup_dfs.add, dfs, /extract
endforeach
subgroup_df_strs = strarr(subgroup_dfs.length)
foreach df, subgroup_dfs, ii do subgroup_df_strs[ii] = df.probe+'_'+time_string(df.obs_time)
uniq_subgroup_df_strs = sort_uniq(subgroup_df_strs)
lprmsg, '# of subgroup DFs: '+string(n_elements(uniq_subgroup_df_strs),format='(I0)')
stop


end
