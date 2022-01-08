;+
; Check how many DFs are excluded after forming subgroups.
;-


search_types = ['around_midn','post_midn','pre_midn']


ndfs = list()
foreach search_type, search_types do begin
    ndf = 0
    routine = 'azim_df_search_'+search_type+'_events'
    candidates = call_function(routine, search_step='df')
    foreach candidate, candidates do ndf += candidate.df_list.length
    ndfs.add, ndf
endforeach
ndfs = ndfs.toarray()


nlarge_dfs = list()
foreach search_type, search_types do begin
    nlarge_df = 0
    routine = 'azim_df_search_'+search_type+'_events'
    candidates = call_function(routine, search_step='large_df')
    foreach candidate, candidates do nlarge_df += candidate.df_list.length
    nlarge_dfs.add, nlarge_df
endforeach
nlarge_dfs = nlarge_dfs.toarray()

nsubgroup_dfs = list()
foreach search_type, search_types do begin
    subgroup_dfs = list()
    routine = 'azim_df_search_'+search_type+'_events'
    candidates = call_function(routine, search_step='subgroup')
    foreach candidate, candidates do begin
        foreach df, candidate.df_list do begin
            subgroup_dfs.add, df.probe+'_'+time_string(df.arrival_time)
        endforeach
    endforeach
    subgroup_dfs = sort_uniq(subgroup_dfs.toarray())
    nsubgroup_df = n_elements(subgroup_dfs)
    nsubgroup_dfs.add, nsubgroup_df
endforeach
nsubgroup_dfs = nsubgroup_dfs.toarray()


foreach search_type, search_types, ii do begin
    nlarge_df = nlarge_dfs[ii]
    nsubgroup_df = nsubgroup_dfs[ii]
    print, nlarge_df, nsubgroup_df
    print, nsubgroup_df/float(nlarge_df)
endforeach

nlarge_df = total(nlarge_dfs)
nsubgroup_df = total(nsubgroup_dfs)
print, nlarge_df, nsubgroup_df
print, nsubgroup_df/float(nlarge_df)



end
