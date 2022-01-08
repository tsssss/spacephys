;+
; Check the error in edge timing.
;-

project = azim_df_load_project()
search_step = 'uniq_subgroup'
candidates = azim_df_search_all_events(project=project, search_step=search_step)

time_lag_dfs = list()
time_lag_xcors = list()
time_lag_diffs = list()
time_lag_ratios = list()
foreach candidate, candidates do begin
    foreach edge, candidate.edge_list do begin
        time_lag_dfs.add, edge.time_lag_df
        time_lag_xcors.add, edge.time_lag_xcor
        time_lag_diffs.add, edge.time_lag_diff
        time_lag_ratios.add, edge.time_lag_ratio
    endforeach
endforeach
time_lag_dfs = time_lag_dfs.toarray()
time_lag_xcors = time_lag_xcors.toarray()
time_lag_diffs = time_lag_diffs.toarray()
time_lag_ratios = time_lag_ratios.toarray()

yys = histogram(time_lag_xcors-time_lag_dfs, location=xxs, binsize=20)
yys = float(yys)
index = where(xxs ne 0)
yys = interpol(yys[index],xxs[index], xxs)
smooth_yys = smooth(yys,3)
fit_yys = gaussfit(xxs, smooth_yys, nterm=3, coef)
plot, xxs, smooth_yys, psym=1
oplot, xxs, fit_yys, color=sgcolor('red')
hmfw = 2*sqrt(2*alog(2))*coef[2]
print, hmfw


end