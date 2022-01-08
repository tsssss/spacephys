;+
; Analyze time lag diff between df_time and xcor_time
;-

if n_elements(project) eq 0 then project = azim_df_load_project()
;events = azim_df_preprocess_df_subgroup(project=project)
events = azim_df_subgroup_analyze_overlap(project=project)
;events = azim_df_filter_triads(events)
;events = azim_df_filter_v2d(events)


all_dfs = list()
df_times = list()
foreach subgroup, events do begin
    df_list = subgroup.df_list
    foreach df, df_list do begin
        all_dfs.add, df
        df_times.add, df.arrival_time
    endforeach
endforeach

df_times = df_times.toarray()
index = uniq(df_times, sort(df_times))
uniq_dfs = all_dfs[index]


df_times = list()
xcor_times = list()
time_lag_df = list()
time_lag_xcor = list()
foreach subgroup, events do begin
    df_list = subgroup.df_list
    foreach df, df_list, ii do begin
        if ii eq 0 then continue
        all_dfs.add, df
        df_times.add, df.arrival_time
        xcor_times.add, df.xcor_time
        
        df_pre = df_list[ii-1]
        time_lag_df.add, df.arrival_time-df_pre.arrival_time
        time_lag_xcor.add, df.xcor_time-df_pre.arrival_time
    endforeach
endforeach
df_times = df_times.toarray()
xcor_times = xcor_times.toarray()
time_lag_df = time_lag_df.toarray()
time_lag_xcor = time_lag_xcor.toarray()


index = uniq(df_times, sort(df_times))
df_times = df_times[index]
xcor_times = xcor_times[index]
time_lag_df = time_lag_df[index]
time_lag_xcor = time_lag_xcor[index]

time_lag_diff = time_lag_df-time_lag_xcor
xcor_yys = histogram(time_lag_diff, binsize=10, locations=xcor_xxs)
index = where(xcor_xxs ne 0)
xcor_yys = interpol(xcor_yys[index], xcor_xxs[index], xcor_xxs)
yys = gaussfit(xcor_xxs, smooth(xcor_yys,3), nterm=3, coef)
plot, xcor_xxs,  smooth(xcor_yys,3), psym=1
oplot, xcor_xxs, yys, color=sgcolor('red')
hmfw = 2*sqrt(2*alog(2))*coef[2]
print, hmfw
;xyouts, tx,ty,/normal, 'FWHM='+sgnum2str(hmfw,ndec=0)+' sec'



df_times = list()
xcor_times = list()
foreach subgroup, events do begin
    df_list = subgroup.df_list
    foreach df, df_list do begin
        df_times.add, df.arrival_time
        xcor_times.add, df.xcor_time
    endforeach
endforeach

df_times = df_times.toarray()
xcor_times = xcor_times.toarray()
time_lag_diff = df_times-xcor_times
xcor_yys = histogram(time_lag_diff, binsize=10, locations=xcor_xxs)
index = where(xcor_xxs ne 0)
xcor_yys = interpol(xcor_yys[index], xcor_xxs[index], xcor_xxs)
yys = gaussfit(xcor_xxs, smooth(xcor_yys,3), nterm=3, coef)
plot, xcor_xxs,  smooth(xcor_yys,3), psym=1
oplot, xcor_xxs, yys, color=sgcolor('red')
hmfw = 2*sqrt(2*alog(2))*coef[2]
print, hmfw
stop



df_time_lags = list()
xcor_time_lags = list()
foreach subgroup, events do begin
    df_list = subgroup.df_list
    probes = list()
    foreach df, df_list do probes.add, df.probe
    probes = probes.toarray()
    combos = choose_from(probes, 2)
    foreach combo, combos do begin
        df1 = df_list[(where(probes eq combo[0]))[0]]
        df2 = df_list[(where(probes eq combo[1]))[0]]
        df_time_lags.add, df1.arrival_time-df2.arrival_time
        xcor_time_lags.add, df1.xcor_time-df2.xcor_time
    endforeach
endforeach

df_time_lags = df_time_lags.toarray()
xcor_time_lags = xcor_time_lags.toarray()
time_lag_diff = df_time_lags-xcor_time_lags
xcor_yys = histogram(time_lag_diff, binsize=10, locations=xcor_xxs)
index = where(xcor_xxs ne 0)
xcor_yys = interpol(xcor_yys[index], xcor_xxs[index], xcor_xxs)
yys = gaussfit(xcor_xxs, smooth(xcor_yys,3), nterm=3, coef)
plot, xcor_xxs,  smooth(xcor_yys,3), psym=1
oplot, xcor_xxs, yys, color=sgcolor('red')
hmfw = 2*sqrt(2*alog(2))*coef[2]
print, hmfw



end