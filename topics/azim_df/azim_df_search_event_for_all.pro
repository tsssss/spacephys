
reset = 0
search_step = 'subgroup'
events = azim_df_search_around_midn_events(reset=reset, search_step=search_step)
events = azim_df_search_pre_midn_events(reset=reset, search_step=search_step)
events = azim_df_search_post_midn_events(reset=reset, search_step=search_step)

;file = '/Volumes/GoogleDrive/My Drive/works/works/azim_df/data/azim_df_analyze_subgroup_good_ones.txt'
;events = azim_df_subgroup_read_file(file)
;events = azim_df_subgroup_filter_by_dtime(events)
;the_events = list()
;min_dt = 30*60  ; sec.
;project = azim_df_load_project()
;foreach event, events do begin
;    df_times = list()
;    foreach df, event.df_list do df_times.add, df.arrival_time
;    df_times = sort_uniq(df_times.toarray())
;    dtimes = df_times[1:*]-df_times[0:-2]
;    index = where(dtimes ge min_dt, count)
;    if count ne 0 then continue
;    
;    event = azim_df_analyze_subgroup(event, project=project)
;    if n_elements(event) eq 0 then continue
;    the_events.add, event
;endforeach
;dirname = '2020_04_subgroups_dtime'
;azim_df_subgroup_gen_diagnostic_plot, the_events, dirname=dirname, project=project

end