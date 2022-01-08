
project = azim_df_load_project()
routines = 'azim_df_search_'+['around_midn','post_midn','pre_midn']+'_events'
search_names = ['beyond_15Re','within_15Re']
search_step = 'triad'
foreach routine, routines do begin
    candidates = call_function(routine, search_step=search_step, project=project)
    foreach search_name, search_names do begin
        count = 0
        hours = list()
        foreach candidate, candidates do begin
            if candidate.search_type ne search_name then continue
            count += 1
            hours.add, total(candidate.time_range*[-1,1])/3600
        endforeach
        if n_elements(hours) eq 0 then continue
        hours = hours.toarray()
        lprmsg, routine+', '+search_name
        lprmsg, 'count: '+string(count)
        lprmsg, 'hours: '+string(total(hours))
        lprmsg, 'mean (hr): '+string(mean(hours))
    endforeach
endforeach

end
