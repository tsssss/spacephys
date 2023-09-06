;+
; Check the Pdyn of all the detected DFs.
;-

test = 0
    project = azim_df_load_project()
    events = azim_df_search_df_group_filter(project=project)

    df_times = list()
    foreach event, events do foreach df, event.df_list do begin
        df_times.add, df.arrival_time
    endforeach
    df_times = df_times.toarray()
    
    
    secofday = constant('secofday')
    time_range = minmax(df_times)+[-1,1]*secofday
    p_var = 'p'
    if check_if_update(p_var, time_range=time_range) then omni_read, time_range, id='pdyn'
    
    pdyns = get_var_data(p_var, at=df_times)
    index = where_pro(pdyns, '()', [10,50], count=count)
    the_times = df_times[index]

    stop

end