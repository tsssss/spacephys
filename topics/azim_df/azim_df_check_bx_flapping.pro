;+
; Check Bx to see if a DF event is a flapping event.
;-

test = 1
    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_load_coherent_events(project=project)

    pad_time = 3600.
    secofday = 86400.
    foreach event, events do begin
        time_range = event.time_range+[1,-1]*secofday
        data_file = join_path([project.data_dir,'event_data',$
            'azim_df_event_data_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmmss'),'_')+'_v01.tplot'])
        tplot_restore, filename=data_file
        
       
        timing_info = event.timing_info.df_time
        probes = timing_info.mlt_sorted_probes
        width = 600.
        foreach probe, probes do begin
            get_data, probe+'_b_sm', times, b_sm
            store_data, probe+'_bmag', times, snorm(b_sm), limits={$
                ytitle: '(nT)', $
                labels: strupcase(probe)}
                
            get_data, probe+'_theta', times, theta
            data_rate = times[1]-times[0]
            theta -= smooth(theta, width/data_rate)
            store_data, probe+'_theta', times, theta, limits={$
                labels: strupcase(probe)}
        endforeach
        
        
        vars = probes+'_bmag'
        options, vars, 'constant', 0
        
        vars = [probes+'_bmag',probes+'_theta']
        tplot, vars, trange=time_range

        stop
    endforeach
end