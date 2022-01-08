;+
; Analysis overlapping events.
;-

    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_preprocess_df_subgroup(project=project)
    foreach event, events, ii do events[ii] = azim_df_analysis_event(event, project=project)
    

    
    
    
;    test_times = time_double([$
;        '2007-11-20/17:02', $
;        '2007-12-04/11:45', $
;        '2008-02-29/08:30', $
;        
;        ; Some good examples.
;        '2014-08-28/10:20', $
;        '2016-10-13/12:30', $
;        '2016-11-27/17:05', $
;        '2016-12-11/09:55'])
;    ntest_time = n_elements(test_times)
;;ntest_time = 0
;    if ntest_time ne 0 then begin
;        event_index = intarr(ntest_time)
;        for ii=0, ntest_time-1 do begin
;            foreach event, events do begin
;                time_range = event.time_range+[1,-1]*constant('secofday')
;                if test_time ne time_range[0] then continue
;                event_index[ii] = event.id-1
;                break
;            endforeach
;        endfor
;        events = events[event_index]
;    endif
    

;    ; Check widths of all DFs.
;    df_times = list()
;    df_widths = list()
;    foreach event, events do begin
;        foreach df, event.df_list do begin
;            df_times.add, df.arrival_time
;            df_widths.add, df.width
;        endforeach
;    endforeach
;    df_times = df_times.toarray()
;    df_widths = df_widths.toarray()
;    index = uniq(df_times, sort(df_times))
;    df_times = df_times[index]
;    df_widths = df_widths[index]
;    ; 70 sec to 2390 sec.
;    ; 90-850 sec for some good examples.
;    stop
    

;   Remove events with long DF duration.
;    event_index = bytarr(events.length)
;    max_df_duration = 60*20
;    foreach event, events, ii do begin
;        foreach df, event.df_list do begin
;            ;if df.arrival_time eq time_double('2007-11-30/12:18') then stop
;            if df.width ge max_df_duration then begin
;                event_index[ii] = 1
;                break
;            endif
;        endforeach
;    endforeach
;    events = events[where(event_index eq 0)]
;    print, events.length

;    events = azim_df_filter_triads(events, project=project)
;    print, events.length
;    foreach event, events do print, time_string(event.time_range)
;    events = azim_df_filter_v2d(events, project=project)
;    print, events.length


;    df_duration = list()
;    test_time = time_double('2007-12-04/12:42')
;    foreach event, events do begin
;        if product(test_time-event.time_range) le 0 then stop
;        foreach df, event.df_list do begin
;            df_duration.add, df.width
;        endforeach
;    endforeach
;    df_duration = df_duration.toarray()
;    stop
    
    
;---Divide them into isolated and overlap events.
    foreach event, events, event_id do events[event_id].id = event_id+1
    overlap_list = list()
    current_list = list()
    foreach the_group, events do begin
        if current_list.length ne 0 then begin
            pre_group = current_list[-1]
            pre_time = []
            foreach df, pre_group.df_list do pre_time = [pre_time,df.arrival_time]
            
            the_time = []
            foreach df, the_group.df_list do the_time = [the_time,df.arrival_time]
            
            if min(the_time) gt max(pre_time) then begin
                if current_list.length gt 1 then overlap_list.add, current_list
                current_list = list()
            endif
        endif
        current_list.add, the_group
    endforeach
    
    nevent = n_elements(events)
    isolated_flag = intarr(nevent)
    foreach event_list, overlap_list do begin
        foreach event, event_list do isolated_flag[event.id-1] = 1
    endforeach
    index = where(isolated_flag eq 0, count)
    isolated_list = events[index]
    
    
; Generate survey plot.
;    dirname = ['event_clusters','isolated_events']
;    azim_df_search_df_group_gen_diagnostic_plot, isolated_list, dirname=dirname, project=project
;    foreach event_list, overlap_list, list_id do begin
;        dirname = ['event_clusters','group_'+string(list_id+1,format='(I0)')]
;        azim_df_search_df_group_gen_diagnostic_plot, event_list, dirname=dirname, project=project
;    endforeach




;---Combine isolated events and overlapping events.
    cluster_list = list()
    cluster_list.add, overlap_list, /extract
    foreach event, isolated_list do cluster_list.add, list(event)
    
    ;test_events = [1,12,18,21,23,30,33,35,45,53,64,68,70,71,72,73,77]-1
    ;test_events = [12,35]-1
    cluster_list = cluster_list[test_events]
    
    
    
;---Analysis overlap events.
    secofday = constant('secofday')
    event_list = list()
    foreach cluster, cluster_list, cluster_id do begin
        ;if cluster_id ne 52 then continue else stop
    
        ; All probes.
        df_probes = list()
        foreach event, cluster do df_probes.add, event.probes, /extract
        df_probes = sort_uniq(df_probes.toarray())
        
        ; Full time range.
        df_time_ranges = list()
        foreach event, cluster do df_time_ranges.add, event.time_range+[1,-1]*secofday
        full_time_range = minmax(df_time_ranges.toarray())
        data_time_range = full_time_range+[-1,1]*secofday
        plot_time_range = full_time_range+[-1,1]*600
        
        ; Load data.
        data_vars = ['theta','r_sm','b_sm']
        time_step = project.time_step
        probe_infos = project.probe_infos
        foreach probe, df_probes do begin
            lprmsg, 'Processing probe: '+probe+' ...'
            prefix = probe_infos[probe].prefix
            foreach var, data_vars do begin
                azim_df_read_data, var, time_range=data_time_range, probe=probe, project=project
                the_var = prefix+var
                uniform_time, the_var, time_step
            endforeach
        endforeach
        
        
        
        ; Check scatter of DF times.
        df_times = list()
        foreach event, cluster do begin
            times = list()
            foreach df, event.df_list do times.add, df.arrival_time
            df_times.add, times.toarray()
        endforeach
        df_time_stddev = fltarr(cluster.length)
        df_time_mean = fltarr(cluster.length)
        foreach times, df_times, ii do begin
            dtimes = times[1:-1]-times[0:-2]
            df_time_stddev[ii] = stddev(dtimes)            
            df_time_mean[ii] = mean(dtimes)
        endforeach
        df_time_scatters = df_time_stddev/df_time_mean
        
        
        ; Check the shape of DFs.
        df_widths = list()
        df_heights = list()
        foreach event, cluster do begin
            widths = list()
            heights = list()
            foreach df, event.df_list do begin
                widths.add, df.width
                heights.add, df.height
            endforeach
            df_widths.add, widths.toarray()
            df_heights.add, heights.toarray()
        endforeach
                
        
        ; Check the coherence of DFs.
        df_xcors = list()
        foreach event, cluster do begin
            xcors = list()
            foreach df, event.df_list do xcors.add, df.xcor_max
            df_xcors.add, (xcors.toarray())[1:*]
        endforeach

        mean_xcors = fltarr(cluster.length)
        foreach xcors, df_xcors, ii do begin
            mean_xcors[ii] = mean(xcors)
            index = where(xcors gt 1, count)
            if count ne 0 then mean_xcors[ii] = 0   ; I see events contain xcor=2, which is for nan.
        endforeach
        
        ; Check r_square.
        df_r2 = list()
        timing_key = 'df_time'
        foreach event, cluster do begin
            timing_info = event.timing_info[timing_key]
            df_r2.add, timing_info.linear_fit_info.r_square
        endforeach
        df_r2 = df_r2.toarray()
        

        min_r2 = 0.5
        min_xcor = 0.6
        index = where(df_r2 ge min_r2 and mean_xcors ge min_xcor, count)
        if count gt 0 then begin
            tmp = max(df_r2[index], index)
            index = where(df_r2 eq tmp)
        endif else begin
            tmp = max(mean_xcors, index)
        endelse
        
        ;stop
        event_list.add, cluster[index[0]]
        
        event_list = list(event)
        dirname = ['test_analyze_overlap_cluster','group_'+string(cluster_id+1,format='(I0)')]
        dirname = ['test_analyze_overlap_cluster','test_event']
        azim_df_search_df_group_gen_diagnostic_plot, event_list, dirname=dirname, project=project
        
        ;stop
        continue
        

        ; AE and Dst.
        omni_read_index, data_time_range
        
        ; Plot each group.
        plot_vars = ['ae',df_probes+'_theta']
        nvar = n_elements(plot_vars)
        foreach probe, df_probes do begin
            var = probe+'_theta'
            data = get_var_data(var, in=plot_time_range)
            yrange = [-1,1]*max(abs(sg_autolim(data)))
            options, var, 'yrange', yrange
            options, var, 'labels', strupcase(probe)
        endforeach
        tplot_options, 'ystyle', 1
        
        foreach event, cluster do begin
            
        ;---Check width and height.
            df_widths = list()
            df_heights = list()
            foreach df, event.df_list do begin
                df_widths.add, df.width
                df_heights.add, df.height
            endforeach
            df_widths = df_widths.toarray()
            df_heights = df_heights.toarray()
            print, minmax(df_widths)
            print, minmax(df_heights)
            
            
        ;---Check arrow scatter.
            timing_key = 'df_time'
            timing_info = event.timing_info[timing_key]

            ; Get the good triads.
            triad_timing = timing_info.triad_timing
            good_triads = azim_df_filter_check_good_triads(triad_timing)
            v2ds = list()
            foreach probe_combo, good_triads do begin
                combo_info = triad_timing[probe_combo]
                v2ds.add, combo_info.v_hat*combo_info.v_mag
            endforeach
            foreach v2d, v2ds do print, v2d
            vmag = list()
            foreach v2d, v2ds do vmag.add, snorm(v2d)
            vmag = max(vmag.toarray())
            plot, [-1,1]*vmag, [-1,1]*vmag, /iso, /nodata
            foreach v2d, v2ds do begin
                plots, [0,v2d[0]], [0,v2d[1]], psym=-1
            endforeach
            stop
            
            
        ;---Plot DF.
            sgopen, 0, xsize=8, ysize=8, /inch
            poss = sgcalcpos(nvar, tmargin=2)
            tplot, plot_vars, position=poss, trange=plot_time_range
            foreach df, event.df_list do begin
                tvar = df.probe+'_theta'
                pos_id = where(tvar eq plot_vars)
                tpos = poss[*,pos_id]
                xrange = plot_time_range
                yrange = get_setting(tvar, 'yrange')
                plot, xrange, yrange, /noerase, /nodata, $
                    xstyle=5, xrange=xrange, $
                    ystyle=5, yrange=yrange, $
                    position=tpos
                xs = df.arrival_time+[-1,1]*0.5*df.width
                ys = df.value_range
                plots, xs[[0,1,1,0,0]], ys[[0,0,1,1,0]], /data
            endforeach
            stop
        endforeach
        stop
    endforeach
    
    
    

end