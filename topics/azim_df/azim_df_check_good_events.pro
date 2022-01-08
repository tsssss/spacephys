event_list = list()

; post_midn, beyond_15Re.
region_name = 'post_midn'
search_name = 'beyond_15Re'
    ; Good ones.
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2007-11-20/17:01','2007-11-20/17:16']), $
        'probes', ['tha','thc','the','thd'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2008-01-09/01:09','2008-01-09/11:39']), $
        'probes', ['thc','tha','the','thd'])
    ; Need to fix.
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2007-12-20/09:27','2007-12-20/09:38']), $
        'probes', ['thb','thc','thd','the','tha'])

; post_midn, within_15Re.
region_name = 'post_midn'
search_name = 'within_15Re'
    ; Good ones.
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2014-08-28/10:11','2014-08-28/10:47']), $
        'probes', ['thd','the','g15','tha','rbspb','g13'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2016-10-13/12:22','2016-10-13/12:45']), $
        'probes', ['rbspa','rbspb','g15','thd','g14','g13'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2016-11-26/16:23','2016-11-26/16:57']), $
        'probes', ['tha','thd','the','g15'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2016-12-11/09:46','2016-12-11/10:03']), $
        'probes', ['tha','thd','the','g13'])
    ; Need to fix.
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2013-05-01/14:11','2013-05-01/14:16']), $
        'probes', ['g15','tha','thd','the'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2014-07-07/14:22','2014-07-07/14:55']), $
        'probes', ['the','rbspa','rbspb','thd'])

; pre_midn, beyond_15Re.
region_name = 'pre_midn'
search_name = 'beyond_15Re'
    ; Good ones.
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2008-02-04/01:16','2008-02-04/01:28']), $
        'probes', ['tha','thd','the','thc'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2008-03-20/12:06','2008-03-20/12:17']), $
        'probes', ['thc','thb','thd','the'])
    event_list.add, dictionary($
            'region', region_name, $
            'search_name', search_name, $
            'time_range', time_double(['2009-04-07/07:05','2009-04-07/07:08']), $
            'probes', ['thd','tha','thb','thc'])
    ; Need fix.
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2009-03-30/01:48','2009-03-30/02:23']), $
        'probes', ['tha','thd','thc','the','thb'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2009-04-13/03:47','2009-04-13/04:06']), $
        'probes', ['tha','thc','thd','the'])


    ; pre_midn, within_15Re.
    region_name = 'pre_midn'
    search_name = 'within_15Re'
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2014-12-26/01:06','2014-12-26/01:19']), $
        'probes', ['g13','tha','the','g15'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2016-06-15/04:19','2016-06-15/04:39']), $
        'probes', ['g13','mms1','g14','g15','tha'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2017-03-28/03:00','2017-03-28/03:14']), $
        'probes', ['tha','the','thd','rbspa','g15'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2017-03-31/02:47','2017-03-31/03:01']), $
        'probes', ['g13','tha','the','thd','g15'])



    project = azim_df_load_project()
    ; Read all events.
    search_step = 'df_group'
    events = list()
    events.add, azim_df_search_post_midn_events(search_step=search_step), /extract
    events.add, azim_df_search_pre_midn_events(search_step=search_step), /extract
    good_events = list()
    foreach event, event_list do begin
        lprmsg, 'Finding event :'+strjoin(time_string(event.time_range),' to ')+' ...'
        foreach tmp, events do begin
            if event.region ne tmp.region then continue
            if event.search_name ne tmp.search_name then continue
            dtimes = event.time_range-tmp.time_range
            if total(dtimes) gt 120 then continue
            if n_elements(event.probes) eq n_elements(tmp.probes) then begin
                good_events.add, tmp
                lprmsg, 'Found the event ...'
                break
            endif
        endforeach
    endforeach



;    foreach event, good_events do begin
;        mlts = list()
;        foreach df, event.df_list do mlts.add, df.obs_mlt
;        mlts = mlts.toarray()
;        print, min(abs(mlts)), total(minmax(mlts)*[-1,1])
;    endforeach
;    stop


    scaled_height = list()
    foreach event, good_events do begin
        heights = list()
        foreach df, event.df_list do heights.add, df.height
        print, minmax(heights.toarray())
        scaled_height.add, heights, /extract
    endforeach
    print, minmax(scaled_height.toarray())
stop



min_scaled_height = 10
scale_width = 3
tab = constant('4space')
foreach event, event_list, event_id do begin
    df_list = list()
    time_range = event.time_range+[-1,1]*60.
    probes = event.probes
    foreach probe, probes do begin
        df = azim_df_read_df(time_range, probe=probe, scale_width=scale_width)
        ndf = n_elements(df)
        heights = fltarr(df.length)
        foreach tmp, df, ii do heights[ii] = tmp.scaled_height
        index = where(heights ge min_scaled_height, count)
        if count ne 1 then tmp = max(heights, index)
        df_list.add, df[index[0]]
    endforeach
    event['df_list'] = df_list
    event['id'] = event_id+1


    lprmsg, '', log_file
    lprmsg, 'Event '+strjoin(time_string(event.time_range),' to ')+' ...', log_file
    foreach df, event.df_list do azim_df_vertex_write, df, filename=log_file, prefix=tab


    settings = dictionary($
        'obs_time_diff_range', [-60,1800], $
        'obs_time_max_diff', 900., $
        'obs_time_max_diff_count', 2, $
        'xcor_section_ratio', [-1,1]*1.5, $
        'xcor_min_value', 0.5, $
        'file_suffix', !null )
    event = azim_df_subgroup_analyze(event, project=project, log_file=log_file, settings=settings)
    if event.haskey('df_list') then foreach df, event.df_list do azim_df_vertex_write, df, filename=log_file, prefix=tab
    if event.haskey('edge_list') then foreach edge, event.edge_list do azim_df_edge_write, edge, filename=log_file, prefix=tab
    if event.haskey('triad_list') then foreach triad, event.triad_list do azim_df_triad_write, triad, filename=log_file, prefix=tab
endforeach


stop


scaled_heights = list()
foreach event, event_list do begin
    foreach df, event.df_list do scaled_heights.add, df.scaled_height
endforeach
scaled_heights = scaled_heights.toarray()
print, minmax(scaled_heights)
stop

durations = list()
foreach event, event_list do begin
    durations.add, total(event.time_range*[-1,1])
endforeach
durations = durations.toarray()
print, 'Durations (min): '
print, durations/60.

nprobes = list()
foreach event, event_list do begin
    nprobes.add, event.df_list.length
endforeach
nprobes = nprobes.toarray()
print, 'Avg durations (min): '
print, durations/(nprobes-1)/60.
stop

dtime_list = list()
foreach event, event_list do begin
    df_times = list()
    foreach df, event.df_list do df_times.add, df.obs_time
    df_times = df_times.toarray()
    dtime_list.add, df_times[1:*]-df_times[0:-2]
endforeach
stop

widths = list()
foreach event, event_list do begin
    foreach df, event.df_list do widths.add, df.width
endforeach
widths = widths.toarray()
print, minmax(widths)

heights = list()
foreach event, event_list do begin
    foreach df, event.df_list do heights.add, df.height
endforeach
heights = heights.toarray()
print, minmax(heights)


project = azim_df_load_project()
dirname = 'check_good_events'
;azim_df_subgroup_gen_diagnostic_plot, event_list, dirname=dirname, project=project
log_file = join_path([project.plot_dir,'diagnostic_plot',dirname,dirname+'.log'])
if file_test(log_file) eq 1 then file_delete, log_file
ftouch, log_file
foreach event, event_list do begin
    lprmsg, '', log_file
    lprmsg, 'Event '+strjoin(time_string(event.time_range),' to ')+' ...', log_file
    foreach df, event.df_list do azim_df_vertex_write, df, filename=log_file, prefix='    '
endforeach
end
