;+
; Update event info.
;-

pro azim_prop_update_event_info, project

    event_info_file = join_path([project.data_dir,'event_info.txt'])
    project.event_info_file = event_info_file
    nheader = 2
    lines = read_all_lines(event_info_file,skip_header=nheader)
    events = hash()
    foreach line, lines do begin
        tinfo = strsplit(line,' ',/extract)
        event_id = tinfo[0]
        time_range = time_double(tinfo[1:2])
        probes = strsplit(tinfo[3],',',/extract)
        cc_time_range = time_double(tinfo[4:5])
        cc_time_lag = float(tinfo[6])
        ref_time = time_double(tinfo[7])
        events[event_id] = dictionary($
            'id', event_id, $
            'time_range', time_range, $
            'probes', probes, $
            'cc_time_range', cc_time_range, $
            'cc_time_lag', cc_time_lag, $
            'ref_time', ref_time)
    endforeach
    project.events = events
    azim_prop_update_project, project
end
