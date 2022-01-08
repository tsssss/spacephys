;+
; Update event info.
;-

function azim_df_load_event_info, project=project, reset=reset

    if n_elements(project) eq 0 then project = azim_df_load_project()
reset = 1
    the_key = 'events'
    if keyword_set(reset) then if project.haskey(the_key) then project.remove, the_key
    if project.haskey(the_key) then return, project[the_key]

    the_key = 'event_info_file'
    if ~project.haskey(the_key) then project[the_key] = join_path([project.data_dir,'event_info.txt'])
    event_info_file = project[the_key]

    nheader = 2
    lines = read_all_lines(event_info_file,skip_header=nheader)
    events = hash()
    foreach line, lines do begin
        tinfo = strsplit(line,' ',/extract)
        ntinfo = n_elements(tinfo)
        
        event_id = tinfo[0]
        time_range = time_double(tinfo[1:2])
        cc_time_range = time_double(tinfo[3:4])
        cc_time_lag = float(tinfo[5])
        ref_time = time_double(tinfo[6])
        
        if ntinfo ge 7 then begin
            mission_probes = strsplit(tinfo[7],',',/extract)
        endif else mission_probes = !null
        
        events[event_id] = dictionary($
            'id', event_id, $
            'time_range', time_range, $
            'cc_time_range', cc_time_range, $
            'cc_time_lag', cc_time_lag, $
            'ref_time', ref_time, $
            'mission_probes', mission_probes)
    endforeach
    project.events = events
    update_project, project
    
    return, events

end