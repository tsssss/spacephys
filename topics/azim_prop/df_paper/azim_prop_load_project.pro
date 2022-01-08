;+
; Load project info, if no file exists, then initialize it.
;-

pro azim_prop_init_project, project


;---What are the probes to search over.
    probes = ['rbspa','rbspb','g13','g14','g15','tha','thd','the','mms1']
    probe_colors = sgcolor(['red','tomato','orange','gold','lime_green','blue','deep_sky_blue','purple','black'])
    probe_symbols = replicate(6,n_elements(probes))

    project.probes = probes
    foreach probe, project.probes, ii do begin
        mission = resolve_probe(probe)
        project[probe] = dictionary()
        project[probe].probe = mission.probe
        project[probe].prefix = mission.name+mission.probe
        project[probe].short_name = mission.short_name+mission.probe
        project[probe].routine_name = mission.routine_name
        project[probe].color = probe_colors[ii]
        project[probe].psym = probe_symbols[ii]
    endforeach


;---Constants.
    constant = dictionary()
    constant.rgb = sgcolor(['red','green','blue'])
    constant.xyz = ['x','y','z']
    constant.deg = 180d/!dpi
    constant.rad = !dpi/180d
    constant.re = 6378d
    constant.secofday = 86400d
    project.constant = constant


;---Event info.
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


function azim_prop_load_project

    root_dir = join_path([googledir(),'works','works','azim_prop','ten_event_paper'])

;---Basic info.
    project = dictionary()
    project.name = 'azim_prop'
    project.var = project.name+'_project_info'
   ;project.root_dir = join_path([shomedir(),'azim_prop'])
    project.root_dir = root_dir
    project.data_dir = join_path([project.root_dir,'data'])
    project.plot_dir = join_path([project.root_dir,'plot'])
    project.file = join_path([project.data_dir,project.name+'_project_info.tplot'])

    if file_test(project.file) eq 0 then azim_prop_init_project, project
    if file_test(project.file) eq 0 then message, 'Something wrong in loading the project info ...'

    tplot_restore, filename=project.file
    project = get_var_data(project.var)

    ; This fix a problem when switching between windows and mac.
    if file_test(project.root_dir,/directory) eq 0 then begin
        project.root_dir = root_dir
        project.data_dir = join_path([project.root_dir,'data'])
        project.plot_dir = join_path([project.root_dir,'plot'])
        project.file = join_path([project.data_dir,project.name+'_project_info.tplot'])
        azim_prop_update_project, project
    endif

    return, project

end

project = azim_prop_load_project()
end
