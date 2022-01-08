;+
; Based on results from pflux_grant.
;-
function pflux_survey_load_project

    project_name = 'pflux_survey'
    project = load_project(project_name)
    project_updated = 0

    if keyword_set(reset) then begin
        del_data, project.var
        file_delete, project.file, /allow_nonexistent
        project = load_project(project_name)
    endif

    the_key = 'data_file_suffix'
    if ~project.haskey(the_key) then begin
        base_name = 'pflux_survey_pflux_data.cdf'
        project[the_key] = base_name
        project_updated = 1
    endif

    the_key = 'all_probes'
    if ~project.haskey(the_key) then begin
        all_probes = 'rbsp'+['a','b']
        project[the_key] = all_probes
        project_updated = 1
    endif

    the_key = 'time_range'
    if ~project.haskey(the_key) then begin
        time_range = time_double(['2012-10-01','2015-10-01'])
        project[the_key] = time_range
        project_updated = 1
    endif
    time_range = project[the_key]

    the_key = 'probe_infos'
    if ~project.haskey(the_key) then begin
        all_probes = project.all_probes
        probe_infos = dictionary()
        foreach probe, all_probes, ii do begin
            probe_info = resolve_probe(probe)
            probe_info['time_range'] = time_range
            probe_infos[probe] = probe_info
        endforeach
        project[the_key] = probe_infos
        project_updated = 1
    endif

    if project_updated then update_project, project
    return, project

end
