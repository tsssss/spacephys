;+
; Return a dictionary of settings of the project.
;-

function pflux_grant_load_project, reset=reset

    project_name = 'pflux_grant'
    project = load_project(project_name)
    project_updated = 0

    if keyword_set(reset) then begin
        del_data, project.var
        file_delete, project.file, /allow_nonexistent
        project = load_project(project_name)
    endif

;---Overall settings.
    the_key = 'all_probes'
    if ~project.haskey(the_key) then begin
        all_probes = 'rbsp'+['a','b']
        project[the_key] = all_probes
        project_updated = 1
    endif

    the_key = 'probe_infos'
    if ~project.haskey(the_key) then begin
        all_probes = project.all_probes
        probe_infos = dictionary()
        foreach probe, all_probes, ii do begin
            probe_info = resolve_probe(probe)
            probe_infos[probe] = probe_info
        endforeach
        project[the_key] = probe_infos
        project_updated = 1
    endif

    the_key = 'common_time_step'
    if ~project.haskey(the_key) then begin
        project[the_key] = 1d/16
        project_updated = 1
    endif

    the_key = 'pflux_calc_setting'
    if ~project.haskey(the_key) then begin
        pflux_calc_setting = dictionary()
        all_probes = project.all_probes
        start_time = '2012-10-01'
        end_time = '2015-10-01'
        foreach probe, all_probes do begin
            ;end_time = (probe eq 'rbspa')? '2015-10-01': '2018-04-01'
            time_range = time_double([start_time,end_time])
            pflux_calc_setting[probe] = dictionary('time_range', time_range)
        endforeach
        
        period_range = [1d,1800]
        p2s = (wavelet_info('morlet'))[7]   ; t2s.
        scale_range = period_range*p2s
        pflux_calc_setting['period_range'] = period_range
        pflux_calc_setting['scale_range'] = scale_range
        dj = 1d/8
        scales = smkgmtrc(scale_range[0], scale_range[1], 2^dj, 'dx')
        periods = scales/p2s
        freqs = 1/periods
        pflux_calc_setting['scales'] = scales
        pflux_calc_setting['periods'] = periods
        pflux_calc_setting['freqs'] = freqs
        pflux_calc_setting['scale_info'] = $
            {s0:scale_range[0], s1:scale_range[1], dj:dj, ns:0}

        project[the_key] = pflux_calc_setting
        project_updated = 1
    endif

    if project_updated then update_project, project
    return, project

end