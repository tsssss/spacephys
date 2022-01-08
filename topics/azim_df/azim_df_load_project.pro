;+
; Return a dictionary of settings of the project.
;-

function azim_df_load_project, reset=reset


    project = load_project('azim_df')
    project_updated = 0

    if keyword_set(reset) then begin
        del_data, project.var
        file_delete, project.file, /allow_nonexistent
        project = load_project('azim_df')
    endif

;---Overall settings.
    the_key = 'time_step'
    if ~project.haskey(the_key) then begin
        project[the_key] = 10. ; sec.
        project_updated = 1
    endif

    the_key = 'scale_width'
    if ~project.haskey(the_key) then begin
        project[the_key] = 4.   ; MLT.
        project_updated = 1
    endif

    the_key = 'triad_min_count'
    if ~project.haskey(the_key) then begin
        project[the_key] = 2.   ; #.
        project_updated = 1
    endif

    the_key = 'triad_angle_range'
    if ~project.haskey(the_key) then begin
        project[the_key] = [15.,165]    ; deg.
        project_updated = 1
    endif

    the_key = 'probe_min_count'
    if ~project.haskey(the_key) then begin
        project[the_key] = 4.   ; #.
        project_updated = 1
    endif

    the_key = 'roi_min_duration'
    if ~project.haskey(the_key) then begin
        project[the_key] = 3600.    ; sec.
        project_updated = 1
    endif

    the_key = 'overall_roi'
    if ~project.haskey(the_key) then begin
        project[the_key] = dictionary($
            'mlt_range', [-1,1]*9, $
            'rxy_range', [4.,30], $
            'pdyn', 10.)
        project_updated = 1
    endif


;---Regions in MLT.
    the_key = 'regions'
    if ~project.haskey(the_key) then begin
        regions = dictionary($
            'pre_midn', dictionary('mlt_range', [-9,0]), $
            'post_midn', dictionary('mlt_range', [0,9]), $
            'around_midn', dictionary('mlt_range', [-1,1]*3) )
        project[the_key] = regions
        project_updated = 1
    endif


;---Search names.
    the_key = 'search_types'
    if ~project.haskey(the_key) then begin
        search_types = orderedhash()

    ;---Beyond 15 Re. THB/C must be beyond_15Re, >4 probes in ROI.
        search_types['beyond_15Re'] = dictionary($
            'name', 'beyond_15Re', $
            'time_range', time_double(['2007-03-01','2009-09-01']), $
            'probes', 'th'+letters('e'), $
            'roi_min_count', 4, $
            'data_file_suffix', 'azim_df_primitive_data_themis.cdf')

    ;---Within 15 Re. >5 probes in ROI.
        search_types['within_15Re'] = dictionary($
            'name', 'within_15Re', $
            'time_range', time_double(['2012-10-01','2017-10-01']), $
            'probes', ['rbsp'+letters('b'),'th'+['a','d','e'],'mms1','g'+['13','14','15']], $   ; turns out no thb/thc candidate.
            'roi_min_count', 5, $
            'data_file_suffix', 'azim_df_primitive_data.cdf')

        project[the_key] = search_types
    endif


    the_key = 'all_probes'
    if ~project.haskey(the_key) then begin
        all_probes = []
        foreach search_type, project.search_types do all_probes = [all_probes, search_type.probes]
        all_probes = reverse(sort_uniq(all_probes)) ; reverse is a better order: themis are the first.
        project[the_key] = all_probes
        project_updated = 1
    endif

    the_key = 'probe_infos'
    if ~project.haskey(the_key) then begin
        all_probes = project.all_probes
        probe_colors = smkarthm(50,250,n_elements(all_probes), 'n')
        probe_color_ct = 40
        foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)
        probe_infos = dictionary()
        foreach probe, all_probes, ii do begin
            probe_info = resolve_probe(probe)
            probe_info['color'] = probe_colors[ii]
            probe_infos[probe] = probe_info
        endforeach
        project[the_key] = probe_infos
        project_updated = 1
    endif

    if project_updated then update_project, project
    return, project

end
