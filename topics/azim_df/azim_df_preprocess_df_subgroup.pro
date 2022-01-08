;+
; Load all subgroups, sort by arrival time.
;
; Check triad.
; Check ROI.
;-


function azim_df_preprocess_df_subgroup, project=project, reset=reset

    if n_elements(project) eq 0 then project = azim_df_load_project()

;---Settings for DF group filter.
    the_key = 'preprocess_df_subgroup'
    if ~project.haskey(the_key) then project[the_key] = dictionary()
    search_df_group_filter = project[the_key]
    min_triad_count = 2
    hr_per_sec_2_deg_per_min = 15.*60
    default_settings = dictionary($
        'min_triad_count', min_triad_count, $
        'file_suffix', project.name+'_sorted_subgroup.txt')
    foreach key, default_settings.keys() do if ~search_df_group_filter.haskey(key) then search_df_group_filter[key] = default_settings[key]


;---Check if file exists, to avoid search again.
    file_suffix = search_df_group_filter.file_suffix
    file = join_path([project.data_dir,file_suffix])
    if keyword_set(reset) then begin
        lprmsg, 'Resetting to preprocess subgroups...'
        file_delete, file, /allow_nonexistent
        lprmsg, 'Clear memory ...'
        del_data, '*'
    endif

    ; Remove empty file.
    if file_test(file) then begin
        if file_lines(file) eq 0 then file_delete, file
    endif

    ; Read and return results.
    if ~keyword_set(test_time) then begin
        if file_test(file) then return, azim_df_subgroup_read_file(file)
    endif


;---Loop through DF groups.
    if file_test(file) eq 0 then ftouch, file
    events = azim_df_search_df_subgroup(project=project)
    tab = constant('4space')
    ndim = 3
    rxy_range = project.search_candidate.search_roi.rxy_range
    regions = project.search_candidate.search_roi.regions
    nvertex = 3
    pdyn = project.search_candidate.search_roi.pdyn
    triad_angle_range = project.search_candidate.search_triad.triad_angle_range
    max_probe_length = 5

    good_index = list()
    foreach df_group, events, event_id do begin
        df_list = df_group.df_list
        min_probe_count = df_list.length

    ;---Check duration.
        df_times = dblarr(df_list.length)
        foreach df, df_list, ii do df_times[ii] = df.arrival_time
        duration = (max(df_times)-min(df_times))/60

    ;---Check in ROI.
        rxys = fltarr(df_list.length)
        foreach df, df_list, ii do rxys[ii] = snorm(df.arrival_r_sm[0:1])
        index = lazy_where(rxys, '[]', rxy_range, count=count)
        if count lt min_probe_count then begin
            lprmsg, tab+tab+'Not enough probes in ROI (Rxy), skip ...'
            continue
        endif
        df_list = df_list[index]

        ; MLT.
        mlts = fltarr(df_list.length)
        foreach df, df_list, ii do mlts[ii] = azim_df_calc_pseudo_mlt(df.arrival_r_sm)
        region_name = df_group.region
        mlt_range = regions[region_name].mlt_range
        index = lazy_where(mlts, '[]', mlt_range, count=count)
        if count lt min_probe_count then begin
            lprmsg, tab+tab+'No enough probes in ROI (MLT), skip ...'
            continue
        endif
        df_list = df_list[index]

        ; Magnetopause.
        r_sms = fltarr(df_list.length,ndim)
        foreach df, df_list, ii do r_sms[ii,*] = df.arrival_r_sm
        index = where(check_if_in_magn(r_sms, dynamic_pressure=pdyn) eq 1, count)
        if count lt min_probe_count then begin
            lprmsg, tab+tab+'No enough probes in ROI (M/pause), skip ...'
            continue
        endif
        df_list = df_list[index]


        ; Triad.
        r_sms = fltarr(df_list.length,ndim)
        times = fltarr(df_list.length)
        probes = strarr(df_list.length)
        foreach df, df_list, ii do begin
            r_sms[ii,*] = df.arrival_r_sm
            times[ii] = df.arrival_time
            probes[ii] = df.probe
        endforeach
        combos = choose_from(probes, nvertex)
        ncombo = n_elements(combos)
        triad_flags = intarr(ncombo)+1  ; 1 for good triad.
        foreach combo, combos, jj do begin
            probe_index = intarr(nvertex)
            foreach probe, combo, kk do probe_index[kk] = where(probes eq probe)
            the_r_sms = transpose(r_sms[probe_index,*])
            the_r_sms[2,*] = 0
            the_r_sms = reform(the_r_sms, [1,3,3])
            triad_angles = triangle_angles(the_r_sms)
            index = lazy_where(triad_angles, '[]', triad_angle_range, count=count)
            if count ne nvertex then triad_flags[jj] = 0
            
            
        endforeach
        index = where(triad_flags eq 1, count)
        if count lt min_triad_count then begin
            lprmsg, tab+tab+'No enough good triad, skip ...'
            continue
        endif
        combos = combos[index]
        combo_probes = list()
        foreach combo, combos do combo_probes.add, combo, /extract
        combo_probes = sort_uniq(combo_probes.toarray())
        nprobe = n_elements(df_group.probes)
        if n_elements(combo_probes) lt nprobe then begin
            lprmsg, tab+tab+'Not all probes form good triads, skip ...'
            continue
        endif

        good_index.add, event_id
    endforeach
    good_index = good_index.toarray()
    events = events[good_index]
    
;---Filter by duration.
    events = azim_df_filter_duration(events, project=project)
    nevent = n_elements(events)

;---Sort by time.
    event_start_times = dblarr(nevent)
    foreach df_group, events, ii do event_start_times[ii] = df_group.time_range[0]
    index = sort(event_start_times)
    events = events[index]
    foreach df_group, events, ii do begin
        df_group['event_id'] = ii+1
        azim_df_subgroup_write_file, df_group, filename=file
    endforeach
    
    
    return, azim_df_subgroup_read_file(file)
end

reset = 0
if n_elements(project) eq 0 then project = azim_df_load_project()
events = azim_df_preprocess_df_subgroup(project=project, reset=reset)
end
