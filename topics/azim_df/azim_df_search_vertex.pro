;+
; Search vertex for ROI candidates.
;-
function azim_df_search_vertex, search_setting, project=project, $
    reset=reset, test_time=test_time

;test_time = time_double('2007-11-16/20:00')
    retval = list()
    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()

;---Settings for search DF in ROI.
    the_key = 'search_vertex'
    if ~search_setting.haskey(the_key) then message, 'No settings for '+the_key+'...'
    search_info = search_setting[the_key]
    mlt_range = search_setting.search_roi.mlt_range
    rxy_range = project.overall_roi.rxy_range
    probe_min_count = project.probe_min_count


;---Check if file exists, to avoid search again.
    file_suffix = search_info.file_suffix
    out_file = join_path([project.data_dir,file_suffix])
    if keyword_set(reset) then begin
        lprmsg, 'Resetting '+the_key+'...'
        file_delete, out_file, /allow_nonexistent
    endif
    if file_test(out_file) eq 1 then begin
        if file_lines(out_file) eq 0 then file_delete, out_file
    endif
    log_file = -1    ; output to console.


;---Do search.
    if file_test(out_file) eq 0 or keyword_set(test_time) then begin
        event_id = 0
        if file_test(out_file) eq 0 then ftouch, out_file

    ;---Load candidates.
        candidates = azim_df_search_roi(search_setting, project=project)
        df_count = 0
        foreach candidate, candidates do begin
            ; Collect DFs.
            probe_list = candidate.probe_list
            time_range_list = candidate.time_range_list
            df_list = list()
            foreach sector_time_range, time_range_list, sector_id do begin
                probes = probe_list[sector_id]
                nprobe = n_elements(probes)
                foreach probe, probes do begin
                    dfs = azim_df_read_df(sector_time_range, probe=probe, project=project)
                    if dfs.length eq 0 then continue
                    df_list.add, dfs, /extract
                endforeach
            endforeach
            ndf = df_list.length
            lprmsg, 'Found '+string(ndf,format='(I0)')+' DFs ...', log_file
            if ndf lt probe_min_count then begin
                lprmsg, 'No enough DF, skip ...', log_file
                continue
            endif

            ; Check ROI.
            obs_mlts = fltarr(ndf)
            obs_rxys = fltarr(ndf)
            foreach df, df_list, df_id do begin
                obs_mlts[df_id] = df.obs_mlt
                obs_rxys[df_id] = df.obs_rxy
            endforeach
            index = where(obs_mlts gt mlt_range[0] and obs_mlts lt mlt_range[1] and $
                obs_rxys gt rxy_range[0] and obs_rxys lt rxy_range[1], ndf)
            lprmsg, 'Found '+string(ndf,format='(I0)')+' DFs in ROI ...', log_file
            if ndf lt probe_min_count then begin
                lprmsg, 'No enough DF in ROI, skip ...', log_file
                continue
            endif
            df_list = df_list[index]

            ; Check probe.
            ndf = df_list.length
            probes = strarr(ndf)
            foreach df, df_list, ii do probes[ii] = df.probe
            probes = sort_uniq(probes)
            nprobe = n_elements(probes)
            if nprobe lt probe_min_count then begin
                lprmsg, 'No enough probe in ROI, skip ...', log_file
                continue
            endif

            ; Sort by obs_time.
            obs_times = dblarr(ndf)
            foreach df, df_list, ii do obs_times[ii] = df.obs_time
            index = sort(obs_times)
            df_list = df_list[index]

            ; Write to file.
            event_id += 1
            df_group = dictionary($
                'id', event_id, $
                'time_range', candidate.time_range, $
                'search_name', candidate.search_name, $
                'region', candidate.region, $
                'probes', probes, $
                'df_list', df_list, $
                'triad_list', list(), $
                'edge_list', list())

            azim_df_subgroup_write_file, df_group, filename=out_file
        endforeach
    endif

    lprmsg, ''
    lprmsg, 'Read DFs from file ...'
    events = azim_df_subgroup_read_file(out_file)
    return, events

end
