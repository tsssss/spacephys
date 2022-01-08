;+
; Search DFs.
;-

function azim_df_search_df, search_setting, project=project, $
    reset=reset, test_time=test_time

;test_time = time_double('2007-11-16/20:00')
    retval = list()
    if n_elements(project) eq 0 then project = azim_df_load_project()
    tab = constant('4space')

;---Settings for selecting DF.
    the_key = 'search_df'
    if ~search_setting.haskey(the_key) then message, 'No settings for DF ...'
    search_info = search_setting[the_key]


;---Check if file exists, to avoid search again.
    file_suffix = search_info.file_suffix
    out_file = join_path([project.data_dir,file_suffix])
    log_file = out_file+'.log'
    if keyword_set(reset) then begin
        lprmsg, 'Resetting DF search ...'
        file_delete, out_file, /allow_nonexistent
        file_delete, log_file, /allow_nonexistent
        lprmsg, 'Clear memory ...'
        del_data, '*'
    endif
    if file_test(out_file) eq 1 then begin
        if file_lines(out_file) eq 0 then file_delete, out_file
    endif
    if file_test(log_file) eq 1 then begin
        if file_lines(log_file) eq 0 then file_delete, log_file
    endif

    if file_test(out_file) eq 0 or keyword_set(test_time) then begin
        event_id = 0

        if file_test(out_file) eq 0 then ftouch, out_file
        if file_test(log_file) eq 0 then ftouch, log_file

        lprmsg, 'Search dipolarizations ...', log_file
        lprmsg, 'Writing results to file: '+out_file+' ...', log_file
        lprmsg, '', log_file


        candidates = azim_df_search_triad(search_setting, project=project)
        foreach candidate, candidates do begin
            time_range = candidate.time_range
            if keyword_set(test_time) then if product(time_range-test_time) ge 0 then continue else stop
            lprmsg, '', log_file
            lprmsg, 'Processing candidate '+string(candidate.id,format='(I0)')+' ...', log_file
            lprmsg, strjoin(time_string(time_range),' to ')+' ...', log_file

            ; Collect DFs.
            probe_list = candidate.probe_list
            time_range_list = candidate.time_range_list
            df_list = list()
            foreach sector_time_range, time_range_list, sector_id do begin
                probes = probe_list[sector_id]
                nprobe = n_elements(probes)
                foreach probe, probes do begin
                    dfs = azim_df_read_df(sector_time_range, probe=probe)
                    if dfs.length eq 0 then continue
                    df_list.add, dfs, /extract
                endforeach
            endforeach
            ndf = df_list.length
            lprmsg, 'Found '+string(ndf,format='(I0)')+' DFs ...', log_file
            foreach df, df_list do begin
                scaled_height = azim_df_scale_theta(df.height, df.obs_mlt, width=scale_width)
                df['scaled_height'] = scaled_height
            endforeach

            ; Check ROI.
            mlts = fltarr(ndf)
            rxys = fltarr(ndf)
            foreach df, df_list, df_id do begin
                mlts[df_id] = df.obs_mlt
                rxys[df_id] = df.obs_rxy
            endforeach
            mlt_range = search_setting.search_roi.mlt_range
            rxy_range = search_setting.search_roi.rxy_range
            index = where(mlts gt mlt_range[0] and mlts lt mlt_range[1] and $
                rxys gt rxy_range[0] and rxys lt rxy_range[1], ndf)
            lprmsg, 'Found '+string(ndf,format='(I0)')+' DFs in ROI ...', log_file
            if ndf eq 0 then begin
                lprmsg, 'No DF in ROI, skip ...', log_file
                continue
            endif
            df_list = df_list[index]
            
            
            probes = strarr(ndf)
            foreach df, df_list, ii do probes[ii] = df.probe
            probes = sort_uniq(probes)


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
