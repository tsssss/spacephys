
    if n_elements(project) eq 0 then project = azim_df_load_project()
    search_settings = project.search_settings
    roi_candidates = azim_df_search_candidate(project=project, /get_roi_candidates)
    time_step = project.orbit_time_step
    tab = constant('4space')
    file = join_path([project.data_dir,'azim_df_analysis_roi.txt'])
    if file_test(file) then file_delete, file
    ftouch, file

    foreach search_setting, search_settings do begin
        full_time_range = search_setting.time_range
        search_name = search_setting.name
        lprmsg, '', file
        lprmsg, 'Search name: '+search_name+' ...', file

        combined_candidates = list()
        foreach candidate, roi_candidates do begin
            if candidate.search_type ne search_name then continue
            ;lprmsg, time_string(candidate.time_range)
            if n_elements(current_candidates) eq 0 then begin
                current_candidates = dictionary($
                    'search_type', search_name, $
                    'region', candidate.region, $
                    'time_ranges', list(candidate.time_range), $
                    'probes', list(candidate.probes))
            endif else begin
                last_time_range = current_candidates.time_ranges[-1]
                time_diff = candidate.time_range[0]-last_time_range[1]
                if abs(time_diff) le time_step then begin
                ;---Current candidate should be joined to previous ones.
                    current_candidates.time_ranges.add, candidate.time_range
                    current_candidates.probes.add, candidate.probes
                endif else begin
                ;---Current candidate starts a new section.
                    combined_candidates.add, current_candidates
                    current_candidates = dictionary($
                    'search_type', search_name, $
                    'region', candidate.region, $
                    'time_ranges', list(candidate.time_range), $
                    'probes', list(candidate.probes))
                endelse
            endelse
        endforeach


        foreach combined_candidate, combined_candidates, ii do begin
            time_ranges = combined_candidate.time_ranges
            probes = combined_candidate.probes
            lprmsg, '', file
            lprmsg, string(ii+1,format='(I5)')+' combined candidate '+combined_candidate.region+' '+string(time_ranges.length,format='(I2)')+' sections ...', file
            foreach time_range, time_ranges, jj do begin
                lprmsg, tab+string(jj+1,format='(I2)')+' section: '+strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmm'),' to ')+$
                    tab+string(total(time_range*[-1,1])/60,format='(I5)')+$
                    tab+strjoin(probes[jj],','), file
            endforeach
        endforeach
    endforeach
end
