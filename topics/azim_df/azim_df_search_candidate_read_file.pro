
function azim_df_search_candidate_read_file, file
    lines = read_all_lines(file)
    nline = n_elements(lines)

    candidates = list()
    for ii=0, nline-1 do begin
        tline = lines[ii]
        if keyword_set(test) then lprmsg, tline
        infos = strsplit(tline,' ',/extract)
        ninfo = n_elements(infos)
        if ninfo eq 8 then begin
            ; The major line: candidate_id | start_time | end_time | region | search | duration | nsection | probes
            candidate_id = fix(infos[0])
            time_range = time_double(infos[1:2])
            region_name = infos[3]
            search_name = infos[4]
            duration = float(infos[5])
            nsection = fix(infos[6])
            all_probes = strsplit(infos[7],',',/extract)
            candidate = dictionary($
                'id', candidate_id, $
                'time_range', time_range, $
                'duration', duration, $  ; in min.
                'nsection', nsection, $
                'search_name', search_name, $
                'region', region_name, $
                'all_probes', all_probes, $
                'time_range_list', list(), $
                'probe_list', list())

            for jj=ii+1, ii+nsection do begin
                tline = lines[jj]
                if keyword_set(test) then lprmsg, tline
                infos = strsplit(tline,' ',/extract)
                ninfo = n_elements(infos)
                ; The section line: section id | start_time | end_time | probes
                candidate.time_range_list.add, time_double(infos[1:2])
                candidate.probe_list.add, strsplit(infos[3],',',/extract)
            endfor

            ; Reached next candidate.
            candidates.add, candidate
            ii = ii+nsection
        endif else continue
    endfor

    return, candidates
end
