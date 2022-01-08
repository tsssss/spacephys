;+
; Read roi_list from file.
; Adopted from azim_df_search_candidate_read_file
;-

function azim_dp_roi_read, filename=file

    lines = read_all_lines(file)
    nline = n_elements(lines)

    roi_list = list()
    for line_id=0,nline-1 do begin
        tline = lines[line_id]
        if keyword_set(test) then lprmsg, tline

        infos = strsplit(tline,' ',/extract)
        ninfo = n_elements(infos)
        if ninfo eq 5 then begin
            time_range = time_double(infos[0:1])
;            duration = float(infos[2])
            nsection = fix(infos[3])
            all_probes = strsplit(infos[4],',',/extract)
            roi = dictionary($
                'time_range', time_range, $
                'nsection', nsection, $
                'all_probes', all_probes, $
                'time_range_list', list(), $
                'probe_list', list() )
            for subline_id=line_id+1, line_id+nsection do begin
                tline = lines[subline_id]
                if keyword_set(test) then lprmsg, tline
                infos = strsplit(tline,' ',/extract)
                ninfo = n_elements(infos)
                roi.time_range_list.add, time_double(infos[1:2])
                roi.probe_list.add, strsplit(infos[3],',',/extract)
            endfor

            ; Reached next roi.
            roi_list.add, roi
            line_id += nsection
        endif else continue
    endfor

    return, roi_list
end
