;+
; Read candidate_list from file.
; Adopted from azim_df_search_candidate_read_file
;-

function azim_dp_candidate_read, filename=file

    lines = read_all_lines(file)
    nline = n_elements(lines)

    candidate_list = list()
    for line_id=0,nline-1 do begin
        tline = lines[line_id]
        if keyword_set(test) then lprmsg, tline

        infos = strsplit(tline,' ',/extract)
        ninfo = n_elements(infos)
        if ninfo eq 6 then begin
            time_range = time_double(infos[0:1])
;            duration = float(infos[2])
            mlt_range = float(infos[3:4])
            probes = strsplit(infos[5],',',/extract)
            candidate = dictionary($
                'time_range', time_range, $
                'mlt_range', mlt_range, $
                'probes', probes )

            ; Reached next candidate.
            candidate_list.add, candidate
        endif else continue
    endfor

    return, candidate_list
end
