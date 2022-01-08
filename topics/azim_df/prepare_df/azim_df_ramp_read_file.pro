;+
; Read DFs.
;-
function azim_df_ramp_read_file, file

    lines = read_all_lines(file)
    nline = n_elements(lines)

    ramps = list()
    for line_id=0,nline-1 do begin
        tline = lines[line_id]
        if keyword_set(test) then lprmsg, tline
        ramps.add, azim_df_ramp_read(tline)
    endfor

    return, ramps

end

file = '/Volumes/GoogleDrive/My Drive/works/works/azim_df/data/azim_df_detect_ramp_beyond_15Re_tha.txt'
ramps = azim_df_ramp_read_file(file)
end
