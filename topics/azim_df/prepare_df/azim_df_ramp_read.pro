;+
; Read DFs.
;-
function azim_df_ramp_read, tline
    infos = strsplit(tline, ', ', /extract)
    ; The DF line: probe | time_range | value_range | obs_time | obs_r_sm | width | height
    probe = infos[0]
    time_range = time_double(infos[1:2])
    value_range = float(infos[3:4])
    obs_time = time_double(infos[5])
    obs_r_sm = float(infos[6:8])
    width = float(infos[9])
    height = float(infos[10])
    the_ramp = dictionary($
        'probe', probe, $
        'time_range', time_range, $
        'value_range', value_range, $
        'obs_time', obs_time, $
        'obs_r_sm', obs_r_sm, $
        'width', width, $
        'height', height )
    return, the_ramp

end
