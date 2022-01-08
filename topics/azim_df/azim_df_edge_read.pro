
function azim_df_edge_read, tline
    infos = strsplit(tline,', ', /extract)
    the_edge = dictionary($
        'probes', infos[0:1], $
        'time_lag_df', float(infos[2]), $   ; probe_2-probe_1, must be +.
        'time_lag_xcor', float(infos[3]), $ ; probe_2-probe_1, can be +/-.
        'xcor_max', float(infos[4]), $
        'xcor_err', float(infos[5]), $
        'time_lag_diff', float(infos[6]), $
        'time_lag_ratio', float(infos[7]) )
    return, the_edge
end
