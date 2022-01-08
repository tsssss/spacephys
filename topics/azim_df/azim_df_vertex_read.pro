
function azim_df_vertex_read, tline
    infos = strsplit(tline, ', ', /extract)
    the_df = dictionary($
        'probe', infos[0], $
        'obs_time', time_double(infos[1]), $
        'width', float(infos[2]), $
        'height', float(infos[3]), $
        'scaled_height', float(infos[4]), $
        'obs_mlt', float(infos[5]), $
        'obs_rxy', float(infos[6]), $
        'time_range', time_double(infos[7:8]), $
        'theta_range', float(infos[9:10]), $
        'obs_r_sm', float(infos[11:13]))
    return, the_df
end