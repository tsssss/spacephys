
function azim_df_triad_read, tline
    infos = strsplit(tline,', ', /extract)
    the_triad = dictionary($
        'probes', infos[0:2], $
        'times', time_double(infos[3:5]), $
        'angles', float(infos[6:8]), $
        'center_r_sm', float(infos[9:11]), $
        'omega_obs_time', float(infos[12]), $           ; in deg/s, +:westward, -:eastward.
        'omega_xcor', float(infos[13]), $
        'vmag_obs_time', float(infos[14]), $            ; in km/s, abs value.
        'vmag_xcor', float(infos[15]), $
        'vhat_obs_time', float(infos[16:17]), $
        'vhat_xcor', float(infos[18:19]), $
        'r_sms', reform(float(infos[20:28]),[3,3]) )    ; [0,*] is [x,y,z] of probe_1, etc.
    return, the_triad
end
