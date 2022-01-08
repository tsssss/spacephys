pro azim_df_ramp_write, ramp, filename=file, prefix=prefix
    tab = constant('4space')
    max_probe_length = 5
    tformat='YYYY-MM-DD/hh:mm:ss'
    if n_elements(prefix) eq 0 then prefix = ''
    msg = prefix[0]+extend_string(ramp.probe, length=max_probe_length)+tab+$
        strjoin(time_string(ramp.time_range,tformat=tformat),' ')+tab+$
        strjoin(string(ramp.value_range,format='(F7.2)'),' ')+tab+$
        time_string(ramp.obs_time,tformat=tformat)+tab+$
        strjoin(string(ramp.obs_r_sm,format='(F6.2)'),',')+tab+$
        string(ramp.width,format='(I5)')+tab+$
        string(ramp.height,format='(F7.2)')
    lprmsg, msg, file

end
