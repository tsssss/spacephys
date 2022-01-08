
pro azim_df_vertex_write, df, filename=file, prefix=prefix
    tab = constant('4space')
    max_probe_length = 5
    tformat='YYYY-MM-DD/hh:mm:ss'
    if n_elements(prefix) eq 0 then prefix = ''
    msg = prefix[0]+extend_string(df.probe,length=max_probe_length)+tab+$
        time_string(df.obs_time,tformat=tformat)+tab+$
        string(df.width,format='(I5)')+tab+$
        string(df.height,format='(F7.2)')+tab+$
        string(df.scaled_height,format='(F10.1)')+tab+$
        string(df.obs_mlt,format='(F5.2)')+tab+$
        string(df.obs_rxy,format='(F5.2)')+tab+$
        strjoin(time_string(df.time_range,tformat=tformat),' ')+tab+$
        strjoin(string(df.theta_range,format='(F7.2)'),' ')+tab+$
        strjoin(string(df.obs_r_sm,format='(F6.2)'),',')+tab
    lprmsg, msg, file
end
