
pro azim_df_edge_write, edge, filename=file, prefix=prefix
    tab = constant('4space')
    max_probe_length = 5
    if n_elements(prefix) eq 0 then prefix = ''
    msg = prefix[0]+extend_string(strjoin(edge.probes,','),length=max_probe_length*2+2)+tab+$
        string(edge.time_lag_df,format='(F9.2)')+tab+$
        string(edge.time_lag_xcor,format='(F9.2)')+tab+$
        string(edge.xcor_max,format='(F5.2)')+tab+$
        string(edge.xcor_err,format='(F5.1)')+tab+$
        string(edge.time_lag_diff,format='(F7.1)')+tab+$
        string(edge.time_lag_ratio,format='(F5.2)')
    lprmsg, msg, file
end
