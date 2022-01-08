
pro azim_df_triad_write, triad, filename=file, prefix=prefix
    tab = constant('4space')
    max_probe_length = 5
    if n_elements(prefix) eq 0 then prefix = ''
    msg = tab+extend_string(strjoin(triad.probes,','),length=max_probe_length*3+2)+tab+$
        strjoin(time_string(triad.times,tformat=tformat),',')+tab+$
        strjoin(string(triad.angles,format='(F5.1)'),',')+tab+$
        strjoin(string(triad.center_r_sm,format='(F6.2)'),',')+tab+$
        string(triad.omega_obs_time,format='(F8.2)')+' '+$
        string(triad.omega_xcor,format='(F8.2)')+tab+$
        string(triad.vmag_obs_time,format='(F6.1)')+' '+$
        string(triad.vmag_xcor,format='(F6.1)')+tab+$
        strjoin(string(triad.vhat_obs_time,format='(F5.2)'),',')+' '+$
        strjoin(string(triad.vhat_xcor,format='(F5.2)'),',')+tab+$
        strjoin(string(triad.r_sms[*],format='(F6.2)'),',')
    lprmsg, msg, file
end
