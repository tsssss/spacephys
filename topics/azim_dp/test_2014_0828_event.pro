;+
; Test the new data and algorithm for the 2014_0828 event.
;-

time_range = time_double(['2014-08-28','2014-08-29'])
probes = ['rbsp'+letters('b'),'th'+letters('e'),'g'+['13','14','15','16','17'],'mms1']



time_range = time_double(['2008-03-05','2008-03-06'])
probes = ['th'+letters('e'),'g'+['10','11','12','13','14','15','16']]


;;---Check if theta is consistent with before.
;;   Yes, diff is within the error of float number.
;    time_range = time_double(['2014-08-28/09:30','2014-08-28/11:30'])
;    probe = 'tha'
;
;    prefix = probe+'_'
;    theta_var1 = prefix+'theta'
;    theta_var2 = prefix+'theta_new'
;
;    azim_dp_read_theta, time_range, probe=probe
;    stplot_renew, theta_var1, newname=theta_var2
;
;    if n_elements(project) eq 0 then project = azim_df_load_project()
;    azim_df_load_basic_data, project=project
;
;
;    theta1 = get_var_data(theta_var1, in=time_range)
;    theta2 = get_var_data(theta_var2, in=time_range)
;    print, max(abs(theta1-theta2))
;
;    stop
;
;
;;---Check if the ramp times are the same as before.
;    file = join_path([homedir(),'test.cdf'])
;    file_delete, file, /allow_nonexistent
;    date = time_double('2014-08-28')
;    probes = 'tha'
;    foreach probe, probes do azim_dp_read_ramp_gen_file, date, probe=probe, filename=file


    mlt_range = [0,9]
    event_list = azim_dp_search_event(time_range, probe=probes, mlt_range=mlt_range)



stop

end


;time_range = time_double(['2019-08-01','2019-09-10'])
;probes = ['rbsp'+letters('b'),'th'+letters('e'),'g'+['14','15','16','17'],'mms1']
;foreach probe, probes do azim_dp_read_theta, time_range, probe=probe
