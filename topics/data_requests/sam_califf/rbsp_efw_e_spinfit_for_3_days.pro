

; The wanted days.
probes = ['a','b']
event_list = list()
event_list.add, dictionary($
    'time_range', time_double(['2013-03-17/09:24','2013-03-17/10:53']), $
    'probe', 'b')
event_list.add, dictionary($
    'time_range', time_double(['2013-03-17/10:39','2013-03-17/12:04']), $
    'probe', 'a')
event_list.add, dictionary($
    'time_range', time_double(['2015-07-13/03:07','2015-07-13/04:38']), $
    'probe', 'b')
event_list.add, dictionary($
    'time_range', time_double(['2015-07-13/07:49','2015-07-13/09:36']), $
    'probe', 'a')
event_list.add, dictionary($
    'time_range', time_double(['2015-09-07/19:32','2015-09-07/21:01']), $
    'probe', 'a')
event_list.add, dictionary($
    'time_range', time_double(['2015-09-07/22:37','2015-09-07/23:55']), $
    'probe', 'b')


secofday = 86400d
root_dir = join_path([rbsp_efw_phasef_local_root()])


foreach event, event_list do begin
    probe = event.probe
    time_range = event.time_range
    day = time_range[0]-(time_range[0] mod secofday)

    rbspx = 'rbsp'+probe
    prefix = 'rbsp'+probe+'_'
    year = time_string(day,tformat='YYYY')
    base = rbspx+'_efw-l3_'+time_string(day,tformat='YYYYMMDD')+'_v04.cdf'
    file = join_path([root_dir,rbspx,'l3',year,base])

    ;rbsp_efw_phasef_read_density, time_range, probe=probe
    ;rbsp_efw_phasef_read_l3, time_range, probe=probe
    rbsp_efw_phasef_read_e_uvw, time_range, probe=probe
    cdf2tplot, file
    rbsp_efw_phasef_read_b_mgse, time_range, probe=probe
    get_data, prefix+'b_mgse', times, b_mgse
    store_data, 'spinfit_angle', times, acos(b_mgse[*,0]/snorm(b_mgse))*constant('deg'), $
        limits={constant: 90+[-1,1]*15}
    
    tplot, [prefix+'e_uvw','efield_in_corotation_frame_spinfit_'+['mgse','edotb_mgse'],$
        'spacecraft_potential','density','lshell','spinfit_angle'], trange=time_range
    stop
endforeach
end
