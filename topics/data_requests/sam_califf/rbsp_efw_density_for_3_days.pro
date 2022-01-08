

; The wanted days.
days = time_double(['2013-03-17', '2015-07-13', '2015-09-07', '2015-09-08'])
probes = ['a','b']

secofday = 86400d
root_dir = join_path([diskdir('Research'),'data','rbsp'])

target_dir = join_path([srootdir(),'data'])
if file_test(target_dir) eq 0 then file_mkdir, target_dir

foreach day, days do begin
    foreach probe, probes do begin
        rbspx = 'rbsp'+probe
        time_range = day+[0,secofday]
        year = time_string(day,tformat='YYYY')
        base = rbspx+'_efw-l3_'+time_string(day,tformat='YYYYMMDD')+'_v03.cdf'
        file = join_path([root_dir,rbspx,'l3',year,base])
        target_file = join_path([target_dir,base])
        file_copy, file, target_file, /overwrite
        ;rbsp_efw_phasef_read_density, time_range, probe=probe
        ;rbsp_efw_phasef_read_l3, time_range, probe=probe
        cdf2tplot, target_file
        tplot, ['efield_in_corotation_frame_spinfit_'+['mgse','edotb_mgse'],$
            'spacecraft_potential','density','position_gse'], trange=time_range
        stop
    endforeach
endforeach
end