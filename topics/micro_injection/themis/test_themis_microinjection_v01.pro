test = 1
high_res = 1

event_list = list()
;event_list.add, dictionary($
;    'tr', time_double('2010-09-08/'+['00:00','24:00']), $
;    'probe', 'e' )

;    event_list.add, dictionary($
;;        'tr', time_double('2010-09-08/'+['17:10','18:20']), $
;        'tr', time_double('2010-09-08/'+['16:00','19:00']), $
;        'probe', 'e' )
event_list.add, dictionary($
    'tr', time_double('2010-09-06/'+['21:00','22:30']), $
    'probe', 'e' )
event_list.add, dictionary($
    'tr', time_double(['2024-03-05/06:00','2024-03-05/16:00']), $
    'probe', 'e' )
event_list.add, dictionary($
    'tr', time_double(['2024-03-07/18:00','2024-03-07/23:00']), $
    'probe', 'e' )
plot_dir = join_path([googledir(),'works','2024_microinjection','plot','themis'])

foreach event, event_list do begin
    time_range = event['tr']
    probe = event['probe']
    prefix = 'th'+probe+'_'
    print, 'Processing '+probe+' for '+time_string(time_range[0])+' ...'
    
    if keyword_set(high_res) then begin
        en_high_var = themis_read_en_spec_esa_sst_integrate(time_range, probe=probe, $
            species='e', id='sst')
        en_high_var = prefix+'kev_e_flux'
        
;        time_range = time_double('2010-09-08/'+['16:00','19:00'])
;        probe = 'e'
        
        timespan, time_range[0], total(time_range*[-1,1]), second=1
        datatype = 'psef'

;        thm_sst_load_calibrate,probe=probe,datatype=datatype,trange=time_range,dist_data=dist
        ;thm_part_moments,inst=datatype,probe=probe,dist_array=dist
;        thm_part_getspec,data_type=datatype,probe=probe,dist_array=dist,outputs='pa gyro'
        thm_part_load, data_type=datatype, probe=probe, trange=time_range
        thm_part_getspec, data_type=datatype, probe=probe, trange=time_range, outputs='pa'
        tplot,prefix+'psef_'+['eflux_pa','en_eflux']
        en_high_var = prefix+datatype+'_flux_energy'
        pa_high_var = prefix+datatype+'_an_eflux_pa'
        options, pa_high_var, zrange=[1e3,1e6]
    endif else begin
        en_high_var = themis_read_kev_electron(time_range, probe=probe, spec=1)
    endelse

    options, en_high_var, color_table=40, zrange=[1e-3,1e1], no_interp=1
    b_var = themis_read_bfield(time_range, probe=probe)
    u_var = themis_read_ion_vel(time_range, probe=probe)
    en_low_var = themis_read_en_spec(time_range, probe=probe, species='e', id='esa_l2')
;    en_low_var = themis_read_en_spec(time_range, probe=probe, species='e')
    options, en_low_var, color_table=40, no_interp=1
    plot_vars = [en_high_var,en_low_var,b_var,u_var]
    if keyword_set(high_res) then plot_vars = [en_high_var,pa_high_var,en_low_var,b_var,u_var]

    plot_file = join_path([plot_dir,'microinjection_test_plot_th'+probe+'_'+time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_v01.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, size=[12,6]
    tplot, plot_vars, trange=time_range
    if keyword_set(test) then stop
    sgclose
endforeach
end