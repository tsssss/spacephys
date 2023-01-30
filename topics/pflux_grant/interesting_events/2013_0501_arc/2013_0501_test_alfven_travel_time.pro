;+
; Trace Alfven waves to ionosphere.
;-


    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()

;---Settings.
    time_range = time_double(['2013-05-01/07:20','2013-05-01/07:50'])
    probe = event_info['probe']
    prefix = event_info['prefix']
    model_time = event_info['snapshot_time']
    model_setting = event_info['model_setting']
    model = model_setting['external_model']
    igrf = model_setting['igrf']
    plasma_param = event_info['plasma_param']
    stop_dis = 0


    ; Trace to ionosphere.
    trace_input_list = list()
    trace_input_list.add, dictionary($
        'time', model_time, $
        'stop_dis', stop_dis, $
        'mod_time', model_time, $
        'model', model, $
        'ion_mass', plasma_param['avg_ion_mass'], $
        'igrf', igrf )

    r_gsm_var = prefix+'r_gsm'
    trace_output_list = list()
    foreach info, trace_input_list do begin        
        time = time_double(info['time'])
        info['r_gsm'] = get_var_data(r_gsm_var, at=time)
        info['time'] = time
        par_var = info['model']+'_var'
        info['par'] = reform(get_var_data(par_var, at=time))
    
        trace_output = trace_alfven_wave_to_ionosphere(time, _extra=info.tostruct())
        trace_output_list.add, trace_output
        stop
    endforeach
    

end