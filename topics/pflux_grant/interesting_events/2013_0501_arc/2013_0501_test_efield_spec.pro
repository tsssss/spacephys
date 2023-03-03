;+
; Test E field freq spectrogram.
;-


    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()

    prefix = event_info['prefix']
    probe = event_info['probe']
    time_range = event_info['time_range']
    snapshot_time = event_info['snapshot_time']

    e_var = prefix+'edot0_fac'
    e_vars = e_var+['b','w','o']
    stplot_split, e_var, newnames=e_vars
    pflux_setting = event_info['pflux_setting']
    scale_info = pflux_setting['scale_info']
    foreach var, e_vars do begin
        stplot_mor, var, scaleinfo=scale_info, frequency=1
    endforeach
    
    mor_vars = e_vars+'_mor'
    options, mor_vars, 'yrange', [1e-3,5]

end
