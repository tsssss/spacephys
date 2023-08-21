;+
; We know the mlat of the arc in northern and southern hemispheres.
;-

function check_2013_0501_0850_mapping_uncertainty, event_info=event_info

    if n_elements(event_info) eq 0 then event_info = _2013_0501_0850_load_data()
    test = 0

    probe = event_info['probe']
    prefix = event_info['prefix']

    time_range = event_info['time_range']

    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,'check_2013_0501_0850_mapping_uncertainty.pdf'])

    model_setting = event_info['model_setting']
    models = model_setting['models']
    var = prefix+'fmlt_dipole_north'
    get_data, var, times
    ntime = n_elements(times)
    
    foreach type, ['dipole','igrf'] do begin
        foreach data_type, ['fmlat','fmlt'] do begin
            foreach model, models, model_id do begin
                the_var = prefix+data_type+'_'+type+'_'+model+'_ns'
                the_data = fltarr(ntime,2)
;if ~check_if_update(the_var, time_range) then continue
                foreach var, prefix+data_type+'_'+type+'_'+['north','south'], var_id do begin
                    get_data, var, times, data, limits=lim
                    if data_type eq 'fmlat' then data = abs(data)
                    the_data[*,var_id] = data[*,model_id]
                endforeach
                store_data, the_var, times, the_data, limits={$
                    ytitle: strupcase(model+'+'+type)+'!C'+strupcase(data_type)+'!C('+lim.unit+')', $
                    labels: ['north','south'], $
                    colors: sgcolor(['red','blue']), $
                    ynozero: 1, $
                    labflag: -1 }
            endforeach
        endforeach
    endforeach
    
    
    vars = tnames(prefix+'fmlat*_ns')
    options, vars, 'yrange', [59,65]
    sgopen, plot_file, size=[6,10], test=test
    tplot, vars, trange=time_range
    sgclose


end

print, check_2013_0501_0850_mapping_uncertainty(event_info=event_info)
end