function low_lshell_outflow_load_data, id;, event_info=event_info

;    if n_elements(event_info) ne 0 then if event_info.id eq id then return, event_info

    plot_dir = join_path([googledir(),'works','low_lshell_outflow','plot',id])
    data_dir = join_path([googledir(),'works','low_lshell_outflow','data'])
    version = 'v01'

    if id eq '2015_0317' then begin
        
    ;---Overall setting.
        event_time_range = time_double(['2015-03-17/00:00','2015-03-18/12:00'])
        rbsp_probes = ['a','b']
        themis_probes = ['a','d','e']
        goes_probes = ['13','15']
        mms_probes = ['1','2','3','4']

        event_info = dictionary($
            'time_range', event_time_range, $
            'id', id, $
            'plot_dir', plot_dir, $
            'data_dir', data_dir, $
            'version', version, $
            'rbsp_probes', rbsp_probes, $
            'themis_probes', themis_probes, $
            'goes_probes', goes_probes, $
            'mms_probes', mms_probes, $
            'themis', dictionary(), $
            'rbsp', dictionary(), $
            'dmsp', dictionary(), $
            'goes', dictionary(), $
            'mms', dictionary(), $
            'ground', dictionary() )
    
    ;---Ground settings.
        ground_file = join_path([data_dir,id+'_ground_data_'+version+'.cdf'])
        ground_time_range = time_double(['2015-03-17/06:00','2015-03-17/10:00'])
        event_info.ground_time_range = ground_time_range

        weygand_setting = dictionary($
            'data_file', ground_file, $
            'time_range', ground_time_range )
        j_vars = event_study_read_weygand(weygand_setting, time_var='weygand_time')
        event_info['weygand_setting'] = weygand_setting

        asi_setting = dictionary($
            'data_file', ground_file, $
            'time_range', ground_time_range, $
            'sites', ['inuv','whit','atha','fsim','fsmi','pina','kapu','snkq','gbay','nrsq'], $
            'min_elevs', [10,5,5,10,10,5,5,10,5,5]/2, $
            'merge_method', 'merge_elev', $
            'calibration_method', 'simple' )
        ;mlt_image_var = event_study_read_asi_mlt_image(asi_setting, time_var='asi_time')
        event_info['asi_setting'] = asi_setting


        
    ;---RBSP.
        rbsp_time_range = event_time_range
        rbsp_colors = sgcolor(['magenta','purple'])
        foreach probe, rbsp_probes, probe_id do begin
            rbsp_file = join_path([data_dir,id+'_rbsp'+probe+'_data_'+version+'.cdf'])
            sc_info = low_lshell_outflow_load_rbsp_data($
                rbsp_time_range, filename=rbsp_file, probe=probe, bad_e_time_ranges=bad_e_time_ranges)
            sc_info['sc_name'] = 'RBSP'
            sc_info['sc_color'] = rbsp_colors[probe_id]
            key = 'rbsp'+probe
            event_info.rbsp[key] = sc_info
        endforeach


    ;---THEMIS.
        foreach probe, themis_probes, probe_id do begin
            key = 'th'+probe
            sc_info = dictionary()
            sc_info['sc_name'] = 'TH'
            event_info.themis[key] = sc_info
        endforeach
    endif

    event_info['pp'] = dictionary($
        'ct', 63, $
        'psym', 8, $
        'symsize', 0.75 )

    return, event_info


end


id = '2015_0317'
event_info = low_lshell_outflow_load_data(id)
end