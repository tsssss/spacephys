function alfven_arc_load_data, id, event_info=event_info

    if n_elements(id) eq 0 then id = event_info.id
    if n_elements(event_info) ne 0 then if even_info.id eq id then return, event_info

    plot_dir = join_path([googledir(),'works','pflux_grant','alfven_arc','plot',id])
    data_dir = join_path([googledir(),'works','pflux_grant','alfven_arc','data'])
    version = 'v01'

    if id eq '2017_0309_0700' then begin
        version = 'v01'

        themis_time_range = time_double(['2017-03-09/06:30','2017-03-09/09:00'])
        ground_time_range = themis_time_range
        themis_probes = ['d','e']
        themis_colors = sgcolor(['magenta','purple'])
        dmsp_probes = 'f'+['18']
        
        event_info = dictionary($
            'time_range', themis_time_range, $
            'id', id, $
            'plot_dir', plot_dir, $
            'data_dir', data_dir, $
            'version', version, $
            'external_model', 't89', $
            'internal_model', 'dipole', $
            'themis', dictionary(), $
            'rbsp', dictionary(), $
            'dmsp', dictionary(), $
            'ground', dictionary() )

    ;---Ground settings.
        ground_file = join_path([data_dir,id+'_ground_data_'+version+'.cdf'])
        asi_setting = dictionary($
            'sites', ['fsim','fsmi','gill','fykn'], $
            'best_site', 'fsim', $
            'min_elevs', [5,10,2.5,2.5], $
            'merge_method', 'merge_elev', $
            'calibration_method', 'moon', $
            'mlt_range', [-1,1]*6, $
            'mlat_range', [60d,75] )
        event_info.ground = alfven_arc_load_ground_data($
            filename=ground_file, $
            ground_time_range, asi_setting=asi_setting)

        
    ;---Themis settings.
        foreach probe, themis_probes, probe_id do begin
            bad_e_time_ranges = !null
            themis_file = join_path([data_dir,id+'_th'+probe+'_data_'+version+'.cdf'])
            event_info.themis['th'+probe] = alfven_arc_load_themis_data($
                filename=themis_file, themis_time_range, probe=probe, bad_e_time_ranges=bad_e_time_ranges)
            (event_info.themis['th'+probe])['sc_name'] = 'TH'
            (event_info.themis['th'+probe])['sc_color'] = themis_colors[probe_id]
        endforeach

    
    
    ;---DMSP settings.
        foreach probe, dmsp_probes do begin
            event_info.dmsp['dmsp'+probe] = dictionary($
                'sc_name', 'DMSP', $
                'probe', probe, $
                'prefix', 'dmsp'+probe+'_', $
                'sc_color', sgcolor('teal') )
        endforeach
    endif

    if id eq '2015_0302_1100' then begin
        version = 'v01'

        rbsp_time_range = time_double(['2015-03-02/10:30','2015-03-02/11:30'])
        ground_time_range = rbsp_time_range
        rbsp_probes = ['a']
        dmsp_probes = ['f19']

        event_info = dictionary($
            'time_range', rbsp_time_range, $
            'id', id, $
            'plot_dir', plot_dir, $
            'data_dir', data_dir, $
            'version', version, $
            'external_model', 't04s', $
            'internal_model', 'dipole', $
            'themis', dictionary(), $
            'rbsp', dictionary(), $
            'dmsp', dictionary(), $
            'ground', dictionary() )
        
    ;---Ground settings.
        ground_file = join_path([data_dir,id+'_ground_data_'+version+'.cdf'])
        asi_setting = dictionary($
            'sites', ['mcgr','fykn','inuv'], $
            'best_site', 'fykn', $
            'min_elevs', float([1,1,1]*2), $
            'merge_method', 'merge_elev', $
            'calibration_method', 'moon', $
            'mlt_range', [-3d,3], $
            'mlat_range', [55d,75] )
        event_info.ground = alfven_arc_load_ground_data($
            filename=ground_file, $
            ground_time_range, asi_setting=asi_setting)

    ;---RBSP settings.
        foreach probe, rbsp_probes do begin
            bad_e_time_ranges = !null
            if probe eq 'a' then begin
                bad_e_time_ranges = list($
                    ['2015-03-02/11:14','2015-03-02/11:22'] )
            endif
            rbsp_file = join_path([data_dir,id+'_rbsp'+probe+'_data_'+version+'.cdf'])
            event_info.rbsp['rbsp'+probe] = alfven_arc_load_rbsp_data($
                filename=rbsp_file, $
                rbsp_time_range, probe=probe, bad_e_time_ranges=bad_e_time_ranges)
        endforeach
    
    
    ;---DMSP settings.
        foreach probe, dmsp_probes do begin
            event_info.dmsp['dmsp'+probe] = dictionary($
                'probe', probe, $
                'prefix', 'dmsp'+probe+'_' )
        endforeach
    endif
    
    if id eq '2015_0105_0100' then begin
        rbsp_time_range = time_double(['2015-01-05/00:00','2015-01-05/02:00'])
        ground_time_range = rbsp_time_range
        rbsp_probes = ['a']
        dmsp_probes = ['f17']
        event_info = dictionary($
            'time_range', rbsp_time_range, $
            'id', id, $
            'plot_dir', plot_dir, $
            'data_dir', data_dir, $
            'version', version, $
            'external_model', 't89', $
            'internal_model', 'dipole', $
            'rbsp', dictionary(), $
            'dmsp', dictionary(), $
            'ground', dictionary() )
        
    ;---Ground settings.
        ground_file = join_path([data_dir,id+'_ground_data_'+version+'.cdf'])
        asi_setting = dictionary($
            'sites', ['nrsq','gill','rank','snkq'], $
            'best_site', 'nrsq', $
            'min_elevs', float([2.5,2.5,2.5,2.5]), $
            'merge_method', 'merge_elev', $
            'calibration_method', 'moon', $
            'mlt_range', [-8d,2], $
            'mlat_range', [60d,75] )
        event_info.ground = alfven_arc_load_ground_data($
            filename=ground_file, $
            ground_time_range, asi_setting=asi_setting)

    ;---RBSP settings.
        foreach probe, rbsp_probes do begin
            bad_e_time_ranges = !null
            if probe eq 'a' then begin
                bad_e_time_ranges = list($
                    ['2015-01-05/00:00:00','2015-01-05/00:02'], $
                    ['2015-01-05/00:05:00','2015-01-05/00:10'], $
                    ['2015-01-05/00:13','2015-01-05/00:14:30'])
            endif
            rbsp_file = join_path([data_dir,id+'_rbsp'+probe+'_data_'+version+'.cdf'])
            event_info.rbsp['rbsp'+probe] = alfven_arc_load_rbsp_data($
                filename=rbsp_file, $
                rbsp_time_range, probe=probe, bad_e_time_ranges=bad_e_time_ranges)
        endforeach
    
    
    ;---DMSP settings.
        foreach probe, dmsp_probes do begin
            event_info.dmsp['dmsp'+probe] = dictionary($
                'probe', probe, $
                'prefix', 'dmsp'+probe+'_' )
        endforeach
    endif

    if id eq '2015_0317_0900' then begin
        time_range = time_double(['2015-03-17/08:20','2015-03-17/09:20'])
        event_info = dictionary($
            'time_range', time_range, $
            'id', id, $
            'plot_dir', plot_dir, $
            'external_model', 't89', $
            'internal_model', 'dipole', $
            'rbsp', dictionary(), $
            'dmsp', dictionary(), $
            'ground', dictionary() )
        
        asi_setting = dictionary($
            'sites', ['whit','atha'], $
            'best_site', 'whit', $
            'min_elevs', float([1,2.5]), $
            'merge_method', 'max_elev', $
            'calibration_method', 'simple', $
            'mlt_range', [-6d,6], $
            'mlat_range', [55d,70] )
        event_info.ground = alfven_arc_load_ground_data(time_range, asi_setting=asi_setting)

        rbsp_probes = ['a']
        foreach probe, rbsp_probes do begin
            bad_e_time_ranges = !null
            if probe eq 'a' then begin
                bad_e_time_ranges = list(['2015-03-17/09:05:50','2015-03-17/09:41'])
            endif
            event_info.rbsp['rbsp'+probe] = alfven_arc_load_rbsp_data(time_range, probe=probe, bad_e_time_ranges=bad_e_time_ranges)
        endforeach
    endif


    if id eq '2015_0312_0900' then begin
        event_info = dictionary($
            'time_range', time_double(['2015-03-12/08:00','2015-03-12/11:00']), $
            'id', id, $
            'plot_dir', plot_dir, $
            'external_model', 't89', $
            'internal_model', 'dipole', $
            'rbsp', dictionary(), $
            'dmsp', dictionary(), $
            'ground', dictionary() )
        time_range = event_info.time_range
        
        asi_setting = dictionary($
            'sites', ['fykn','mcgr'], $
            'best_site', 'fykn', $
            'min_elevs', float([5,10]), $
            'merge_method', 'max_elev', $
            'calibration_method', 'simple', $
            'mlt_range', [-6d,0], $
            'mlat_range', [55d,75] )
        event_info.ground = alfven_arc_load_ground_data(time_range, asi_setting=asi_setting)

        rbsp_probes = ['b']
        foreach probe, rbsp_probes do begin
            if probe eq 'b' then bad_e_time_ranges = (probe eq 'b')? list(['2015-03-12/09:04:51','2015-03-12/09:04:54']): !null
            event_info.rbsp['rbsp'+probe] = alfven_arc_load_rbsp_data(time_range, probe=probe, bad_e_time_ranges=bad_e_time_ranges)
        endforeach
        
        ;---DMSP settings.
        dmsp_probes = ['f19']
        foreach probe, dmsp_probes do begin
            event_info.dmsp['dmsp'+probe] = dictionary($
                'probe', probe, $
                'prefix', 'dmsp'+probe+'_' )
        endforeach
    endif

    return, event_info

end


event_info = alfven_arc_load_data('2017_0309_0700')
stop
event_info = alfven_arc_load_data('2015_0302_1100')
stop
event_info = alfven_arc_load_data('2015_0105_0100')
stop
event_info = alfven_arc_load_data('2015_0317_0900')
end