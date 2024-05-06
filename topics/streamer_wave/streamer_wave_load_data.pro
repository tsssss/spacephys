function streamer_wave_load_data, id, version=version

    plot_dir = join_path([googledir(),'works','streamer_wave','plot',id])
    data_dir = join_path([googledir(),'works','streamer_wave','data'])
    if n_elements(version) eq 0 then version = 'v01'

    if id eq '2008_0328' then begin
        event_time_range = time_double(['2008-03-28/06:25','2008-03-28/06:45'])
        themis_probes = ['d','e']

        event_info = dictionary($
            'time_range', event_time_range, $
            'id', id, $
            'plot_dir', plot_dir, $
            'data_dir', data_dir, $
            'version', version, $
            'themis_probes', themis_probes, $
            'themis', dictionary(), $
            'ground', dictionary() )
        
    ;---Ground settings.
        ground_file = join_path([data_dir,id+'ground_data_'+version+'.cdf'])
        ground_time_range = event_time_range
        event_info.ground_time_range = ground_time_range

        weygand_setting = dictionary($
            'data_file', ground_file, $
            'time_range', ground_time_range )
        event_info['weygand_setting'] = weygand_setting
        j_vars = event_study_read_weygand(weygand_setting, time_var='weygand_time')

        asi_sites = ['fsim','fsmi']
        asi_elevs = !null
        mlat_range = [65d,75]
        asi_setting = dictionary($
            'data_file', ground_file, $
            'time_range', ground_time_range, $
            'sites', asi_sites, $
            'min_elevs', asi_elevs, $
            'merge_method', 'merge_elev', $
            'calibration_method', 'simple' )
        event_info['asi_setting'] = asi_setting
        mlt_image_var = themis_asf_read_mlt_image_rect(ground_time_range, site=asi_sites, mlat_range=mlat_range)

    ;---Themis.
        themis_time_range = event_time_range+[-1,1]*600
        themis_colors = sgcolor(['red','green'])
        foreach probe, themis_probes, probe_id do begin
            data_file = join_path([data_dir,id+'_themis_'+probe+'_data_'+version+'.cdf'])

;            sc_info = streamer_wave_load_data_themis($
;                themis_time_range, filename=data_file, probe=probe)
            sc_info = alfven_arc_load_themis_data($
                filename=data_file, themis_time_range, probe=probe, bad_e_time_ranges=!null)
            sc_info['sc_name'] = 'TH'
            sc_info['sc_color'] = themis_colors[probe_id]

            mission_probe = 'th'+probe
            event_info.themis[mission_probe] = sc_info
        endforeach

    endif

    return, event_info

end

id = '2008_0328'
event_info = streamer_wave_load_data(id)
end