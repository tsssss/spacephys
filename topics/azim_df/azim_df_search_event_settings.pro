
function azim_df_search_event_settings, project=project

    if n_elements(project) eq 0 then project = azim_df_load_project()

;---Some common settings.
    mlt_limit = 9.  ; hour.
    mlt_range = [-1,1]*mlt_limit
    rxy_range = [4.,30] ; Re.
    pdyn = 10.  ; nPa, ~8 Re at nose.

    small_angle_limit = 15.
    triad_angle_range = [0,180]+[1,-1]*small_angle_limit

    search_names = ['beyond_15Re','within_15Re']
    search_settings = list()
    foreach search_name, search_names do begin
        if search_name eq 'beyond_15Re' then begin
            search_time_range = time_double(['2007-03-01','2009-09-01'])
            probes = 'th'+letters('e')
            data_file_suffix = 'azim_df_primitive_data_themis.cdf'
            roi_min_count = 4
            triad_min_count = 2
            roi_probes = !null
        endif else begin
            search_time_range = time_double(['2012-10-01','2017-10-01'])
            probes = ['rbsp'+letters('b'),'th'+['a','d','e'],'mms1','g'+['13','14','15']]   ; turns out no thb/thc candidate.
            data_file_suffix = 'azim_df_primitive_data.cdf'
            roi_min_count = 5
            triad_min_count = 5
            roi_probes = !null
        endelse

        search_setting = dictionary($
            'name', search_name, $
            'probes', probes, $
            'time_range', search_time_range, $
            'data_file_suffix', data_file_suffix)


    ;---Search times within ROI.
        search_roi = dictionary($
            'name', search_name, $
            'mlt_range', mlt_range, $
            'rxy_range', rxy_range, $
            'pdyn', pdyn, $
            'roi_min_count', roi_min_count, $
            'roi_min_duration', constant('secofhour'), $
            'roi_probes', !null, $
            'file_suffix', !null )
        search_setting['search_roi'] = search_roi


    ;---Search times with enough triads.
        search_triad = dictionary($
            'name', search_name, $
            'small_angle_limit', small_angle_limit, $
            'triad_angle_range', triad_angle_range, $
            'triad_min_count', triad_min_count, $
            'triad_min_duration', constant('secofhour'), $
            'file_suffix', !null )
        search_setting['search_triad'] = search_triad


    ;---Search DF within the candidate periods.
        search_df = dictionary($
            'name', search_name, $
            'file_suffix', !null )
        search_setting['search_df'] = search_df


    ;---Search large DF within the candidate times.
        probe_min_count = 4.
        clear_window_size = 60.*10          ; sec.
        height_range = [0.,80]              ; deg.
        width_range = [1.,25]*60.           ; sec.
        min_scaled_height = 8.              ; deg.
        scale_width = 4.                    ; MLT.
        section_time_range = [0.,10]*60.    ; sec.
        section_min_scaled_height = 2.      ; deg.
        section_min_ratio = 0.5             ; #.
        search_large_df = dictionary($
            'name', search_name, $
            'probe_min_count', probe_min_count, $
            'clear_window_size', clear_window_size, $
            'height_range', height_range, $
            'width_range', width_range, $
            'min_scaled_height', min_scaled_height, $
            'scale_width', scale_width, $
            'section_time_range', section_time_range, $
            'section_min_scaled_height', section_min_scaled_height, $
            'section_min_ratio', section_min_ratio, $
            'file_suffix', !null )
        search_setting['search_large_df'] = search_large_df


    ;---Select subgroups.
        min_consecutive_df_count = 4
        search_subgroup = dictionary($
            'name', search_name, $
            'min_consecutive_df_count', min_consecutive_df_count, $

            'file_suffix', !null )
        search_setting['search_subgroup'] = search_subgroup


   ;---Select DF groups.
        search_df_group = dictionary($
            'name', search_name, $
            ; Used to filter dtimes.
            'dtime_min_mean', 8.*60, $          ; sec.
            'dtime_range', [0.,30]*60, $        ; sec.
            ; Used to check distribution.
            'min_abs_mlt', 4., $                ; hr.
            'min_mlt_count', 1, $           	; #.
            ; Used to filter geometry.
            'triad_angle_range', triad_angle_range, $
            'min_triad_count', 2, $
            'triad_dtime_range', [1.,1e9]*60, $  ; sec.

            'file_suffix', !null )
        search_setting['search_df_group'] = search_df_group


    ;---Select uniq DF.
        search_uniq_subgroup = dictionary($
            'name', search_name, $
            'file_suffix', !null )
        search_setting['search_uniq_subgroup'] = search_uniq_subgroup


    ;---Select coherent DF.
        search_coherent_df = dictionary($
            'name', search_name, $
            'file_suffix', !null )
        search_setting['search_coherent_df'] = search_coherent_df


        search_settings.add, search_setting
    endforeach

    return, search_settings
end
