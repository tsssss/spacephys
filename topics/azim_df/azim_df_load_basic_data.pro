;+
; Load scaled_theta for all probes.
; project=.
; scale_width=.
;-

pro azim_df_load_basic_data_gen_cdf, cdf_file, project=project

;---Check input.
    cdf_id = (file_test(cdf_file) eq 0)? cdf_create(cdf_file): cdf_open(cdf_file)


;---The time range for all probes.
    lprmsg, 'Loading '+project.name+' basic data ...'
    all_probes = list()
    search_types = project.search_types
    foreach search_type, search_types do all_probes.add, search_type.probes, /extract
    all_probes = sort_uniq(all_probes.toarray())

    full_time_ranges = list()
    foreach probe, all_probes do begin
        full_time_range = list()
        foreach search_type, search_types do begin
            index = where(search_type.probes eq probe, count)
            if count eq 0 then continue
            full_time_range.add, search_type.time_range
        endforeach
        full_time_range = minmax(full_time_range.toarray())
        full_time_ranges.add, full_time_range
    endforeach


;---Dst/AE.
    ; This is pretty slow, so avoid reload.
    load_omni = 0
    foreach var, ['omni_ut','dst','ae'] do begin
        if ~cdf_has_var(var, filename=cdf_id) then begin
            load_omni = 1
            break
        endif
    endforeach

    if load_omni then begin
        times = []
        ae = []
        dst = []
        foreach search_type, search_types do begin
            time_range = search_type.time_range
            omni_read_index, time_range
            ae = [ae,get_var_data('ae',times=uts)]
            dst = [dst,get_var_data('dst')]
            times = [times,uts]
        endforeach
        index = sort(times)
        times = times[index]
        ae = float(ae[index])
        dst = float(dst[index])

        time_var = 'omni_ut'
        if ~cdf_has_var(time_var, filename=cdf_id) then begin
            cdf_save_var, time_var, value=times, filename=cdf_id
        endif

        the_var = 'dst'
        if ~cdf_has_var(the_var, filename=cdf_id) then begin
            settings = dictionary($
                'depend_0', time_var, $
                'display_type', 'scalar', $
                'unit', 'nT', $
                'short_name', 'Dst')
            cdf_save_var, the_var, value=dst, filename=cdf_id
            cdf_save_setting, settings, filename=cdf_id, varname=the_var
        endif

        the_var = 'ae'
        if ~cdf_has_var(the_var, filename=cdf_id) then begin
            settings = dictionary($
                'depend_0', time_var, $
                'display_type', 'scalar', $
                'unit', 'nT', $
                'short_name', 'AE')
            cdf_save_var, the_var, value=ae, filename=cdf_id
            cdf_save_setting, settings, filename=cdf_id, varname=the_var
        endif
    endif


;---Load data for each probe.
    foreach probe, all_probes, probe_id do begin
        prefix = probe+'_'
        full_time_range = full_time_ranges[probe_id]
        ;msg = 'Processing '+probe+': '+strjoin(time_string(full_time_range),' to ')+' ...'
        ;lprmsg, msg

    ;---theta.
        theta_var = prefix+'theta'
        if ~cdf_has_var(theta_var, filename=cdf_id) then begin
            times = []
            theta = []
            foreach search_type, search_types do begin
                index = where(search_type.probes eq probe, count)
                if count eq 0 then continue
                data_file_suffix = search_type.data_file_suffix
                data_file = join_path([project.data_dir,data_file_suffix])
                theta = [theta,cdf_read_var(theta_var, filename=data_file)]
                times = [times,cdf_read_var('bfield_ut', filename=data_file)]
            endforeach
            index = sort(times)
            times = times[index]
            theta = theta[index]
            store_data, theta_var, times, theta
        endif

    ;---r_sm.
        r_sm_var = prefix+'r_sm'
        if ~cdf_has_var(r_sm_var, filename=cdf_id) then begin
            times = []
            r_sms = []
            foreach search_type, search_types do begin
                index = where(search_type.probes eq probe, count)
                if count eq 0 then continue
                data_file_suffix = search_type.data_file_suffix
                data_file = join_path([project.data_dir,data_file_suffix])
                r_sms = [r_sms,cdf_read_var(r_sm_var, filename=data_file)]
                times = [times,cdf_read_var('orbit_ut', filename=data_file)]
            endforeach
            index = sort(times)
            times = times[index]
            r_sms = r_sms[index,*]
            store_data, r_sm_var, times, r_sms
        endif
    endforeach


;---Save to file.
    foreach probe, all_probes, probe_id do begin
        prefix = probe+'_'

    ;---Theta.
        theta_var = prefix+'theta'
        get_data, theta_var, times, theta
        theta = float(theta)
        
        time_var = prefix+'ut'
        if ~cdf_has_var(time_var, filename=cdf_id) then begin
            cdf_save_var, time_var, value=times, filename=cdf_id
        endif

        the_var = theta_var
        settings = dictionary($
            'depend_0', time_var, $
            'display_type', 'scalar', $
            'unit', 'deg', $
            'short_name', tex2str('theta'))
        cdf_save_var, the_var, value=theta, filename=cdf_id
        cdf_save_setting, settings, filename=cdf_id, varname=the_var


    ;---R SM.
        r_sm_var = prefix+'r_sm'
        get_data, r_sm_var, times, r_sms
        r_sms = float(r_sms)

        orbit_time_var = prefix+'orbit_ut'
        if ~cdf_has_var(orbit_time_var, filename=cdf_id) then begin
            cdf_save_var, orbit_time_var, value=times, filename=cdf_id
        endif

        the_var = r_sm_var
        settings = dictionary($
            'depend_0', orbit_time_var, $
            'display_type', 'vector', $
            'unit', 'Re', $
            'coord', 'SM', $
            'coord_labels', constant('xyz'), $
            'short_name', 'R')
        cdf_save_var, the_var, value=r_sms, filename=cdf_id
        cdf_save_setting, settings, filename=cdf_id, varname=the_var


    ;---mlt.
        mlt_var = prefix+'pseudo_mlt'
        if ~cdf_has_var(mlt_var, filename=cdf_id) then begin
            times = cdf_read_var(orbit_time_var, filename=cdf_id)
            r_sms = cdf_read_var(r_sm_var, filename=cdf_id)
            mlt = pseudo_mlt(r_sms)
            mlt = float(mlt)

            the_var = mlt_var
            settings = dictionary($
                'depend_0', orbit_time_var, $
                'display_type', 'scalar', $
                'unit', 'hr', $
                'short_name', 'MLT')
            cdf_save_var, the_var, value=mlt, filename=cdf_id
            cdf_save_setting, settings, filename=cdf_id, varname=the_var
        endif


    ;---scaled_theta.
        scale_width = project.scale_width
        scaled_theta_var = prefix+'scaled_theta'
        if ~cdf_has_var(scaled_theta_var, filename=cdf_id) then begin
            times = cdf_read_var(time_var, filename=cdf_id)
            theta = cdf_read_var(theta_var, filename=cdf_id)
            orbit_times = cdf_read_var(orbit_time_var, filename=cdf_id)
            mlt = cdf_read_var(mlt_var, filename=cdf_id)
            mlt = interpol(mlt, orbit_times, times)
            scaled_theta = azim_df_scale_theta(theta, mlt, width=scale_width)
            scaled_theta = float(scaled_theta)

            the_var = scaled_theta_var
            settings = dictionary($
                'depend_0', time_var, $
                'display_type', 'scalar', $
                'unit', 'deg', $
                'short_name', 'scaled '+tex2str('theta'))
            cdf_save_var, the_var, value=scaled_theta, filename=cdf_id
            cdf_save_setting, settings, filename=cdf_id, varname=the_var
        endif
    endforeach

    cdf_close, cdf_id


end

pro azim_df_load_basic_data, project=project, scale_width=scale_width

    if n_elements(project) eq 0 then project = azim_df_load_project()
    search_types = azim_df_search_event_settings(project=project)
    status_var = 'azim_df_have_basic_data'
    if tnames(status_var) ne '' then begin
        if get_var_data(status_var) then return
    endif

    file_suffix = project.name+'_basic_data.cdf'
    cdf_file = join_path([project.data_dir,file_suffix])
    if file_test(cdf_file) eq 0 then azim_df_load_basic_data_gen_cdf, cdf_file, project=project

    foreach var, ['dst','ae'] do cdf_load_var, var, filename=cdf_file
    all_probes = project.all_probes
    vars = ['scaled_theta','theta','r_sm','pseudo_mlt']
    foreach probe, all_probes do begin
        prefix = probe+'_'
        foreach var, vars do cdf_load_var, prefix+var, filename=cdf_file
    endforeach


;---To avoid loading data again.
    store_data, status_var, 0, 1
end

azim_df_load_basic_data
end