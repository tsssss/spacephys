;+
; Load storm list, and orbit data.
; Pick out storms when enough # of sc in the desired region.
;-

pro azim_df_filter_storm, project=project


;---Prepare settings.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if ~project.haskey('done_find_storm') then azim_df_find_storm, project=project
    storm_times = project.storm_times
    storm_ids = project.storm_id

    mlt_limit = 7.  ; hr away from midnight.
    settings = dictionary($
        'mission_probes', ['rbsp'+letters('b'),'th'+['a','d','e']], $
        'all_mission_probes', ['rbsp'+letters('b'),'th'+letters('e')], $
        'data_file_suffix', 'apogee_history.tplot', $
        'search_time_range', time_double(['2013-01-01','2018-01-01']), $
        'mlt_limit', mlt_limit, $
        'mlt_range', [-1,1]*mlt_limit, $
        'file_suffix', 'filter_storm.txt')
    file = join_path([project.data_dir,settings.file_suffix])
    log_file = join_path([project.data_dir,settings.file_suffix+'.log'])
    if file_test(log_file) ne 0 then file_delete, log_file
    ftouch, log_file
    if keyword_set(test) then log_file = -1

;---We first get an apogee history of all spacecraft.
    full_time_range = []
    foreach time_range, storm_times do full_time_range = minmax([full_time_range, time_range])
    epoch = stoepoch(full_time_range, 'unix')
    epoch = [sepochfloor(epoch[0],'mo'),sepochceil(epoch[1],'mo')]
    full_time_range = sfmepoch(epoch, 'unix')

    apogee_file = join_path([project.data_dir,settings.data_file_suffix])
    all_mission_probes = settings.all_mission_probes
    if file_test(apogee_file) eq 1 then begin
        tplot_restore, filename=apogee_file
    endif else begin
        the_vars = list()
        foreach mission_probe, all_mission_probes do begin
            mission_info = resolve_probe(mission_probe)
            routine_name = mission_info.routine_name
            probe = mission_info.probe
            prefix = mission_info.prefix
            rvar = prefix+'r_gsm'

            data_file = join_path([googledir(),'works','works','global_efield','data','combined_data',$
                'global_efield_'+mission_probe+'_combined_data_v01.cdf'])
            cdf_load_var, rvar, filename=data_file, time_var='ut_sec'
            ; Should do this but it's faster to use what I already got.
            ; call_procedure, routine_name+'_read_orbit', full_time_range, probe=probe

            find_apogee, rvar, apogee_times=times
            rgsm = get_var_data(rvar, at=times)
            the_var = prefix+'apogee_history'
            store_data, the_var, times, rgsm
            add_setting, the_var, /smart, {$
                display_type: 'vector', $
                unit: 'Re', $
                short_name: 'R', $
                coord: 'GSM', $
                coord_labels: ['x','y','z']}
            the_vars.add, the_var
        endforeach
        tplot_save, the_vars.toarray(), filename=apogee_file
    endelse


;---Remove data outside the magnetopause.
    all_prefixes = strarr(n_elements(all_mission_probes))
    foreach mission_probe, all_mission_probes, ii do begin
        mission_info = resolve_probe(mission_probe)
        all_prefixes[ii] = mission_info.prefix
    endforeach
    deg = 180d/!dpi
    foreach prefix, all_prefixes do begin
        the_var = prefix+'apogee_history'
        get_data, the_var, times, rgsm
        flags = check_if_in_magn(rgsm)
        index = where(flags eq 0, count)
        if count ne 0 then begin
            rgsm[index,*] = !values.f_nan
            store_data, the_var, times, rgsm
        endif
        rmag = cotran(rgsm, times, 'gsm2mag')
        mlon = atan(rmag[*,1],rmag[*,0])*deg
        mlt = mlon2mlt(mlon, times)
        store_data, prefix+'apogee_mlt', times, mlt
    endforeach


;---Filter out storms.
    search_time_range = settings.search_time_range
    lprmsg, 'Search for storms from '+strjoin(time_string(search_time_range),' to ')+' ...', log_file
    mission_probes = settings.mission_probes
    lprmsg, 'Search for spacecraft among '+strjoin(strupcase(mission_probes),',')+' ...', log_file
    lprmsg, 'Search for all apogees on one side of midnight, and within '+sgnum2str(mlt_limit)+' hours from midnight ...', log_file
    prefixes = strarr(n_elements(mission_probes))
    foreach mission_probe, mission_probes, ii do begin
        mission_info = resolve_probe(mission_probe)
        prefixes[ii] = mission_info.prefix
    endforeach


    candidate = hash()
    nprefix = n_elements(prefixes)
    foreach time_range, storm_times do begin
        lprmsg, '', log_file
        lprmsg, 'Processing storm from '+strjoin(time_string(time_range),' to ')+' ...', log_file

    ;---Throw away storms out of the search time range.
        mid_time = mean(time_range)
        index = lazy_where(mid_time, '][', search_time_range, count=count)
        if count ne 0 then begin
            lprmsg, 'Outside search time ...', log_file
            continue
        endif

    ;---Throw away storms when some probe is around noon.
        mlt_range = settings.mlt_range
        post_mlt_flags = bytarr(nprefix)
        pre_mlt_flags = bytarr(nprefix)
        foreach prefix, prefixes, ii do begin
            the_var = prefix+'apogee_mlt'
            sc = strupcase(strmid(prefix,0,strpos(prefix,'_')))
            mlt = get_var_data(the_var, in=time_range)
            index = where(finite(mlt), count)
            if count eq 0 then begin
                lprmsg, sc+' outside magnetopause ...', log_file
                continue
            endif else mlt = mlt[index]
            index = lazy_where(mlt, '][', mlt_range, count=count)
            if count ne 0 then begin
                lprmsg, sc+' round noon ...', log_file
                continue
            endif
            index = where(mlt le 0, count)
            pre_mlt_flags[ii] = count ne 0
            index = where(mlt ge 0, count)
            post_mlt_flags[ii] = count ne 0
        endforeach

        the_flag = 0
        pre_mlt_index = where(pre_mlt_flags eq 1, pre_mlt_count)
        lprmsg, sgnum2str(pre_mlt_count)+' spacecraft in pre-midnight ...', log_file
        post_mlt_index = where(post_mlt_flags eq 1, post_mlt_count)
        lprmsg, sgnum2str(post_mlt_count)+' spacecraft in post-midnight ...', log_file
        if max([pre_mlt_count,post_mlt_count]) eq nprefix then begin
            candidate_id = time_string(mean(time_range),tformat='YYYY_MMDD_hh')
            candidate[candidate_id] = dictionary($
                'id', candidate_id, $
                'time_range', time_range)
            lprmsg, 'Found a candidate ...', log_file
            continue
        endif
    endforeach


;---Save the result to project.
    lprmsg, 'Save result to project.candidate ...'
    project['candidate_storm'] = candidate
    project['filter_storm_setting'] = settings
    project['done_filter_storm'] = 1
    update_project, project


;---Output to a text file.
    lprmsg, 'Save result to '+file+' ...'
    openw, lun, file, /get_lun
    keys = candidate.keys()
    keys = keys.sort()
    foreach key, keys do begin
        tinfo = candidate[key]
        printf, lun, tinfo.id+'    '+strjoin(time_string(tinfo.time_range), ' to ')
    endforeach
    free_lun, lun
end

azim_df_filter_storm
end
