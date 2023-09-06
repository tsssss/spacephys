;+
; Get the middle day of a storm, get the apogee MLT for each probe, only include those storms when:
;   1. THEMIS and RBSP apogee MLTs are away from noon by 3 hours;
;   2. THEMIS and RBSP apogee MLTs are in the same side of the midnight.
;-

pro azim_prop_filter_by_apogee, project

    ; Get the storm list.
    if n_elements(project) eq 0 then project = azim_prop_load_project()
    if ~project.haskey('storm_times') then azim_prop_find_storm
    storm_times = project.storm_times

    ; Settings.
    apogee_filter_setting = dictionary()
    apogee_filter_setting.probes = ['rbspa','rbspb','tha','thd','the']
    apogee_filter_setting.mlt_around_noon = 5.
    apogee_filter_setting.mlt_range = [-1,1]*(12-apogee_filter_setting.mlt_around_noon)  ; [-12,12] hour.
    apogee_filter_setting.file = join_path([project.data_dir,'candidate_list.txt'])
    apogee_filter_setting.log_file = join_path([project.data_dir,'candidate_list_'+sgnum2str(apogee_filter_setting.mlt_around_noon)+'hr_around_noon.log'])

    if file_test(apogee_filter_setting.log_file) ne 0 then file_delete, apogee_filter_setting.log_file
    stouch, apogee_filter_setting.log_file
    nprobe = n_elements(apogee_filter_setting.probes)

    candidate_list = list()
    foreach time_range, storm_times do begin
        lprmsg, '', apogee_filter_setting.log_file
        lprmsg, 'Processing the storm from '+time_string(time_range[0])+' to '+time_string(time_range[1])+' ...', apogee_filter_setting.log_file

        date = mean(time_range) & date = date-(date mod project.constant.secofday)
        orbit_time_range = date+[-1,1]*project.constant.secofday

    ;---Get apogee MLT.
        apogee_mlt = list()
        foreach probe, apogee_filter_setting.probes do begin
            call_procedure, project[probe].routine_name+'_read_orbit', orbit_time_range, probe=project[probe].probe, errmsg=errmsg
            rvar = project[probe].prefix+'_r_gsm'
            find_apogee, rvar, apogee_times=times
            rgsm = get_var_data(rvar, at=times)
            rmag = cotran(rgsm, times, 'gsm2mag')
            if n_elements(times) eq 1 then rmag = transpose(rmag)   ; in rare cases, there is only 1 apogee over 2 days for themis.
            mlon = atan(rmag[*,1],rmag[*,0])*project.constant.deg
            mlt = mlon2mlt(mlon, times)
            apogee_mlt.add, mean(mlt)
        endforeach
        apogee_mlt = apogee_mlt.toarray()   ; should already in [-12,12].
        apogee_mlt = round(apogee_mlt*10)*0.1   ; round to the first decimal.
        foreach probe, apogee_filter_setting.probes, ii do lprmsg, string(strupcase(probe),format='(A10)')+string(apogee_mlt[ii],format='(F5.1)')+' MLT', apogee_filter_setting.log_file


    ;---Exclude events when probes are around noon.
        index = where_pro(apogee_mlt, apogee_filter_setting.mlt_range, count=count)
        if count ne nprobe then begin
            lprmsg, 'Only '+sgnum2str(count)+' probes in ['+strjoin(string(apogee_filter_setting.mlt_range,format='(I0)'),',')+'] MLT ...', apogee_filter_setting.log_file
            continue
        endif


    ;---Exclude events if not all probes are in the same side of the midnight.
        signs = apogee_mlt gt 0
        index = where(signs eq 1, count)
        if count ne 0 and count ne nprobe then begin
            lprmsg, 'Not all probes are on one side of the midnight ...', apogee_filter_setting.log_file
            continue
        endif

;    ;---Exclude events if the max MLT separation is too large.
;        index = sort(apogee_mlt)
;        apogee_mlt_diff = apogee_mlt[index]
;        apogee_mlt_diff = apogee_mlt_diff[1:-1]-apogee_mlt_diff[0:-2]
;        max_apogee_mlt_diff = max(apogee_mlt_diff)
;        if max_apogee_mlt_diff gt apogee_filter_setting.mlt_separation then begin
;            lprmsg, 'Max MLT separation: '+sgnum2str(max_apogee_mlt_diff,ndec=1)+' > '+sgnum2str(apogee_filter_setting.mlt_separation)+' hr ...'
;            continue
;        endif

        lprmsg, 'Found a candidate ...', apogee_filter_setting.log_file
        candidate_list.add, time_range
    endforeach


;---Save the result to project.
    candidate_id = list()
    foreach time_range, candidate_list do candidate_id.add, time_string(mean(time_range),tformat='YYYY_MMDD_hh')
    project.apogee_filter_setting = apogee_filter_setting
    project.candidate_list = candidate_list
    project.candidate_id = candidate_id.toarray()
    azim_prop_update_project, project


;---Output to a text file.
    openw, lun, apogee_filter_setting.file, /get_lun
    foreach time_range, candidate_list, ii do printf, lun, candidate_id[ii]+'    '+strjoin(time_string(time_range),'    ')
    free_lun, lun
end
