;+
; Load Dst and find storms for given dst_threshold and time_range.
; Adopted from azim_prop_find_storm.
;-

pro azim_df_find_storm_themis, project=project


;---Prepare settings.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    settings = project.storm_setting
    the_key = 'themis_search_time_range'
    if ~project.haskey(the_key) then project[the_key] = time_double(['2007-03-01','2009-09-01'])
    time_range = project[the_key]
    dst_threshold = settings.dst_threshold
    pad_time = settings.pad_time
    main_phase_min_duration = settings.main_phase_min_duration
    file = join_path([project.data_dir,'storm_list_themis.txt'])


;---Load the Dst data.
    lprmsg, 'Loading Dst from '+time_string(time_range[0])+' to '+time_string(time_range[1])+' ...'
    dst_var = 'dst'
    dst_data_rate = 60.     ; 1 min.
    omni_read_index, time_range
    uniform_time, dst_var, dst_data_rate

;---Loop through each candidate.
    lprmsg, 'Loop through candidates ...'
    storms = list()
    get_data, dst_var, times, dst
    index = where(dst le dst_threshold, count)
    for ii=0, count-1 do begin
        time_start = times[index[ii]]
        for jj=ii+1, count-1 do begin
            if round((times[index[jj]]-times[index[jj-1]])/dst_data_rate) eq 1 then continue
            break
        endfor
        time_end = times[index[jj-1]]
        ii = jj
        if time_end-time_start lt main_phase_min_duration then continue
        lprmsg, 'Found a main phase: '+strjoin(time_string([time_start,time_end]),' to ')+' ...'
        storms.add, [time_start,time_end]
    endfor


    ; Extend the main phase to full storm.
    lprmsg, 'Expand main phase by +/-'+string(pad_time,format='(I0)')+' sec ...'
    foreach time_range, storms, ii do storms[ii] += [-1,1]*pad_time


    ; Throw away storms when Dst never above 0 nT.
    lprmsg, 'Throw away storms when Dst remained <0 ...'
    nstorm = n_elements(storms)
    flags = bytarr(nstorm)
    for ii=0, nstorm-1 do begin
        time_range = storms[ii]
        dst = get_var_data(dst_var, in=time_range)
        index = where(dst gt 0, count)
        flags[ii] = count gt 0
    endfor
    index = where(flags eq 1, nstorm)
    if nstorm eq 0 then begin
        lprmsg, 'No storm found ...'
        return
    endif
    storms = storms[index]


    ; Combine storms.
    lprmsg, 'Combine overlapping storms ...'
    nstorm = n_elements(storms)
    flags = intarr(nstorm)
    for ii=1, nstorm-1 do begin
        flags[ii] = flags[ii-1]
        if max(storms[ii-1]) lt min(storms[ii]) then flags[ii] += 1
    endfor
    nstorm_time = max(flags)+1
    storm_times = list(length=nstorm_time)
    storm_id = strarr(nstorm_time)
    for ii=0, nstorm_time-1 do begin
        storm_times[ii] = minmax(storms[where(flags eq ii)].toarray())
        storm_id[ii] = time_string(mean(storm_times[ii]),tformat='YYYY_MMDD_hh')
    endfor





;---Save the result to project.
    lprmsg, 'Saving results to project.storm_times ...'
    project.done_find_storm_themis = 1
    project.storm_times_themis = storm_times
    project.storm_list_themis = storm_list
    project.storm_id_themis = storm_id
    update_project, project


;---Output to a text file.
    lprmsg, 'Saving results to '+file+' ...'
    openw, lun, file, /get_lun
    foreach time_range, storm_times do printf, lun, time_string(time_range)
    free_lun, lun

end
