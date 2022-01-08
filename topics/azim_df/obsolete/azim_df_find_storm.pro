;+
; Load Dst and find storms for given dst_threshold and time_range.
; Adopted from azim_prop_find_storm.
;-

pro azim_df_find_storm, time_range, project=project


;---Prepare settings.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    the_key = 'storm_setting'
    if ~project.haskey(the_key) then message, 'Need to init '+the_key+' ...'
    settings = project[the_key]
    if n_elements(time_range) ne 2 then time_range = settings.time_range
    dst_threshold = settings.dst_threshold
    pad_time = settings.pad_time
    main_phase_min_duration = settings.main_phase_min_duration
    file = join_path([project.data_dir,settings.file_suffix])


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
    project.done_find_storm = 1
    project.storm_times = storm_times
    update_project, project
    project.storm_list = storm_list
    project.storm_id = storm_id


;---Output to a text file.
    lprmsg, 'Saving results to '+file+' ...'
    openw, lun, file, /get_lun
    foreach time_range, storm_times do printf, lun, time_string(time_range)
    free_lun, lun

end

time_range = time_double(['2012-10-01','2017-10-01'])
azim_df_find_storm, time_range
end
