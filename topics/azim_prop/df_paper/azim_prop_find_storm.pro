;+
; Load Dst and find storms for given dst_threshold and time_range.
;-

pro azim_prop_find_storm, project, time_range=time_range, dst_threshold=dst_threshold, pad_time=pad_time, main_phase_min_duration=main_phase_min_duration


;---Settings.
    if n_elements(project) eq 0 then project = azim_prop_load_project()
    if n_elements(time_range) ne 2 then time_range = time_double(['2012-10-01','2017-10-01'])
    if n_elements(dst_threshold) eq 0 then dst_threshold = -50.
    if n_elements(pad_time) eq 0 then pad_time = project.constant.secofday*0.5
    if n_elements(main_phase_min_duration) eq 0 then main_phase_min_duration = 3600.

    storm_setting = dictionary()
    storm_setting.time_range = time_range
    storm_setting.dst_threshold = dst_threshold
    storm_setting.pad_time = pad_time
    storm_setting.main_phase_min_duration = main_phase_min_duration
    storm_setting.file = join_path([project.data_dir,'storm_list.txt'])


;---Load the Dst data.
    dst_var = 'dst'
    dst_data_rate = 60.     ; 1 min.
    omni = sread_omni(storm_setting.time_range)
    uts = sfmepoch(omni.epoch,'unix')
    times = make_bins(minmax(uts), dst_data_rate)
    dst_data = interpol(omni.sym_h, uts, times)
    store_data, dst_var, times, dst_data


;---Loop through each candidate.
    storms = list()
    get_data, dst_var, times, dst
    index = where(dst le storm_setting.dst_threshold, count)
    for ii=0, count-1 do begin
        time_start = times[index[ii]]
        for jj=ii+1, count-1 do begin
            if round((times[index[jj]]-times[index[jj-1]])/dst_data_rate) eq 1 then continue
            break
        endfor
        time_end = times[index[jj-1]]
        ii = jj
        if time_end-time_start lt main_phase_min_duration then continue
        storms.add, [time_start,time_end]
    endfor

    ; Extend the main phase to full storm.
    foreach time_range, storms, ii do storms[ii] += [-1,1]*pad_time

    ; Combine storms.
    nstorm = n_elements(storms)
    flags = intarr(nstorm)
    for ii=1, nstorm-1 do begin
        flags[ii] = flags[ii-1]
        if max(storms[ii-1]) lt min(storms[ii]) then flags[ii] += 1
    endfor

    nstorm_time = max(flags)+1
    storm_times = list(length=nstorm_time)
    for ii=0, nstorm_time-1 do storm_times[ii] = minmax(storms[where(flags eq ii)].toarray())
    

;---Save the result to project.
    project.storm_setting = storm_setting
    project.storm_times = storm_times
    azim_prop_update_project, project


;---Output to a text file.
    openw, lun, storm_setting.file, /get_lun
    foreach time_range, storm_times do printf, lun, time_string(time_range)
    free_lun, lun

end


azim_prop_find_storm
end
