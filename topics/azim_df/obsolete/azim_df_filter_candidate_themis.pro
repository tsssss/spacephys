;+
; Load candidate list, and AE index.
; Pick out storms when enough # of sc in the desired region.
;-

pro azim_df_filter_candidate_themis, project=project


;---Prepare settings.
    if n_elements(project) eq 0 then project = azim_df_load_project()

    time_range = project['themis_search_time_range']
    time_step = project['themis_search_time_step']
    common_times = make_bins(time_range, time_step)

    candidate_list = join_path([project.data_dir,'azim_df_search_themis_good_triangle_list.txt'])
    lines = read_all_lines(candidate_list)
    candidates = list()
    foreach line, lines do begin
        info = strsplit(line,' ',/extract)
        candidates.add, time_double(info[[0,2]])
    endforeach

    if check_if_update('ae', time_range) then omni_read_index, time_range
    all_ae = get_var_data('ae', at=common_times)

    ae_flags = list()
    min_ae = 500.   ; nT.
    min_ae_duration = 30*60.    ; sec.
    events = list()
    foreach the_time_range, candidates do begin
        the_flag = 0
;        index = where_pro(common_times, '[]', the_time_range+[-1,1]*min_ae_duration, count=count)
        index = where_pro(common_times, '[]', the_time_range, count=count)
        if count eq 0 then begin
            ae_flags.add, the_flag
            continue
        endif
        the_ae = all_ae[index]
        the_times = common_times[index]
        index = where(the_ae ge min_ae, count)
        if count le 1 then begin
            ae_flags.add, the_flag
            continue
        endif

; test_time = time_double('2008-03-08/13:05')
; if test_time gt the_time_range[0] and test_time lt the_time_range[1] then stop

        time_ranges = time_to_range(the_times[index], time_step=time_step)
        ntime_range = n_elements(time_ranges)
        durations = time_ranges[*,1]-time_ranges[*,0]
        index = where(durations ge min_ae_duration, count)
        if count eq 0 then begin
            ae_flags.add, the_flag
            continue
        endif

        the_flag = 1
        ae_flags.add, the_flag
        for ii=0, count-1 do events.add, reform(time_ranges[index[ii],*])
    endforeach

    ae_flags = ae_flags.toarray()
    index = where(ae_flags eq 1, count)
    if count eq 0 then message, 'No candidate found, stop here ...'
    time_ranges = events.toarray()   ; in [N,2]
    ntime_range = n_elements(time_ranges)/2

;---Write result to a file.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    log_file = join_path([project.data_dir,'azim_df_candidate_themis.txt'])
    openw, lun, /get_lun, log_file
    for ii=0, ntime_range-1 do begin
        printf, lun, strjoin(reform(time_string(time_ranges[ii,*])),' to ')
    endfor
    free_lun, lun

end

azim_df_filter_candidate_themis
end
