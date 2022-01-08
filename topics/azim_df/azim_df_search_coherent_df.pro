;+
; Search for overlapping events, and select the best.
;-
function azim_df_search_coherent_df, search_setting, project=project, $
    reset=reset, test_time=test_time

    retval = list()
    tab = constant('4space')

;---Settings for selecting subgroup.
    the_key = 'search_coherent_df'
    if ~search_setting.haskey(the_key) then message, 'No settings for coherent DF ...'
    search_info = search_setting[the_key]


;---Check if file exists, to avoid search again.
    file_suffix = search_info.file_suffix
    out_file = join_path([project.data_dir,file_suffix])
    if keyword_set(reset) then begin
        lprmsg, 'Resetting coherent DF search ...'
        file_delete, out_file, /allow_nonexistent
        lprmsg, 'Clear memory ...'
        del_data, '*'
    endif
    if file_test(out_file) eq 1 then begin
        if file_lines(out_file) eq 0 then file_delete, out_file
    endif


    if file_test(out_file) eq 0 or keyword_set(test_time) then begin
        candidates = azim_df_search_df_group(search_setting, project=project)
        ftouch, out_file

        ; Find overlapping events.
        ncandidate = candidates.length
        time_ranges = dblarr(ncandidate,2)
        foreach candidate, candidates, ii do time_ranges[ii,*] = candidate.time_range
        index = sort(time_ranges[*,0])
        time_ranges = time_ranges[index,*]
        candidates = candidates[index]

        overlap_list = list()
        current_list = list()
        foreach the_group, candidates do begin
            if current_list.length ne 0 then begin
                pre_group = current_list[-1]
                pre_time_range = pre_group.time_range
                the_time_range = the_group.time_range

                if the_time_range[0] gt pre_time_range[1] then begin
                    if current_list.length gt 1 then overlap_list.add, current_list
                    current_list = list()
                endif
            endif
            current_list.add, the_group
        endforeach

        id_list = list()
        foreach cluster, overlap_list do begin
            foreach candidate, cluster do id_list.add, candidate.id
        endforeach
        id_list = id_list.toarray()
        isolated_flag = intarr(ncandidate)
        foreach candidate, candidates, ii do begin
            index = where(id_list eq candidate.id, count)
            if count eq 1 then continue
            isolated_flag[ii] = 1
        endforeach
        index = where(isolated_flag eq 1, count)
        isolated_list = (count eq 0)? list(): candidates[index]


        selected_list = list()
        foreach events, overlap_list do begin
            nevent = events.length
            event_ids = intarr(nevent)
            foreach event, events, ii do event_ids[ii] = event.id
            lprmsg, 'Resolve overlapping events: '+strjoin(string(event_ids,format='(I0)'),',')+' ...'
            

            ; Check coherency in vmag.
            vmag_scatters = fltarr(events.length)
            foreach event, events, ii do begin
                vmags = list()
                foreach triad, event.triad_list do vmags.add, triad.vmag_xcor
                vmag_scatters[ii] = stddev(vmags.toarray())
            endforeach
            min_vmag_scatter = min(vmag_scatters)
            index = where(vmag_scatters eq min_vmag_scatter, count)
            lprmsg, 'Found '+string(count,format='(I0)')+' events with stddev of vmag_xcor '+$
                string(min_vmag_scatter,format='(I0)')+' km/s ...'
            selected_events = events[index]
            if count eq 1 then begin
                selected_list.add, selected_events[0]
                continue
            endif
        endforeach
        
        
gen_plot = 1
        if keyword_set(gen_plot) then begin
            root_dir = ['2020_04_coherent_dfs_test',search_setting.name]
            foreach events, overlap_list, event_id do begin
                dirname = [root_dir,events[0].search_name+string(event_id,format='(I0)')]
                azim_df_subgroup_gen_diagnostic_plot, events, dirname=dirname, project=project
            endforeach
            dirname = [root_dir,'selected_events']
            azim_df_subgroup_gen_diagnostic_plot, selected_list, dirname=dirname, project=project 
            dirname = root_dir           
            azim_df_subgroup_gen_diagnostic_plot, isolated_list, dirname=dirname, project=project
        endif
        stop

    endif

;    if keyword_set(test_time) then stop
;    return, azim_df_subgroup_read_file(out_file)
end
