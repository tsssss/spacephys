
function micro_injection_read_mi_times, mission_probe

    search_trs = micro_injection_load_search_time_range()
    nsearch_tr = n_elements(search_trs[*,0])
    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = mission+probe+'_'
    project_info = micro_injection_stat_load_project()
    common_time_step = project_info['common_time_step']
    data_dir = project_info['data_dir']
    file = join_path([data_dir,'micro_injection_'+prefix+'mi_times.cdf'])

    mi_time_var = prefix+'mi_times'
    if check_if_update(mi_time_var) then begin
        if file_test(file) eq 0 then begin
            mi_times = list()
            for ii=0,nsearch_tr-1 do begin
                tr = reform(search_trs[ii,*])
                tr_id = string(ii+1,format='(I0)')
                tr_suffix = '_'+tr_id

                ; Load MI times.
                mi_trs = micro_injection_stat_read_event_times(tr, mission_probe=mission_probe)
                nmi_tr = n_elements(mi_trs[*,0])
                for tid=0,nmi_tr-1 do begin
                    mi_tr = reform(mi_trs[tid,*])
                    mi_times.add, make_bins(mi_tr,common_time_step), extract=1
                endfor
            endfor
            mi_times = mi_times.toarray()
            cdf_save_var, mi_time_var, value=mi_times, filename=file
        endif
        mi_times = cdf_read_var(mi_time_var, filename=file)
        mi_time_var = var_store(mi_time_var, mi_times)
    endif

    mi_times = var_get_data(mi_time_var)
    return, mi_times

end


mission_probe = ['mms','1']
mi_times = micro_injection_read_mi_times(mission_probe)
end