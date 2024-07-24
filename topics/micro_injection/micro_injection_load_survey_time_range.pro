;+
; Load time range for survey plots.
;-
function micro_injection_load_survey_time_range, root_dir=root_dir

    if n_elements(root_dir) eq 0 then root_dir = srootdir()
    plot_dir = root_dir
    data_dir = root_dir

;---Settings.
    dis_range = [5d,13]
    base = 'micro_injection_load_survey_time_range_v01.cdf'
    data_file = join_path([data_dir,base])

    survey_time_range_var = 'survey_time_ranges'
    probe = '1'
    prefix = 'mms'+probe+'_'
    mission_probe = 'mms'+probe
    

    if file_test(data_file) eq 0 then begin
    ;---Load orbit data over the search time range to get the time for each orbit.
        search_trs = micro_injection_load_search_time_range(root_dir=root_dir)
        nsearch_tr = n_elements(search_trs[*,0])
        survey_trs = []
        for tid=0,nsearch_tr-1 do begin
            search_tr = reform(search_trs[tid,*])
            r_gsm_var = lets_read_this(func='mms_read_orbit',search_tr, probe=mission_probe)
            mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)
            dis_var = mlat_vars['dis']
            diss = get_var_data(dis_var, times=times)
            index = where_pro(diss, '[]', dis_range, count=count)
            if count eq 0 then continue
            my_indexs = time_to_range(index,time_step=1)
            my_trs = times[my_indexs]
            ; This is basically to remove the start and end b/c they are not full orbits.
            del_diss = diss[my_indexs[*,1]]-diss[my_indexs[*,0]]
            index = where(abs(del_diss) le 0.5, count)
            if count ne 0 then begin
                my_indexs = my_indexs[index,*]
                my_trs = my_trs[index,*]
            endif
            survey_trs = [survey_trs,my_trs]
        endfor
        cdf_save_var, survey_time_range_var, value=survey_trs, save_as_one=1, filename=data_file
    endif
    survey_trs = cdf_read_var(survey_time_range_var, filename=data_file)
    
    ; Round to full minutes.
    survey_trs = survey_trs-(survey_trs mod 60)
    survey_trs[*,1] += 60
    

    return, survey_trs

end


survey_trs = micro_injection_load_survey_time_range()
end