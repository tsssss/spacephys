
function micro_injection_load_round1_event_list_time_range, event_list_file

    if n_elements(event_list_file) eq 0 then begin
        project_id = 'micro_injection'
        project_info = project_load(project_id)
        event_list_file = join_path([project_info.data_dir,project_id+'_survey_round1_event_list.txt'])
    endif

    if file_test(event_list_file) eq 0 then stop

    lines = read_all_lines(event_list_file)
    nline = n_elements(lines)
    trs = dblarr(nline,2)
    time_len = 21
    time_format = 'YYYY_MMDDhh:mm'
    secofday = constant('secofday')
    for ii=0,nline-1 do begin
        tinfo = lines[ii]
        infos = strsplit(strmid(tinfo,0,time_len),' -', extract=1)
        if n_elements(infos) ne 3 then stop
        print, strjoin(infos, ' ')
        t0 = time_double(infos[0]+infos[1],tformat=time_format)
        t1 = time_double(infos[0]+infos[2],tformat=time_format)
        if t1 lt t0 then t1 += secofday
        trs[ii,*]= [t0,t1]
    endfor

    return, trs

end

trs = micro_injection_load_round1_event_list_time_range()
end