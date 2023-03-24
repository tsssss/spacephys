
function rbsp_themis_injection_identify_rbsp_injection, filename=txt_file

    pinfo = rbsp_themis_injection_load_project()
    search_time_range = pinfo['search_time_range']
    version = pinfo['latest_version']

    probes = ['a','b']
    my_dir = join_path([pinfo['root_dir'],'identify_rbsp_injection_'+version])
    ;my_dir = !null
    secofday = constant('secofday')
    dates = make_bins(search_time_range+[0,-1]*secofday, secofday)
    if n_elements(txt_file) eq 0 then begin
        txt_file = join_path([my_dir,'identify_rbsp_injection_'+version+'.txt'])
    endif
    if file_test(txt_file) eq 1 then return, txt_file else ftouch, txt_file

    log_file = txt_file+'.log'
    if file_test(log_file) eq 1 then begin
        file_delete, log_file
        ftouch, log_file
    endif

    msg = 'YYYY_MMDD   RBSPA   RBSPB'
    lprmsg, msg, log_file
    foreach date, dates do begin
        msg = time_string(date,tformat='YYYY_MMDD')
        foreach probe, probes do begin
            time_range = date+[0,secofday]
            trs = rbsp_identify_injection($
                time_range, mission_probe='rbsp'+probe)
            ntr = n_elements(trs)*0.5
            msg = msg+'    '+string(ntr,format='(I3)')
            for ii=0,ntr-1 do begin
                txt = 'rbsp'+probe+'    '+$
                    strjoin(time_string(reform(trs[ii,*]),tformat='YYYY-MM-DD/hh:mm'),'  ')
                lprmsg, txt, txt_file
            endfor
        endforeach
        lprmsg, msg, log_file
    endforeach
    
    return, txt_file

end


tmp = rbsp_themis_injection_identify_rbsp_injection()
end