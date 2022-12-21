;+
; Generate a list for rbsp injections.
;-

function rbsp_themis_injection_gen_rbsp_injection_list, filename=injection_file
    

    if n_elements(injection_file) eq 0 then message, 'Require filename ...'
    if file_test(injection_file) eq 1 then return, injection_file
    
    txt_file = rbsp_themis_injection_identify_rbsp_injection()
    lines = read_all_lines(txt_file)

    secofday = 86400d
    coord = 'gse'
    tab = '    '
    sep = '  '
    headers = [$
        'Probe        Start and end times                         GSE X,Y,Z (Re)                      L (#)           MLT (h)           MLat (deg)', $
        '-    ----------------------------------    ------------------------------------------    ------------     -------------    ----------------' ]
    ftouch, injection_file
    foreach header, headers do lprmsg, header, injection_file
    nheader = n_elements(headers)
    format = 
    

    foreach line, lines do begin
        infos = strsplit(line,' ', extract=1)
        probe = strmid(infos[0],4,1)
        time_range = time_double(infos[1:2])
        r_var = rbsp_read_orbit(time_range, probe=probe, coord=coord)
        r_gses = get_var_data(r_var, at=time_range)
        
        lshell_var = rbsp_read_lshell(time_range, probe=probe)
        lshells = get_var_data(lshell_var, at=time_range)
        
        mlt_var = rbsp_read_mlt(time_range, probe=probe)
        mlts = get_var_data(mlt_var, at=time_range)
        
        mlat_var = rbsp_read_mlat(time_range, probe=probe)
        mlats = get_var_data(mlat_var, at=time_range)
        
        msg = ''
        msg += probe+tab
        msg += strjoin(time_string(time_range,tformat='YYYY-MM-DD/hh:mm'),sep)+tab
        msg += strjoin(string(r_gses[0,*],format='(F6.3)'),',')+sep
        msg += strjoin(string(r_gses[1,*],format='(F6.3)'),',')+tab
        msg += strjoin(string(lshells,format='(F5.3)'),sep)+tab
        msg += strjoin(string(mlts,format='(F6.3)'),sep)+tab
        msg += strjoin(string(mlats,format='(F7.3)'),sep)+tab
        lprmsg, msg, injection_file
    endforeach

    return, injection_file

end


injection_file = rbsp_themis_injection_gen_rbsp_injection_list()
end