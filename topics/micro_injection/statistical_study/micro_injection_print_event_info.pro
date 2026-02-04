
function micro_injection_print_event_info, tr

    r_var = lets_read('orbit', tr, source=['mms','4'], coord='sm')
    mlat_vars = lets_read_mlat_vars(orbit_var=r_var)
    mlt_var = mlat_vars['mlt']
    
    dis = mean(snorm(get_var_data(r_var, in=tr)))
    mlt = mean(get_var_data(mlt_var, in=tr))
    if mlt le 0 then mlt += 24
    msg = strjoin(time_string(tr,tformat='YYYY_MMDD_hhmm'),'  ')
    msg += string(dis,format='(F5.1)')+' Re'+'    '
    msg += string(mlt,format='(F4.1)')+' MLT'
    return, msg

end

google_dir = googledir()
data_dir = join_path([google_dir,'works','2024_microinjection','data'])
plot_dir = join_path([google_dir,'works','2024_microinjection','plot'])
base = 'micro_injection_survey_round1_event_list.txt'
manual_tr_file = join_path([data_dir,base])
lines = read_all_lines(manual_tr_file)
nline = n_elements(lines)
secofday = constant('secofday')

trs = dblarr(nline,2)
for ii=0,nline-1 do begin
    info = lines[ii]
    date = strmid(info,0,9)
    times = strmid(info,10,11)
    tr_str = date+'/'+strsplit(times,'-',extract=1)
    the_tr = time_double(tr_str,tformat='YYYY_MMDD/hh:mm')
    if the_tr[1] le the_tr[0] then begin
        the_tr[1] += secofday
    endif
    trs[ii,*] = the_tr
endfor
nevent = n_elements(trs[*,0])

event_info_file = join_path([data_dir,'event_info_manual_time_range_v01.txt'])
if file_test(event_info_file) eq 0 then ftouch, event_info_file
msgs = strarr(nevent)
for ii=0,nevent-1 do begin
    tr = reform(trs[ii,*])
    lprmsg, micro_injection_print_event_info(tr), event_info_file
endfor



end