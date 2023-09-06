;+
; Check when FACs have data, based on when Polar had MLT images.
;-

pro polar_vs_fac_read_fac_time_ranges_gen_file, fac_avail_file

    polar_mltimg_time_ranges = polar_vs_fac_read_polar_mltimg_time_ranges()
    nsection = n_elements(polar_mltimg_time_ranges)*0.5

    if file_test(fac_avail_file) eq 1 then file_delete, fac_avail_file
    tab = '    '
    ftouch, fac_avail_file
    msg = 'FAC start and end times'
    lprmsg, msg, fac_avail_file

    mlt_range = [-1,1]*9
    fac_var = 'thg_j_up_ewo'
    for section_id=0,nsection-1 do begin
        the_time_range = reform(polar_mltimg_time_ranges[section_id,*])
if the_time_range[0] lt time_double('2007-12-12') then continue

        del_data, fac_var
        themis_read_upward_current_ewo, the_time_range, mlt_range=mlt_range, errmsg=errmsg
        if errmsg ne '' then continue

        get_data, fac_var, times, ewo
        index = where_pro(times, '[]', the_time_range, count=count)
        if count eq 0 then continue

        time_range = minmax(times[index])
        duration = total(time_range*[-1,1])
        the_duration = total(the_time_range*[-1,1])
        if duration lt the_duration*0.8 then continue

        ; Pass.
        msg = strjoin(time_string(the_time_range,tformat='YYYY-MM-DD/hh:mm'),tab)
        lprmsg, msg, fac_avail_file
    endfor


end



function polar_vs_fac_read_fac_time_ranges

    ; The overall Polar data availability.
    data_dir = join_path([srootdir(),'data'])
    fac_avail_file = join_path([data_dir,'fac_avail.txt'])
    if file_test(fac_avail_file) eq 0 then begin
        polar_vs_fac_read_fac_time_ranges_gen_file, fac_avail_file
    endif

    lines = read_all_lines(fac_avail_file)
    nline = n_elements(lines)
    if nline le 1 then return, !null

    lines = lines[1:*]
    ntime_range = n_elements(lines)
    fac_time_ranges = strarr(ntime_range,2)
    for ii=0,ntime_range-1 do fac_time_ranges[ii,*] = strsplit(lines[ii],' ',/extract)
    fac_time_ranges = time_double(fac_time_ranges)
    return, fac_time_ranges

end

fac_time_ranges = polar_vs_fac_read_fac_time_ranges()
end
