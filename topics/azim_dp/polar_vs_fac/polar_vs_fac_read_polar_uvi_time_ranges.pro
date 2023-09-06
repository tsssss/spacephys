;+
; Check when SECS and Polar are both available, what are the time ranges Polar have data.
;-

pro polar_vs_fac_read_polar_uvi_time_ranges_gen_file, project_time_range, file=polar_uvi_avail_file

    ; The flag.
    uvi_avail_var = 'po_uvi_avail'
    if check_if_update(uvi_avail_var) then begin
        ; The common times.
        time_step = 60.
        polar_uvi_times = make_bins(project_time_range, time_step)
        ntime = n_elements(polar_uvi_times)
        polar_uvi_avail = fltarr(ntime)

        uvi_var = 'FRAMERATE'
        secofday = constant('secofday')
        days = make_bins(project_time_range+[0,-1]*secofday, secofday, inner=1)
        foreach day, days do begin
            day_time_range = day+[0,secofday]
            lprmsg, strjoin(time_string(day_time_range), ' to '), log_file

            ; Adopted from polar_read_mlt_image.
            del_data, uvi_var
            polar_read_uvi, day_time_range, id='l1'
            get_data, uvi_var, ut1s, rate
            times = ut1s-(rate+4)*9.2

            ; No data.
            index = where_pro(times, '[]', day_time_range, count=count)
            if count eq 0 then continue
            times = times[index]

            ; Merge times to the common times.
            time_index = (times-polar_uvi_times[0])/time_step
            time_index = sort_uniq(floor(time_index))
            polar_uvi_avail[time_index] += 1

            store_data, uvi_avail_var, polar_uvi_times, polar_uvi_avail
        endforeach

        store_data, uvi_avail_var, polar_uvi_times, polar_uvi_avail
    endif
    get_data, uvi_avail_var, polar_uvi_times, polar_uvi_avail

    ; Write to file.
    ftouch, polar_uvi_avail_file
    tab = '    '
    msg = 'Polar UVI start and end times, searched in '+$
        strjoin(time_string(project_time_range), ' to ')
    lprmsg, msg, polar_uvi_avail_file
    index = where(polar_uvi_avail ge 1, count)
    if count eq 0 then begin
        polar_uvi_time_ranges = 0
    endif else begin
        polar_uvi_time_ranges = polar_uvi_times[time_to_range(index, time_step=1)]
    endelse
    ntime_range = n_elements(polar_uvi_time_ranges)*0.5
    for ii=0, ntime_range-1 do begin
        the_time_range = reform(polar_uvi_time_ranges[ii,*])
        msg = strjoin(time_string(the_time_range,tformat='YYYY-MM-DD/hh:mm'),tab)
        lprmsg, msg, polar_uvi_avail_file
    endfor

end

function polar_vs_fac_read_polar_uvi_time_ranges, min_duration=min_duration

    project_time_range = time_double(['2007-01-01','2008-04-22'])
    ;project_time_range = time_double(['2007-01-01','2007-01-05'])

    ; The overall Polar data availability.
    data_dir = join_path([srootdir(),'data'])
    polar_uvi_avail_file = join_path([data_dir,'polar_uvi_avail.txt'])
    if file_test(polar_uvi_avail_file) eq 0 then begin
        polar_vs_fac_read_polar_uvi_time_ranges_gen_file, project_time_range, file=polar_uvi_avail_file
    endif

    lines = read_all_lines(polar_uvi_avail_file)
    nline = n_elements(lines)
    if nline le 1 then return, !null

    lines = lines[1:*]
    ntime_range = n_elements(lines)
    polar_uvi_time_ranges = strarr(ntime_range,2)
    for ii=0,ntime_range-1 do polar_uvi_time_ranges[ii,*] = strsplit(lines[ii],' ',/extract)
    polar_uvi_time_ranges = time_double(polar_uvi_time_ranges)
    durations = polar_uvi_time_ranges[*,1]-polar_uvi_time_ranges[*,0]
    if n_elements(min_duration) eq 0 then min_duration = 30d*60
    index = where(durations ge min_duration, count)
    if count eq 0 then return, !null
    polar_uvi_time_ranges = polar_uvi_time_ranges[index,*]

    return, polar_uvi_time_ranges

end


polar_uvi_time_ranges = polar_vs_fac_read_polar_uvi_time_ranges()
end