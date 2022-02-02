;+
; Check when SECS and Polar are both available, what are the time ranges Polar have good data.
;-

function polar_mltimg_check_quality, mltimgs

    ntime = n_elements(mltimgs[*,0,0])
    flags = fltarr(ntime)

    min_count = 5e3
    min_stddev = 20

    for time_id=0,ntime-1 do begin
        timg = reform(mltimgs[time_id,*,*])

        ; We want enough number of pixels.
        index = where(timg gt 0, count)
        if count le min_count then continue

        ; To remove image with no physical signal.
        stddev = stddev(timg[index])
        if stddev le min_stddev then continue

        flags[time_id] = 1
    endfor

    return, flags

end

pro polar_vs_fac_read_polar_mltimg_time_ranges_gen_file, polar_mltimg_avail_file, min_duration=min_duration

    polar_uvi_time_ranges = polar_vs_fac_read_polar_uvi_time_ranges(min_duration=min_duration)
    ntime_range = n_elements(polar_uvi_time_ranges)*0.5


    if file_test(polar_mltimg_avail_file) eq 1 then file_delete, polar_mltimg_avail_file
    tab = '    '
    ftouch, polar_mltimg_avail_file
    msg = 'Polar MLT image start and end times'
    lprmsg, msg, polar_mltimg_avail_file

    for section_id=0,ntime_range-1 do begin
        the_time_range = reform(polar_uvi_time_ranges[section_id,*])
        polar_read_mltimg, the_time_range

        ; Raw flags.
        get_data, 'po_mltimg', times, mltimgs
        raw_flags = polar_mltimg_check_quality(mltimgs)

        ; Remove small chunks and small gaps.
        time_step = 60.
        common_times = make_bins(the_time_range, time_step)
        ntime = n_elements(common_times)
        flags = interp(raw_flags, times, common_times) eq 1
        step = 6
        flags = interp(flags[0:*:step], common_times[0:*:step], common_times) ne 0

        ; Select big chunks of good data.
        index = where(flags eq 1, count)
        if count eq 0 then continue

        time_ranges = common_times[time_to_range(index, time_step=1)]
        durations = time_ranges[*,1]-time_ranges[*,0]

        index = where(durations ge min_duration, count)
        if count eq 0 then continue

        time_ranges = time_ranges[index,*]
        for ii=0,count-1 do begin
            the_time_range = reform(time_ranges[ii,*])
            msg = strjoin(time_string(the_time_range,tformat='YYYY-MM-DD/hh:mm'),tab)
            lprmsg, msg, polar_mltimg_avail_file
        endfor
    endfor

end


function polar_vs_fac_read_polar_mltimg_time_ranges, min_duration=min_duration

    if keyword_set(min_duration) eq 0 then min_duration = 30d*60

    ; The overall Polar data availability.
    data_dir = join_path([srootdir(),'data'])
    polar_mltimg_avail_file = join_path([data_dir,'polar_mltimg_avail.txt'])
    if file_test(polar_mltimg_avail_file) eq 0 then begin
        polar_vs_fac_read_polar_mltimg_time_ranges_gen_file, polar_mltimg_avail_file, min_duration=min_duration
    endif

    lines = read_all_lines(polar_mltimg_avail_file)
    nline = n_elements(lines)
    if nline le 1 then return, !null

    lines = lines[1:*]
    ntime_range = n_elements(lines)
    polar_mltimg_time_ranges = strarr(ntime_range,2)
    for ii=0,ntime_range-1 do polar_mltimg_time_ranges[ii,*] = strsplit(lines[ii],' ',/extract)
    polar_mltimg_time_ranges = time_double(polar_mltimg_time_ranges)

    ; Combine time ranges that are not far away.
    max_sep = 5*60d
    nsection = n_elements(polar_mltimg_time_ranges)*0.5
    seps = polar_mltimg_time_ranges[1:nsection-1,0]-polar_mltimg_time_ranges[0:nsection-2,1]
    tr_list = list()
    for section_id=1,nsection-1 do begin
        if n_elements(tr_list) eq 0 then curr_tr = polar_mltimg_time_ranges[0,*] else curr_tr = tr_list[-1]
        next_tr = polar_mltimg_time_ranges[section_id,*]
        if seps[section_id-1] gt max_sep then begin
            tr_list.add, next_tr
        endif else begin
            curr_tr[1] = next_tr[1]
            tr_list[-1] = curr_tr
        endelse
    endfor
    polar_mltimg_time_ranges = reform(tr_list.toarray())
    return, polar_mltimg_time_ranges

end


polar_mltimg_time_ranges = polar_vs_fac_read_polar_mltimg_time_ranges()
end
