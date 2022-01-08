;+
; Find substorm for given time range, [min,max]_value, [min,max]_duration.
; Returns the time_ranges, in [N,2].
;
; time_range. The time range in ut sec.
; pad_time=. A number in sec to pad the main phase.
; min_value=. A number in nT for min value.
; min_duration=. A number in sec for the duration below min_value.
; max_value=. A number in nT for max value.
; max_duration=. A number in sec for the duration above max_value.
;-

function find_substorm, time_range, pad_time=pad_time, $
    min_value=min_value, min_duration=min_duration, $
    max_value=max_value, max_duration=max_duration


;---Check input.
    retval = !null
    if n_elements(time_range) ne 2 then begin
        errmsg = handle_error('No input time_range ...')
        return, retval
    endif

    default_duration = 1800.
    the_var = 'ae'
    if n_elements(pad_time) eq 0 then pad_time = 900.
    if n_elements(min_value) eq 0 then min_value = 500.
    if n_elements(min_duration) eq 0 then min_duration = default_duration
    if n_elements(max_value) eq 0 then max_value = 800.
    if n_elements(max_duration) eq 0 then max_duration = !null

;---Load data.
    time_step = 60. ; omni data in 1 min cadence.
    common_times = make_bins(time_range, time_step)
    if check_if_update(the_var, time_range) then omni_read_index, time_range
    the_data = get_var_data(the_var, at=common_times)

;---Min value. Find main phase.
    if n_elements(min_value) ne 0 then begin
        lprmsg, 'Find data above min_value: '+sgnum2str(min_value)+' ...'
        index = where(the_data ge min_value, count)
        if count eq 0 then begin
            errmsg = handle_error('No data found ...')
            return, retval
        endif
        times = common_times[index]
        time_ranges = time_to_range(times, time_step=time_step)
        if n_elements(min_duration) ne 0 then begin
            lprmsg, 'Find data longer than '+sgnum2str(min_duration)+' ...'
            durations = time_ranges[*,1]-time_ranges[*,0]
            index = where(durations ge min_duration, count)
            if count eq 0 then begin
                errmsg = handle_error('No data found ...')
                return, retval
            endif
            time_ranges = time_ranges[index,*]
        endif
    endif

;---Pad time. Expand main phase.
    time_ranges[*,0] -= pad_time
    time_ranges[*,1] += pad_time

;---Max value.
    if n_elements(max_value) ne 0 then begin
        ntime_range = n_elements(time_ranges)/2
        flags = bytarr(ntime_range)
        for ii=0, ntime_range-1 do begin
            the_time_range = reform(time_ranges[ii,*])
            tdata = the_data[lazy_where(common_times, '[]', the_time_range)]
            index = where(tdata ge max_value, count)
            flags[ii] = count ne 0
        endfor
        index = where(flags eq 1, count)
        if count eq 0 then begin
            errmsg = handle_error('No time_range above max_value: '+sgnum2str(max_value)+' ...')
            return, retval
        endif
        time_ranges = time_ranges[index,*]
        if n_elements(max_duration) ne 0 then begin
            durations = time_ranges[*,1]-time_ranges[*,0]
            index = where(durations ge max_duration, count)
            if count eq 0 then begin
                errmsg = handle_error('No data above max_value longer than '+sgnum2str(max_duration)+' ...')
                return, retval
            endif
            time_ranges = time_ranges[index,*]
        endif
    endif


;---Merge overlapping events.
    lprmsg, 'Combine overlapping subtorms ...'
    ntime_range = n_elements(time_ranges)/2
    flags = intarr(ntime_range)
    for ii=1, ntime_range-1 do begin
        flags[ii] = flags[ii-1]
        if max(time_ranges[ii-1,*]) lt min(time_ranges[ii,*]) then flags[ii] += 1
    endfor
    nsubstorm_time = max(flags)+1
    substorm_times = dblarr(nsubstorm_time,2)
    for ii=0, nsubstorm_time-1 do begin
        substorm_times[ii,*] = minmax(time_ranges[where(flags eq ii),*])
    endfor

    return, substorm_times

end

time_range = time_double(['2007-03-01','2009-09-01'])
time_range = time_double(['2008-02-01','2008-03-01'])
time_range = time_double(['2012-10-01','2017-10-01'])
time_ranges = find_substorm(time_range)
ntime_range = n_elements(time_ranges)/2
for ii=0, ntime_range-1 do print, string(ii+1,format='(I0)')+'    '+strjoin(time_string(reform(time_ranges[ii,*]), tformat='YYYY-MM-DD/hh:mm'), ' ')
end
