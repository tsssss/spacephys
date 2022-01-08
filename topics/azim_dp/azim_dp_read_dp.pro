;+
; Read dipolarization for the given time_range and mlt_range.
;-

function azim_dp_read_dp, time_range, probe=probe, mlt_range=mlt_range

    retval = list()


    ramp_list = azim_dp_read_ramp(time_range, probe=probe, mlt_range=mlt_range)
    nramp = ramp_list.length



;---Select isolated ramps.
    if nramp ge 2 then begin
        clear_window_size = 10.*60  ; sec.
        obs_times = dblarr(nramp)
        foreach ramp, ramp_list, ramp_id do obs_times[ramp_id] = ramp.time
        is_isolated_ramp = bytarr(nramp)+1
        dtimes = obs_times[1:nramp-1]-obs_times[0:nramp-2]
        index = where(dtimes lt clear_window_size, count)
        if count ne 0 then is_isolated_ramp[index] = 0
        lprmsg, '', log_file
        lprmsg, 'Found '+string(nramp,format='(I0)')+' ramps that are not isolated ...', log_file
        index = where(is_isolated_ramp eq 1, nramp)
        if nramp eq 0 then return, retval
        ramp_list = ramp_list[index]
    endif
    

;---Select ramps with proper shape.
    df_index = list()
    foreach ramp, ramp_list, ramp_id do begin
        if azim_dp_is_good_ramp(ramp) eq 0 then continue
        df_index.add, ramp_id
    endforeach
    nramp = df_index.length
    if nramp eq 0 then return, retval
    df_index = df_index.toarray()
    ramp_list = ramp_list[df_index]

    return, ramp_list

end


; Test.
mlt_range = [0,9]
time_range = time_double(['2014-08-28/09:30','2014-08-28/11:00'])
probes = ['tha','thd','the','rbspa','rbspb','g13','g14','g15']


; Test.
mlt_range = [-6,0]
time_range = time_double(['2019-08-05','2019-08-06'])
probes = ['tha','thd','the','rbspa','rbspb','g14','g15','g16','g17','mms1']


dp_list = list()
foreach probe, probes do begin
    dp_list.add, azim_dp_read_dp(time_range, probe=probe, mlt_range=mlt_range), /extract
endforeach

; Sort by ramp_time.
ndp = dp_list.length
times = dblarr(ndp)
foreach dp, dp_list, dp_id do times[dp_id] = dp.time
index = sort(times)
dp_list = dp_list[index]

foreach dp, dp_list do print, dp.probe+'  '+time_string(dp.time)

end
