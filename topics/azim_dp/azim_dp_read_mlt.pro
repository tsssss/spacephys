;+
; Read MLT data.
;
; time.
; probe=.
;-

pro azim_dp_read_mlt, time, probe=probe, errmsg=errmsg

    prefix = probe+'_'
    mlt_var = prefix+'mlt'
    if ~check_if_update(mlt_var, time) then return
    
    azim_dp_read_level2_data, time, probe=probe, datatype='mlt', errmsg=errmsg
    if errmsg ne '' then return

    get_data, mlt_var, times, mlts
    index = where(finite(mlts[*,0]), count)
    if count eq 0 then begin
        errmsg = 'No data ...'
        del_data, mlt_var
        return
    endif

end

time_range = time_double(['2019-08-01','2019-09-10'])
probes = ['rbspa','rbspb','mms1','g13','g14','g15','tha','thd','the']
foreach probe, probes do begin
    azim_dp_read_mlt, time_range, probe=probe
endforeach
end