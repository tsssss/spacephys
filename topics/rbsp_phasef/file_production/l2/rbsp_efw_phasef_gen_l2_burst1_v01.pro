;+
; Generate L2 burst1 data (E and B fields) v01 cdfs.
;-

pro rbsp_efw_phasef_gen_l2_burst1_v01_per_day, date, probe=probe, filename=file, log_file=log_file

    on_error, 0
    errmsg = ''

    msg = 'Processing '+file+' ...'
    lprmsg, msg, log_file

;---Check input.
    if n_elements(file) eq 0 then begin
        errmsg = 'cdf file is not set ...'
        lprmsg, errmsg, log_file
        return
    endif

    if n_elements(probe) eq 0 then begin
        errmsg = 'No input probe ...'
        lprmsg, errmsg, log_file
        return
    endif
    if probe ne 'a' and probe ne 'b' then begin
        errmsg = 'Invalid probe: '+probe+' ...'
        lprmsg, errmsg, log_file
        return
    endif
    prefix = 'rbsp'+probe+'_'
    rbspx = 'rbsp'+probe

    secofday = constant('secofday')
    the_date = time_double(date)
    the_date = the_date-(the_date mod secofday)
    time_step = 15*60d
    the_times = make_bins(the_date+[0,secofday], time_step)
    ntime = n_elements(the_times)-1
    for time_id=0, ntime-1 do begin
        time_range = the_times[time_id]+[0,time_step]
        vars = prefix+'efw_'+['eb1','mscb1']+'_mgse'
        del_data, vars
        rbsp_efw_read_l1_burst_efield, time_range, probe=probe, keep_spin_axis=1
        rbsp_efw_read_l1_burst_bfield, time_range, probe=probe
        nodata = 0
        foreach var, vars do begin
            if check_if_update(var) then nodata = 1
        endforeach
        if nodata then continue
        

        ; Save data to file.
        
    endfor


end


date = '2013-06-07'
probe = 'a'
file = join_path([homedir(),'test_l2_burst1_v01.cdf'])
rbsp_efw_phasef_gen_l2_burst1_v01_per_day, date, probe=probe, filename=file
end
