;+
; Test moving cross correlation. Specify the time range of the first spacecraft,
; cross correlate it to the next one, determine the time range of the second
; spacecraft, and so on.
;-

pro test_moving_cross_correlation, project, event_id=event_id, time_range=time_range

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events
    if n_elements(event_id) eq 0 then event_id = '2014_0828_10'
    if n_elements(time_range) eq 0 then time_range = time_double(['2014-08-28/10:13','2014-08-28/10:25'])

    data_file = join_path([project.data_dir,event_id+'_basic_data.tplot'])
    if file_test(data_file) eq 0 then azim_prop_gen_survey_plot, /save_data
    if file_test(data_file) eq 0 then message, 'Something is wrong when loading the data ...'
    events[event_id].file = data_file
    tplot_restore, file = data_file

    probes = events[event_id].probes
    event_time_range = events[event_id].time_range
    mean_time = mean(event_time_range)
    mean_mlts = list()
    foreach probe, probes do mean_mlts.add, get_var_data(probe+'_mlt', at=mean_time)
    mean_mlts = mean_mlts.toarray()
    sorted_probes = probes[sort(abs(mean_mlts))]    ; closest to the midnight comes first.

    ; Unify the time rate for the dB tilt angle.
    vars = sorted_probes+'_db_tilt'
    tilt_data_rate = 5.
    times = make_bins(event_time_range, tilt_data_rate)
    foreach var, vars do interp_time, var, times

    ; Running cross correlate.
    cc_info = hash()
    shift_limit = 40.*60
    shifts = smkarthm(-shift_limit*0.2,shift_limit*0.8,tilt_data_rate,'dx')
    nshift = n_elements(shifts)
    lags = round(shifts/tilt_data_rate)
    foreach probe, sorted_probes, ii do begin
        cc_info[probe] = dictionary()
        if ii eq 0 then begin
            cc_info[probe].time_shift = 0
        endif else begin
            get_data, sorted_probes[ii-1]+'_db_tilt_sector', uts, f0s
            x0s = round((uts-event_time_range[0])/tilt_data_rate)
            nx0 = n_elements(x0s)
            get_data, sorted_probes[ii]+'_db_tilt', uts, f1s
            x1s = round((uts-event_time_range[0])/tilt_data_rate)
            nx1 = n_elements(x1s)
            corr_2d = fltarr(nshift)
                ; Calculate the new time.
                for kk=0, nshift-1 do begin
                    index = x0s[0]+lags[kk]
                    if index lt 0 then continue
                    if index+nx0 ge nx1 then continue
                    tf1 = f1s[index:index+nx0-1]
                    corr_2d[kk] = c_correlate(tf1,f0s,0)
                endfor
            max_corr = max(corr_2d)
            index = where(corr_2d eq max_corr)
            the_shift = shifts[index[0]]

            get_data, sorted_probes[ii-1]+'_db_tilt_sector', uts, f0s
            ; Apply scale.
            uts = smkarthm(uts[0],tilt_data_rate,n_elements(uts),'x0')
            ; Apply shift.
            uts += the_shift
            cc_info[probe].time_shift = the_shift

            store_data, sorted_probes[ii-1]+'_db_tilt_sector_new', uts, f0s
            tplot, [sorted_probes[ii-1]+'_db_tilt_sector_new', sorted_probes[ii]+'_db_tilt'], trange=event_time_range
        endelse
        var = probe+'_db_tilt'
        get_data, var, uts, dat, limits=lim
        index = where_pro(times, time_range)
        store_data, var+'_sector', uts[index], dat[index], limits=lim
        cc_info[probe].var = probe+'_db_tilt_sector'
    endforeach

    time_lag = 0.
    foreach probe, sorted_probes, ii do begin
        time_lag -= cc_info[probe].time_shift
        get_data, probe+'_db_tilt', uts, dat, limits=lim
        store_data, probe+'_db_tilt_new', uts+time_lag, dat, limits=lim
    endforeach
    
    tplot, sorted_probes+'_db_tilt_new', trange=event_time_range
    stop

stop

end
