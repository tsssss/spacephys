;+
; Search events for possible azim_vs_fac events in ROI.
;-

function azim_dp_vs_fac_search_roi, time_range, probes=probes, $
    reset=reset, $
    mlt_range=mlt_range, rxy_range=rxy_range, pdyn=pdyn, roi_min_count=roi_min_count


;---Check input.
    roi_list = list()
    if n_elements(time_range) ne 2 then return, roi_list
    if n_elements(roi_min_count) eq 0 then roi_min_count = 5    ; #.
    if n_elements(probes) lt roi_min_count then probes = ['th'+letters('e'),'rbsp'+letters('b'),'mms1','g'+string(make_bins([10,17],1),format='(I0)')]
    if n_elements(mlt_range) ne 2 then mlt_range = [-1,1]*9     ; hr.
    if n_elements(rxy_range) ne 2 then rxy_range = [4.,30]      ; Re.
    if n_elements(pdyn) eq 0 then pdyn = 10.    ; nPa.
    if n_elements(roi_min_duration) eq 0 then roi_min_duration = 1800.  ; sec.
    errmsg = ''


;---Settings.
    root_dir = join_path([googledir(),'works','azim_dp_vs_fac','data'])
    if file_test(root_dir) eq 0 then file_mkdir, root_dir
    base = 'azim_dp_vs_fac_search_roi_'+time_string(mean(time_range),tformat='YYYY')
    log_file = join_path([root_dir,base+'.log'])
    out_file = join_path([root_dir,base+'.txt'])
    if keyword_set(test) then begin
        log_file = -1
        out_file = -1
    endif

    files = [log_file,out_file]
    if keyword_set(reset) then begin
        foreach file, files do if file_test(file) eq 1 then file_delete, file
    endif


;---Try to read roi_list.
    if file_test(out_file) eq 1 then begin
        roi_list = azim_dp_roi_read(filename=out_file)
        return, roi_list
    endif


    foreach file, files do if file_test(file) eq 0 then ftouch, file

;---Search for ROI times.
    tab = constant('4space')

    lprmsg, 'Search roi_list for azim_dp_vs_fac', log_file
    lprmsg, '', log_file
    lprmsg, 'ROI mlt_range (hr): '+strjoin(string(mlt_range,format='(I0)'),','), log_file
    lprmsg, 'ROI rxy_range (Re): '+strjoin(string(rxy_range,format='(I0)'),','), log_file
    lprmsg, 'ROI Pdyn (nPa): '+string(pdyn,format='(I0)'), log_file
    lprmsg, 'ROI min duraiton (sec): '+string(roi_min_duration,format='(I0)'), log_file
    lprmsg, 'ROI min count (#): '+string(roi_min_count,format='(I0)'), log_file
    lprmsg, 'Probes searched: '+strjoin(strupcase(probes),','), log_file
    lprmsg, 'Time range searched: '+strjoin(time_string(time_range),','), log_file

    time_step = 60.
    common_times = make_bins(time_range, time_step)
    ntime = n_elements(common_times)

    secofday = constant('secofday')
    fillval = !values.f_nan
    days = make_bins(time_range, secofday)
    foreach probe, probes do begin
        lprmsg, '', log_file
        lprmsg, 'Processing '+strupcase(probe)+' ...', log_file

        prefix = probe+'_'
        mlt_var = prefix+'mlt'
        if ~check_if_update(mlt_var, time_range, dtime=time_step) then continue

        mlts = fltarr(ntime)+fillval
        foreach day, days do begin
            lprmsg, tab+'Processing '+time_string(day,tformat='YYYY-MM-DD')+' ...', log_file

            day_time_range = day+[0,secofday]
            time_index = where_pro(common_times, '[]', day_time_range, count=count)
            if count le roi_min_count then continue

            azim_dp_read_mlt, day_time_range, probe=probe, errmsg=errmsg
            if errmsg ne '' then continue

            get_data, mlt_var, times
            if n_elements(times) le roi_min_count then continue

            interp_time, mlt_var, common_times[time_index]
            mlts[time_index] = get_var_data(mlt_var)
        endforeach
        store_data, mlt_var, common_times, mlts
    endforeach


    mlt_counts = fltarr(ntime)
    foreach probe, probes do begin
        prefix = probe+'_'
        mlt_var = prefix+'mlt'
        mlts = get_var_data(mlt_var)
        if n_elements(mlts) ne ntime then continue
        index = where_pro(mlts, '[]', mlt_range, count=count)
        if count eq 0 then continue
        mlt_counts[index] += 1
    endforeach
    index = where(mlt_counts ge roi_min_count, count)
    if count eq 0 then return, roi_list
    roi_time_ranges = common_times[time_to_range(index,time_step=1)]

    durations = roi_time_ranges[*,1]-roi_time_ranges[*,0]
    index = where(durations ge roi_min_duration, count)
    if count eq 0 then return, roi_list


;---Now we break down the time_range into small pieces, use old method.
    roi_time_ranges = roi_time_ranges[index,*]
    nroi_time = count
    for roi_id=0, nroi_time-1 do begin
        roi_time_range = roi_time_ranges[roi_id,*]
        roi_list.add, /extract, azim_dp_search_roi(roi_time_range, $
            probes=probes, mlt_range=mlt_range, rxy_range=rxy_range, $
            pdyn=pdyn, min_dp_count=roi_min_count)
    endfor

;---Save the results.
    foreach roi, roi_list do azim_dp_roi_write, roi, filename=out_file
    lprmsg, '', log_file
    lprmsg, '# of ROI found: '+string(roi_list.length,format='(I0)') , log_file

    return, roi_list


end

time_range = time_double(['2008-01-01','2018-01-01'])
years = 2012
pad_time = 3600d
foreach year, years, year_id do begin
    time_range = time_double(string(year+[0,1],format='(I4)'))+[-1,1]*pad_time
    roi_list = azim_dp_vs_fac_search_roi(time_range, reset=1)
endforeach


end
