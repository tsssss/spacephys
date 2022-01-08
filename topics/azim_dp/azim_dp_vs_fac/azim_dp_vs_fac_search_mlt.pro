;+
; Search candidates with certain MLT distribution.
;-

function azim_dp_vs_fac_search_mlt, roi_list, $
    reset=reset


;---Check input.
    mlt_list = list()
    if n_elements(roi_list) eq 0 then return, mlt_list
    if n_elements(mlt_min_extent) eq 0 then mlt_min_extent = 5. ;hr
    if n_elements(mlt_max_sep) eq 0 then mlt_max_sep = 3.   ; hr
    ;mlt_max_avg_sep = 1.5   ; hr, mlt avg separation must be smaller than this.
    roi_min_duration = 3600.    ; sec.
    min_probe_count = 6


;---Settings.
    root_dir = join_path([googledir(),'works','works','azim_dp_vs_fac','data'])
    if file_test(root_dir) eq 0 then file_mkdir, root_dir
    log_file = join_path([root_dir,'azim_dp_vs_fac_search_mlt.log'])
    out_file = join_path([root_dir,'azim_dp_vs_fac_search_mlt.txt'])
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
        mlt_list = azim_dp_candidate_read(filename=out_file)
        return, mlt_list
    endif

    foreach file, files do if file_test(file) eq 0 then ftouch, file


;---Search for candidate with proper MLT distribution.
    tab = constant('4space')
    time_step = 60. ; sec.

    lprmsg, 'Search candidate_list for azim_dp_vs_fac', log_file
    lprmsg, '', log_file
    lprmsg, 'MLT '


    foreach roi, roi_list do begin


        for section_id=0, roi.nsection-1 do begin
            time_range = roi.time_range_list[section_id]
            probes = roi.probe_list[section_id]

            duration = total(time_range*[-1,1])
            if duration le roi_min_duration then continue


            ; Load MLT.
            times = make_bins(time_range, time_step)
            ntime = n_elements(times)
            nprobe = n_elements(probes)
            mlts = fltarr(ntime,nprobe)
            foreach probe, probes, probe_id do begin
                prefix = probe+'_'
                mlt_var = prefix+'mlt'
                azim_dp_read_mlt, time_range, probe=probe, errmsg=errmsg
                interp_time, mlt_var, times
                mlts[*,probe_id] = get_var_data(mlt_var)
            endforeach

            ; Figure out the flag for each time.
            flags = fltarr(ntime)
            mlt_ranges = fltarr(ntime,2)
            for time_id=0,ntime-1 do begin
                the_mlts = reform(mlts[time_id,*])
                index = sort(the_mlts)
                the_mlts = the_mlts[index]
                the_diffs = the_mlts[1:nprobe-1]-the_mlts[0:nprobe-2]

                ; Remove outlier.
                if the_diffs[0] ge mlt_max_sep then begin
                    the_mlts = the_mlts[1:*]
                    the_diffs = the_diffs[1:*]
                endif
                if the_diffs[-1] ge mlt_max_sep then begin
                    the_mlts = the_mlts[0:-2]
                    the_diffs = the_diffs[0:-2]
                endif

                ndata = n_elements(the_mlts)
                if ndata lt min_probe_count then continue

                ; Check mean, extent, etc.
                the_dmlt_max = max(the_diffs)
                the_mlt_del = the_mlts[-1]-the_mlts[0]
                the_dmlt_sep = mean(the_diffs)
                if the_dmlt_max ge mlt_max_sep then continue
                if the_mlt_del lt mlt_min_extent then continue

                mlt_ranges[time_id,*] = minmax(the_mlts)
                flags[time_id] = 1
            endfor
            index = where(flags eq 1, count)
            if count eq 0 then continue

            time_ranges = times[time_to_range(index,time_step=1)]
            durations = time_ranges[*,1]-time_ranges[*,0]
            index = where(durations ge roi_min_duration, count)
            if count eq 0 then continue
            time_ranges = time_ranges[index,*]

            for ii=0,count-1 do begin
                the_time_range = reform(time_ranges[ii,*])
                index = lazy_where(times, '[]', the_time_range)
                the_mlt_range = minmax(mlt_ranges[index,*])
                candidate = dictionary($
                    'time_range', the_time_range, $
                    'mlt_range', the_mlt_range, $
                    'probes', probes)

                mlt_list.add, candidate
            endfor
        endfor
    endforeach

    return, mlt_list

end

time_range = time_double(['2014-08-28','2014-08-29'])
;time_range = time_double(['2016-10-13','2016-10-14'])
time_range = time_double(['2013-06-07','2013-06-08'])
roi_list = azim_dp_vs_fac_search_roi(time_range)
;foreach roi, roi_list do print, time_string(roi.time_range), n_elements(roi.all_probes)
mlt_list = azim_dp_vs_fac_search_mlt(roi_list)
foreach tmp, mlt_list do print, time_string(tmp.time_range), strjoin(tmp.probes,','), tmp.mlt_range
end
