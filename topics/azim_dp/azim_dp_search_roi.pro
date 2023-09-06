;+
; Search ROI in given time_range, probes, mlt_range, rxy_range.
;-

function azim_dp_search_roi, time_range, probes=probes, $
    mlt_range=mlt_range, rxy_range=rxy_range, pdyn=pdyn, min_dp_count=min_dp_count


;---Check input.
    roi_list = list()
    if n_elements(time_range) ne 2 then return, roi_list
    if n_elements(min_dp_count) eq 0 then min_dp_count = 4  ; #.
    if n_elements(probes) lt min_dp_count then return, roi_list
    if n_elements(mlt_range) ne 2 then return, roi_list
    if n_elements(rxy_range) ne 2 then rxy_range = [4.,30]  ; Re.
    if n_elements(pdyn) eq 0 then pdyn = 10.    ; nPa.

;---Search times in ROI.
    roi_min_duration = 3600.    ; sec.
    roi_min_count = min_dp_count
    probe_min_count = min_dp_count

    time_step = 60.     ; sec.
    common_times = make_bins(time_range, time_step)
    ncommon_time = n_elements(common_times)
    foreach probe, probes do begin
        prefix = probe+'_'
        roi_flag_var = prefix+'roi_flag'
        roi_flags = intarr(ncommon_time)
        store_data, roi_flag_var, common_times, roi_flags
        add_setting, roi_flag_var, /smart, dictionary($
            'display_type', 'scalar', $
            'yrange', [-0.1,1.1], $
            'short_name', strupcase(probe), $
            'unit', '#' )

        orbit_var = prefix+'r_sm'
        azim_dp_read_orbit, time_range, probe=probe
        if tnames(orbit_var) eq '' then begin
            have_data = 0
        endif else begin
            get_data, orbit_var, times, r_sm
            dis = snorm(r_sm)
            index = where(finite(dis), count)
            have_data = count ne 0
        endelse

        if ~have_data then continue
        roi_flags[*] = 1
        interp_time, orbit_var, common_times
        r_sm = get_var_data(orbit_var)

        ; rxy.
        rxy = snorm(r_sm[*,0:1])
        index = where_pro(rxy, '][', rxy_range, count=count)
        if count ne 0 then roi_flags[index] = 0

        ; magnetopause.
        r_gsm = cotran(r_sm, common_times, 'sm2gsm')
        index = where(check_if_in_magn(r_gsm, dynamic_pressure=pdyn) eq 0, count)
        if count ne 0 then roi_flags[index] = 0

        ; mlt.
        r_mag = cotran(r_sm, common_times, 'sm2mag')
        mlons = atan(r_mag[*,1],r_mag[*,0])*constant('deg')
        mlt = mlon2mlt(mlons, times)
        index = where_pro(mlt, '][', mlt_range, count=count)
        if count ne 0 then roi_flags[index] = 0
        index = where(finite(mlt,/nan), count)
        if count ne 0 then roi_flags[index] = 0


        ; Stay in ROI longer than roi_min_duration.
        index = where(roi_flags eq 1, count)
        if count ne 0 then begin
            time_ranges = common_times[time_to_range(index,time_step=1)]
            durations = time_ranges[*,1]-time_ranges[*,0]
            index = where(durations ge roi_min_duration, count)
            if count ne 0 then begin
                roi_flags[*] = 0
                time_ranges = time_ranges[index,*]
                ntime_range = n_elements(time_ranges)/2
                for time_id=0, ntime_range-1 do roi_flags[where_pro(common_times,'[]',reform(time_ranges[time_id,*]))] = 1
            endif
        endif

        store_data, roi_flag_var, common_times, roi_flags
    endforeach

;---Combine all flags into one.
    roi_flag_var = 'azim_dp_roi_flag'
    nprobe = n_elements(probes)
    roi_flags = intarr(ncommon_time, nprobe)
    foreach probe, probes, probe_id do begin
        prefix = probe+'_'
        roi_flags[*,probe_id] = get_var_data(prefix+'roi_flag')
    endforeach
    store_data, roi_flag_var, common_times, roi_flags, limits={$
        ytitle:'ROI flag (#)', yrange:[-0.1,1.1]}

;---Count of # of probes in ROI at each time tag.
    roi_counts = total(roi_flags,2)
    time_index = roi_counts ge probe_min_count
    index = where(time_index eq 1, count)
    if count eq 0 then return, roi_list
    times = common_times[index]

;---Rule out time ranges that are too short.
    time_ranges = time_to_range(times, time_step=time_step)
    durations = time_ranges[*,1]-time_ranges[*,0]
    index = where(durations gt roi_min_duration, count)
    if count eq 0 then return, roi_list
    time_ranges = time_ranges[index,*]

    ntime_range = n_elements(time_ranges)/2
    for ii=0, ntime_range-1 do begin
        time_range = reform(time_ranges[ii,*])
        if keyword_set(test_time) then if product(time_range-test_time) lt 0 then stop
        roi = dictionary($
            'time_range', time_range, $
            'duration', total(time_range*[-1,1])/60, $  ; in min.
            'nsection', 0, $
            'mlt_range', mlt_range, $
            'all_probes', !null, $
            'time_range_list', list(), $
            'probe_list', list())

        index = where_pro(common_times,'[]',time_range)
        the_roi_flags = roi_flags[index,*]
        the_times = common_times[index]
        probe_counts = total(the_roi_flags, 1)
        index = where(probe_counts ne 0, navailable_probe)
        if navailable_probe lt probe_min_count then message, 'Inconsistency, stop here ...'
        the_roi_flags = the_roi_flags[*,index]

        boundarys = time_range
        if navailable_probe gt probe_min_count then begin
            diff = fix(abs(the_roi_flags[1:-1,*]-the_roi_flags[0:-2,*]))
            diff = total(diff,2)
            index = where(diff ge 1, count)
            if count ne 0 then boundarys = [boundarys, the_times[index+1]]
        endif
        boundarys = sort_uniq(boundarys)
        nsection = n_elements(boundarys)-1
        time_range_list = list()
        probe_list = list()
        for jj=0, nsection-1 do begin
            the_time_range = boundarys[jj:jj+1]
            the_roi_flags = roi_flags[where_pro(common_times,'[)',the_time_range),*]
            probe_counts = total(the_roi_flags, 1)
            available_probes = probes[where(total(the_roi_flags,1) ne 0)]
            time_range_list.add, the_time_range
            probe_list.add, available_probes
        endfor
        roi.nsection = nsection
        roi.probe_list = probe_list
        roi.time_range_list = time_range_list

        ; Get all probes for all sections.
        all_section_probes = []
        for jj=0, nsection-1 do all_section_probes = [all_section_probes,roi.probe_list[jj]]
        all_section_probes = sort_uniq(all_section_probes)
        roi['all_probes'] = all_section_probes

        roi_list.add, roi
    endfor

    return, roi_list

end
