;+
; Search for time periods of good triangles.
;-

pro azim_df_search_themis_good_triangle, project=project

;---Settings.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    the_key = 'filter_storm_setting'
    if ~project.haskey(the_key) then azim_df_filter_storm, project
    filter_storm_setting = project[the_key]

    the_key = 'themis_search_time_range'
    if ~project.haskey(the_key) then project[the_key] = time_double(['2007-03-01','2009-09-01'])
    time_range = project[the_key]
    the_key = 'themis_search_time_step'
    if ~project.haskey(the_key) then project[the_key] = 5*60. ; 1 min cadence runs way too slow.
    time_step = project[the_key]
    mission_probes = 'th'+letters('e')
    common_times = make_bins(time_range, time_step)

;---Prepare thx_r_sm on common_times.
    foreach mission_probe, mission_probes do begin
        probe_info = resolve_probe(mission_probe)
        prefix = probe_info.prefix
        probe = probe_info.probe
        if check_if_update(prefix+'r_gsm', time_range) then begin
            themis_read_orbit, time_range, probe=probe
        endif

        the_var = prefix+'r_sm'
        if check_if_update(the_var, time_range) then begin
            get_data, prefix+'r_gsm', times, rgsm
            rsm = cotran(rgsm, times, 'gsm2sm')
            store_data, the_var, times, rsm
            add_setting, the_var, /smart, {$
                display_type: 'vector', $
                short_name: 'R', $
                unit: 'Re', $
                coord: 'SM', $
                coord_labels: ['x','y','z']}
            interp_time, the_var, common_times
        endif

        deg = 180d/!dpi
        the_var = prefix+'mlt'
        if check_if_update(the_var, time_range) then begin
            get_data, prefix+'r_sm', times, rsm
            rmag = cotran(rsm, times, 'sm2mag')
            mlat = asin(rmag[*,2]/snorm(rmag))*deg
            mlon = atan(rmag[*,1],rmag[*,0])*deg
            mlt = mlon2mlt(mlon, times)
            store_data, the_var, common_times, mlt
            add_setting, the_var, /smart, {$
                display_type: 'scalar', $
                short_name: 'MLT', $
                unit: 'hr'}
        endif
    endforeach

;---Search for times when thb,thc were in the same side of midnight.
    possible_times = common_times
    get_data, 'thb_mlt', common_times, thb_mlt
    get_data, 'thc_mlt', common_times, thc_mlt
    index = where(thb_mlt*thc_mlt ge 0, count)
    if count eq 0 then message, 'No possible time, stop here ...'
    possible_times = possible_times[index]
    thb_mlt = thb_mlt[index]
    thc_mlt = thc_mlt[index]


;---Need to be within 7 hr around midnight.
    mlt_range = filter_storm_setting.mlt_range
    index = where_pro(thb_mlt, '[]', mlt_range, count=count)
    if count eq 0 then message, 'No possible time, stop here ...'
    possible_times = possible_times[index]
    thb_mlt = thb_mlt[index]
    thc_mlt = thc_mlt[index]
    index = where_pro(thc_mlt, '[]', mlt_range, count=count)
    if count eq 0 then message, 'No possible time, stop here ...'
    possible_times = possible_times[index]
    thb_mlt = thb_mlt[index]
    thc_mlt = thc_mlt[index]

;---Need to be within the magnetopause.
    flags = check_if_in_magn(get_var_data('thb_r_gsm', at=possible_times))
    index = where(flags eq 1, count)
    if count eq 0 then message, 'No possible time, stop here ...'
    possible_times = possible_times[index]
    thb_mlt = thb_mlt[index]
    thc_mlt = thc_mlt[index]
    flags = check_if_in_magn(get_var_data('thc_r_gsm', at=possible_times))
    index = where(flags eq 1, count)
    if count eq 0 then message, 'No possible time, stop here ...'
    possible_times = possible_times[index]
    thb_mlt = thb_mlt[index]
    thc_mlt = thc_mlt[index]

;---Need to be outside 15 Re.
    min_dis = 15.
    flags = snorm(get_var_data('thb_r_gsm', at=possible_times)) ge min_dis
    index = where(flags eq 1, count)
    if count eq 0 then message, 'No possible time, stop here ...'
    possible_times = possible_times[index]
    thb_mlt = thb_mlt[index]
    thc_mlt = thc_mlt[index]
    flags = snorm(get_var_data('thc_r_gsm', at=possible_times)) ge min_dis
    index = where(flags eq 1, count)
    if count eq 0 then message, 'No possible time, stop here ...'
    possible_times = possible_times[index]
    thb_mlt = thb_mlt[index]
    thc_mlt = thc_mlt[index]


;---Search for times of "good triad" among thb,thc,thx.
    npossible_time = n_elements(possible_times)
    num_good_triads = bytarr(npossible_time)
    other_probes = ['a','d','e']
    small_angle_limit = 15. ; deg.
    angle_range = [small_angle_limit,180-small_angle_limit]
    ndim = 3
    for ii=0, npossible_time-1 do begin
        lprmsg, 'Processing '+time_string(possible_times[ii])+' ...'
        the_time = possible_times[ii]
        num_good_triad = 0
        foreach probe, other_probes do begin
            prefix = 'th'+probe+'_'
            ; Must be in magnetopause.
            if ~check_if_in_magn(get_var_data(prefix+'r_gsm', at=the_time)) then continue
            ; Must be in the MLT range.
            the_mlt = get_var_data(prefix+'mlt', at=the_time)
            index = where_pro(the_mlt, '[]', mlt_range, count=count)
            if count eq 0 then continue
            ; Must be on the same side of midnight.
            if the_mlt*thb_mlt[ii] lt 0 then continue
            xxs = fltarr(ndim)
            yys = fltarr(ndim)
            angles = fltarr(ndim)
            combo = ['b','c',probe]
            for jj=0, ndim-1 do begin
                rsm = get_var_data('th'+combo[jj]+'_r_sm', at=the_time)
                xxs[jj] = rsm[0]
                yys[jj] = rsm[1]
            endfor
            for jj=0,ndim-1 do begin
                xxs = shift(xxs,1)
                yys = shift(yys,1)
                angles[jj] = sang($
                    [xxs[1]-xxs[0],yys[1]-yys[0]], $
                    [xxs[2]-xxs[0],yys[2]-yys[0]], /degree)
            endfor
            index = where_pro(angles, '[]', angle_range, count=count)
            if count eq ndim then begin
                num_good_triad += 1
                lprmsg, 'Found a good triangle with TH-'+strupcase(probe)+' ...'
            endif
        endforeach
        num_good_triads[ii] = num_good_triad
    endfor


;---Get the time ranges of candidates.
    min_num_triad = 1
    index = where(num_good_triads ge min_num_triad, count)
    if count eq 0 then message, 'No possible time, stop here ...'
    good_triad_times = possible_times[index]
    pad_times = [-1,1]*time_step
    time_ranges = time_to_range(good_triad_times, time_step=time_step, pad_times=pad_times)
    ntime_range = n_elements(time_ranges)/2


;---Remove short candidates.
    min_triad_duration = 60*60.     ; sec.
    flags = bytarr(ntime_range)
    for ii=0, ntime_range-1 do begin
        time_range = time_ranges[ii,*]
        flags[ii] = total(time_range*[-1,1]) ge min_triad_duration
    endfor
    index = where(flags eq 1, ntime_range)
    if ntime_range eq 0 then message, 'No candidate, stop here ...'
    time_ranges = time_ranges[index,*]


;---Write result to a file.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    log_file = join_path([project.data_dir,'azim_df_search_themis_good_triangle_list.txt'])
    openw, lun, /get_lun, log_file
    for ii=0, ntime_range-1 do begin
        printf, lun, strjoin(reform(time_string(time_ranges[ii,*])),' to ')
    endfor
    free_lun, lun

end
