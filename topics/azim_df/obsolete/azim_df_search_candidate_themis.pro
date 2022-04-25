;+
; Search for time periods of good triangles.
;-

pro azim_df_search_candidate_themis, project=project

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


;---Remove data outside the region of interest.
    ncommon_time = n_elements(common_times)
    common_flags = bytarr(ncommon_time)+1
    wanted_probes = 'th'+['b','c']
    wanted_min_dis = 15.    ; Re.
    mlt_range = filter_storm_setting.mlt_range
    nwanted_probe = n_elements(wanted_probes)
    mlt_flags = fltarr(ncommon_time)+1
    foreach mission_probe, wanted_probes do begin
        prefix = mission_probe+'_'

        ; Within the magnetopause.
        the_var = prefix+'r_gsm'
        rgsm = get_var_data(the_var, at=common_times)
        flags = check_if_in_magn(rgsm)
        index = where(flags eq 0, count)
        if count ne 0 then common_flags[index] = 0

        ; Beyond min_dis.
        dis = snorm(rgsm)
        flags = where(dis lt wanted_min_dis, count)
        if count ne 0 then common_flags[index] = 0

        ; Around midnight.
        the_var = prefix+'mlt'
        mlt = get_var_data(the_var, at=common_times)
        mlt_flags *= mlt
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then common_flags[index] = 0
    endforeach


;---Wanted probes must be on the same side of midnight.
    index = where(mlt_flags le 0, count)
    if count ne 0 then common_flags[index] = 0


;---Search for times of "good triad" among thb,thc,thx.
    index = where(common_flags eq 1, npossible_time)
    if npossible_time eq 0 then message, 'No possible time, stop here ...'
    possible_times = common_times[index]
    other_probes = 'th'+['a','d','e']
    nother_probe = n_elements(other_probes)
    data_flags = bytarr(npossible_time, nother_probe)+1
    other_min_dis = 4.
    foreach mission_probe, other_probes, ii do begin
        prefix = mission_probe+'_'
        ; Within the magnetopause.
        the_var = prefix+'r_gsm'
        rgsm = get_var_data(the_var, at=possible_times)
        flags = check_if_in_magn(rgsm)
        index = where(flags eq 0, count)
        if count ne 0 then data_flags[index,ii] = 0

        ; Beyond min_dis.
        dis = snorm(rgsm)
        index = where(dis lt other_min_dis, count)
        if count ne 0 then data_flags[index,ii] = 0

        ; Around midnight.
        the_var = prefix+'mlt'
        mlt = get_var_data(the_var, at=possible_times)
        mlt_flags *= mlt
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then data_flags[index,ii] = 0
    endforeach


    ; Prepare combined r_sm and mlt to make later calculation faster.
    ndim = 3
    nprobe = n_elements(mission_probes)
    all_rsm = fltarr(npossible_time, ndim, nprobe)
    all_mlt = fltarr(npossible_time, nprobe)
    foreach mission_probe, mission_probes, ii do begin
        prefix = mission_probe+'_'
        all_rsm[*,*,ii] = get_var_data(prefix+'r_sm', at=possible_times)
        all_mlt[*,ii] = get_var_data(prefix+'mlt', at=possible_times)
    endforeach
    thb_index = where(mission_probes eq 'thb')
    thc_index = where(mission_probes eq 'thc')
    other_indexs = intarr(nother_probe)
    foreach mission_probe, other_probes, ii do other_indexs[ii] = where(mission_probes eq mission_probe)

    ; Search for good triads.
    triad_angle_range = project.triad_angle_range
    num_good_triads = bytarr(npossible_time)
    xxs = fltarr(ndim)
    yys = fltarr(ndim)
    angles = fltarr(ndim)
    for ii=0, npossible_time-1 do begin
        lprmsg, 'Processing '+time_string(possible_times[ii])+' ...'
        ;if possible_times[ii] eq time_double('2008-01-03/13:10') then stop
        num_good_triad = 0
        xxs[0:1] = all_rsm[ii,0,[thb_index,thc_index]]
        yys[0:1] = all_rsm[ii,1,[thb_index,thc_index]]
        foreach probe, other_probes, jj do begin
            prefix = 'th'+probe+'_'
            ; Have data.
            if data_flags[ii,jj] eq 0 then continue
            ; Same side of midnight.
            if all_mlt[ii,other_indexs[jj]]*all_mlt[ii,thb_index] le 0 then continue
            ; Calculate triad.
            xxs[2] = all_rsm[ii,0,other_indexs[jj]]
            yys[2] = all_rsm[ii,1,other_indexs[jj]]
            for kk=0,ndim-1 do begin
                xxs = shift(xxs,1)
                yys = shift(yys,1)
                angles[kk] = sang($
                    [xxs[1]-xxs[0],yys[1]-yys[0]], $
                    [xxs[2]-xxs[0],yys[2]-yys[0]], /degree)
            endfor
            index = lazy_where(angles, '[]', triad_angle_range, count=count)
            if count eq ndim then begin
                num_good_triad += 1
            endif
        endforeach
        num_good_triads[ii] = num_good_triad
        if num_good_triad ne 0 then lprmsg, 'Found '+string(num_good_triad,format='(I0)')+' good triad(s) ...'
    endfor



;---Get the time ranges of candidates.
    min_num_triad = 3
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
    log_file = join_path([project.data_dir,'azim_df_search_themis_candidate.txt'])
    openw, lun, /get_lun, log_file
    for ii=0, ntime_range-1 do begin
        printf, lun, strjoin(reform(time_string(time_ranges[ii,*])),' to ')
    endfor
    free_lun, lun



;---Filter according to AE.
    time_range = project['themis_search_time_range']
    time_step = project['themis_search_time_step']
    common_times = make_bins(time_range, time_step)
    candidate_list = log_file

    lines = read_all_lines(candidate_list)
    candidates = list()
    foreach line, lines do begin
        info = strsplit(line,' ',/extract)
        candidates.add, time_double(info[[0,2]])
    endforeach

    if check_if_update('ae', time_range) then omni_read_index, time_range
    all_ae = get_var_data('ae', at=common_times)

    ae_flags = list()
    min_ae = 500.   ; nT.
    min_ae_duration = 30*60.    ; sec.
    events = list()
    foreach the_time_range, candidates do begin
        the_flag = 0
;        index = lazy_where(common_times, '[]', the_time_range+[-1,1]*min_ae_duration, count=count)
        index = lazy_where(common_times, '[]', the_time_range, count=count)
        if count eq 0 then begin
            ae_flags.add, the_flag
            continue
        endif
        the_ae = all_ae[index]
        the_times = common_times[index]
        index = where(the_ae ge min_ae, count)
        if count le 1 then begin
            ae_flags.add, the_flag
            continue
        endif

; test_time = time_double('2008-03-08/13:05')
; if test_time gt the_time_range[0] and test_time lt the_time_range[1] then stop

        time_ranges = time_to_range(the_times[index], time_step=time_step)
        ntime_range = n_elements(time_ranges)
        durations = time_ranges[*,1]-time_ranges[*,0]
        index = where(durations ge min_ae_duration, count)
        if count eq 0 then begin
            ae_flags.add, the_flag
            continue
        endif

        the_flag = 1
        ae_flags.add, the_flag
        for ii=0, count-1 do events.add, reform(time_ranges[index[ii],*])
    endforeach

    ae_flags = ae_flags.toarray()
    index = where(ae_flags eq 1, count)
    if count eq 0 then message, 'No candidate found, stop here ...'
    time_ranges = events.toarray()   ; in [N,2]
    ntime_range = n_elements(time_ranges)/2

;---Write result to a file.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    candidate_list = join_path([project.data_dir,'azim_df_candidate_themis.txt'])
    openw, lun, /get_lun, candidate_list
    for ii=0, ntime_range-1 do begin
        printf, lun, strjoin(reform(time_string(time_ranges[ii,*])),' to ')
    endfor
    free_lun, lun

    project['done_search_candidate_themis'] = 1
    project['candidate_suffix_themis'] = fgetbase(candidate_list)
    update_project, project

end
