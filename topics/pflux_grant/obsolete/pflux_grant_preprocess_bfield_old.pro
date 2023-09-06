;+
; Preprocess B field for a given time range.
;-

pro pflux_grant_preprocess_bfield, time_range, probe=probe

    ; Decompose into orbits.
    orbit_time_ranges = pflux_grant_calculate_orbit_time_range(time_range, probe=probe)
    norbit = n_elements(orbit_time_ranges)/2
    if norbit eq 0 then return
    orbit_period = mean(orbit_time_ranges[*,1]-orbit_time_ranges[*,0])

    prefix = 'rbsp'+probe+'_'
    bfield_time_step = 1d/64
    bfield_resolution = 'hires'
    rbsp_read_bfield, time_range, probe=probe, resolution=bfield_resolution

    ; Load B field data.
    b_var = prefix+'b_gsm'
    bfield_min_data_ratio = 0.8
    bfield_min_nrec = bfield_min_data_ratio*orbit_period/bfield_time_step
    if tnames(b_var) eq '' then return

    ; Load R data.
    r_var = prefix+'r_gsm'
    rbsp_read_orbit, time_range, probe=probe

    ; Get B model.
    orbit_time_step = 60.
    ndim = 3
    orbit_times = make_bins(time_range, orbit_time_step)
    norbit_time = n_elements(orbit_times)
    r_gsm = get_var_data(r_var, at=orbit_times)
    bmod_gsm = fltarr(norbit_time,ndim)
    model = 't89'
    par = 2.
    xyz = constant('xyz')
    foreach time, orbit_times, ii do begin
        tilt = geopack_recalc(time)
        rx = r_gsm[ii,0]
        ry = r_gsm[ii,1]
        rz = r_gsm[ii,2]
        ; in-situ B field.
        geopack_igrf_gsm, rx,ry,rz, bx,by,bz
        geopack_t89, par, rx,ry,rz, dbx,dby,dbz
        bmod_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
    endforeach
    bmod_gsm_var = prefix+'bmod_gsm'
    store_data, bmod_gsm_var, orbit_times, bmod_gsm
    add_setting, bmod_gsm_var, /smart, {$
        display_type: 'vector', $
        unit: 'nT', $
        short_name: 'B0!S!U'+strupcase(model)+'!N!R', $
        coord: 'GSM', $
        coord_labels: xyz}


;---Loop through orbits.
    common_time_step = 1d/16
    common_times = make_bins(time_range, common_time_step)
    ncommon_time = n_elements(common_times)
    common_b1_gsm = fltarr(ncommon_time,ndim)
    perigee_shell = 4.
    shell_width = 600.  ; sec
    width1 = 3600d/common_time_step
    width2 = 30d/common_time_step
    for orbit_id=0, norbit-1 do begin
        orbit_time_range = reform(orbit_time_ranges[orbit_id,*])
        data = get_var_data(b_var, in=orbit_time_range)
        if n_elements(data) eq 0 then continue
        ndata = n_elements(data[*,0])
        if ndata lt bfield_min_nrec then continue

        ; Get dB = B_meas - B_model.
        time_index = where_pro(common_times, '[]', orbit_time_range)
        times = common_times[time_index]
        b_gsm = get_var_data(b_var, at=times)
        bmod_gsm = get_var_data(bmod_gsm_var, at=times)
        db_gsm = b_gsm-bmod_gsm

        ; Smooth dB in two bands.
        db1_gsm = db_gsm
        for jj=0, ndim-1 do db1_gsm[*,jj] = smooth(db1_gsm[*,jj], width1)
        db2_gsm = db_gsm
        for jj=0, ndim-1 do db2_gsm[*,jj] = smooth(db2_gsm[*,jj], width2)

        ; Weight the two smoothed dB to get the background dB.
        dis = snorm(get_var_data(r_var, at=times))
        perigee_shell_time = times[minmax(where(dis ge perigee_shell))]
        shell_times = perigee_shell_time+[-1,1]*shell_width
        ntime = n_elements(times)
        weight = fltarr(ntime)
        middle_index = ntime/2
        weight[0:middle_index-1] = (tanh((times[0:middle_index-1]-shell_times[0])/shell_width)+1)*0.5
        weight[middle_index:-1] = (1-tanh((times[middle_index:-1]-shell_times[1])/shell_width))*0.5
        for jj=0, ndim-1 do begin
            db_gsm[*,jj] -= db1_gsm[*,jj]*weight+db2_gsm[*,jj]*(1-weight)
        endfor

        ; Weight down the dB within perigee shell.
        weight2 = fltarr(ntime)
        decay_times = (orbit_time_range+perigee_shell_time)*0.5
        decay_width = min(abs(orbit_time_range-decay_times))*0.8
        weight2[0:middle_index-1] = (tanh((times[0:middle_index-1]-decay_times[0])/decay_width)+1)*0.5
        weight2[middle_index:-1] = (1-tanh((times[middle_index:-1]-decay_times[1])/decay_width))*0.5
        for jj=0, ndim-1 do db_gsm[*,jj] *= weight2

        ; Remove spikes within perigee shell.
        perigee_index = where(dis lt perigee_shell)
        perigee_db_gsm = db_gsm[perigee_index,*]
        perigee_times = times[perigee_index]
        bmag = snorm(perigee_db_gsm)
        bmag -= smooth(bmag, width2*10, /edge_truncate)
        max_bmag = 5.
        index = where(abs(bmag) ge max_bmag, count, complement=good_index)
        if count ne 0 then begin
            perigee_db_gsm = sinterpol(perigee_db_gsm[good_index,*], perigee_times[good_index], perigee_times)
            db_gsm[perigee_index,*] = perigee_db_gsm
        endif

        common_b1_gsm[time_index,*] = db_gsm
    endfor

    b1_gsm_var = prefix+'b1_gsm'
    store_data, b1_gsm_var, common_times, common_b1_gsm
    add_setting, b1_gsm_var, /smart, {$
        display_type: 'vector', $
        unit: 'nT', $
        short_name: 'B1', $
        coord: 'GSM', $
        coord_labels: xyz}

    stop
end

time_range = time_double(['2014-02-21','2014-02-23'])
probe = 'a'
time_range = time_double(['2013-06-06','2013-06-09'])
probe = 'b'
pflux_grant_preprocess_bfield, time_range, probe=probe
end
