;+
; Test to preprocess RBSP B field for a given time range and probe.
;-

pro test_preprocess_bfield, time_range, probe=probe

    orbit_time_ranges = pflux_grant_calculate_orbit_time_range(time_range, probe=probe)
    norbit = n_elements(orbit_time_ranges)/2
    if norbit eq 0 then stop
    bfield_time_step = 4.
    bfield_resolution = '4sec'
    rbsp_read_bfield, time_range, probe=probe, resolution=bfield_resolution
    bfield_times = make_bins(time_range, bfield_time_step)
    bfield_ntime = n_elements(bfield_times)
    ndim = 3
    par = 2.
    model = 't89'
    xyz = constant('xyz')
    rgb = constant('rgb')

    ; Calculate B model and dB.
    prefix = 'rbsp'+probe+'_'
    r_var = prefix+'r_gsm'
    r_gsm = get_var_data(r_var, at=bfield_times)
    b0_gsm = fltarr(bfield_ntime,ndim)
    for ii=0, bfield_ntime-1 do begin
        tilt = geopack_recalc(bfield_times[ii])
        ; in-situ position
        rx = r_gsm[ii,0]
        ry = r_gsm[ii,1]
        rz = r_gsm[ii,2]
        ; in-situ B field.
        geopack_igrf_gsm, rx,ry,rz, bx,by,bz
        geopack_t89, par, rx,ry,rz, dbx,dby,dbz
        b0_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
    endfor
    b0_var = prefix+'b0_gsm'
    store_data, b0_var, bfield_times, b0_gsm
    add_setting, b0_var, /smart, {$
        display_type: 'vector', $
        unit: 'nT', $
        short_name: 'B0!S!UT89!N!R', $
        coord: 'GSM', $
        coord_labels: xyz}

    b_var = prefix+'b_gsm'
    b_gsm = get_var_data(b_var, at=bfield_times)
    db_gsm = b_gsm-b0_gsm
    db_var = prefix+'db_gsm'
    store_data, db_var, bfield_times, db_gsm
    add_setting, db_var, /smart, {$
        display_type: 'vector', $
        unit: 'nT', $
        short_name: 'dB', $
        coord: 'GSM', $
        coord_labels: xyz}

    dbmag = snorm(db_gsm)
    dbmag_var = prefix+'dbmag'
    store_data, dbmag_var, bfield_times, dbmag
    add_setting, dbmag_var, /smart, {$
        display_type: 'scalar', $
        unit: 'nT', $
        short_name: '|dB|'}


    width1 = 3600d/bfield_time_step
    width2 = 10.
    perigee_lshell = 4.
    shell_width = 600.
    new_var = prefix+'b1_gsm'
    del_data, new_var
    for ii=0, norbit-1 do begin
        orbit_time_range = reform(orbit_time_ranges[ii,*])
        dis = snorm(get_var_data(r_var, in=orbit_time_range, times=times))
        perigee_lshell_times = times[minmax(where(dis ge perigee_lshell))]
        shell_times = perigee_lshell_times+[-1,1]*shell_width

        db_gsm = get_var_data(db_var, in=orbit_time_range, times=times, limits=lim)
        db1_gsm = db_gsm
        for jj=0, ndim-1 do db1_gsm[*,jj] = smooth(db1_gsm[*,jj], width1)
        db2_gsm = db_gsm
        for jj=0, ndim-1 do db2_gsm[*,jj] = smooth(db2_gsm[*,jj], width2)

        ; Construct a weight which is 1 outside 4 Re, and 0 around perigee.
        ; The width is w, and the exponential decay is right at w away from the time at 4 Re.
        ; The weight is used to calculate the background of dB.
        ntime = n_elements(times)
        weight = fltarr(ntime)
        middle_index = ntime/2
        weight[0:middle_index-1] = (tanh((times[0:middle_index-1]-shell_times[0])/shell_width)+1)*0.5
        weight[middle_index:-1] = (1-tanh((times[middle_index:-1]-shell_times[1])/shell_width))*0.5

        ; Suppress field within 4 Re.
        for jj=0, ndim-1 do begin
            db_gsm[*,jj] -= db1_gsm[*,jj]*weight+db2_gsm[*,jj]*(1-weight)
        endfor
        
        ; Remove spikes within 4 Re.
        index = convert_time

        ; Construct a weight, with a width of half of the perigee length
        weight2 = fltarr(ntime)
        decay_times = (orbit_time_range+perigee_lshell_times)*0.5
        decay_width = min(abs(orbit_time_range-decay_times))*0.8
        weight2[0:middle_index-1] = (tanh((times[0:middle_index-1]-decay_times[0])/decay_width)+1)*0.5
        weight2[middle_index:-1] = (1-tanh((times[middle_index:-1]-decay_times[1])/decay_width))*0.5
        for jj=0, ndim-1 do db_gsm[*,jj] *= weight2

        if tnames(new_var) eq '' then begin
            store_data, new_var, times, db_gsm, limits=lim
        endif else begin
            yys = get_var_data(new_var, times=xxs)
            store_data, new_var, [xxs,times], [yys,db_gsm]
        endelse
    endfor
end

; Low frequ fluctations around perigee.
time_range = time_double(['2014-02-20','2014-02-23'])
probe = 'a'
time_range = time_double(['2013-06-06','2013-06-09'])
probe = 'a'
test_preprocess_bfield, time_range, probe=probe
end
