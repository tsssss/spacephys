;+
; Preprocess B field for the whole mission to:
;   1. Get the perturbation field B1
;   2. Get the background field B0
;-

pro pflux_grant_read_bfield_gen_file, time, probe=probe, filename=file, local_root=local_root;, pad_time=pad_time

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    if n_elements(file) eq 0 then begin
        errmsg = handle_error('No output file ...')
        return
    endif

    if n_elements(probe) eq 0 then begin
        errmsg = handle_error('No input probe ...')
        return
    endif

    if n_elements(time) eq 0 then begin
        errmsg = handle_error('No input time ...')
        return
    endif


;---Settings.
    secofday = constant('secofday')
    date = time[0]-(time[0] mod secofday)
    day_time_range = date+[0,secofday]
    window_sizes = [30d,2400]
    pad_time = max(window_sizes)    ; pad to remove boundary effect due to smooth.
    time_range = day_time_range+[-1,1]*pad_time
    project = pflux_grant_load_project()
    common_time_step = project.common_time_step
    common_times = make_bins(time_range, common_time_step)
    ncommon_time = n_elements(common_times)
    ndim = 3
    prefix = 'rbsp'+probe+'_'
    xyz = constant('xyz')
    fillval = !values.f_nan
    min_ntime = 1800d*64    ; 0.5 hour of data.

;---Load B field.
    b_gsm_var = prefix+'b_gsm'
    pflux_grant_read_b_gsm, time_range, probe=probe, local_root=local_root
    interp_time, b_gsm_var, common_times

;---Load position and get model B field.
    r_gsm_var = prefix+'r_gsm'
    r_gse_var = prefix+'r_gse'
    rbsp_read_orbit, time_range, probe=probe
    r_gse = get_var_data(r_gse_var, times=orbit_times)
    r_gsm = cotran(r_gse, orbit_times, 'gse2gsm')
    store_data, r_gsm_var, orbit_times, r_gsm
    add_setting, r_gsm_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'R', $
        'unit', 'Re', $
        'coord', 'GSM', $
        'coord_labels', xyz )

    par = 2.
    bmod_gsm_var = prefix+'bmod_gsm'
    r_gsm = get_var_data(r_gsm_var, at=orbit_times)
    bmod_gsm = float(r_gsm)
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
    bmod_gsm = sinterpol(bmod_gsm, orbit_times, common_times, /quadratic)
    store_data, bmod_gsm_var, common_times, bmod_gsm
    add_setting, bmod_gsm_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'T89 B', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz )


;---Separation B0 and B1.
    b0_gsm_var = prefix+'b0_gsm'
    b_gsm = get_var_data(b_gsm_var)
    bmod_gsm = get_var_data(bmod_gsm_var)

    db_gsm = b_gsm-bmod_gsm
    index = where(finite(snorm(b_gsm),/nan), count)
    if count ne 0 then db_gsm[index,*] = fillval

    ; Get dB in two bands.
    bg_vars = prefix+'db_bg_'+string(window_sizes,format='(I0)')+'_sec'
    foreach window_size, window_sizes,id do begin
        width = window_size/common_time_step
        bg = fltarr(ncommon_time,ndim)
        for ii=0,ndim-1 do begin
            tdat = db_gsm[*,ii]
            bg[*,ii] = smooth(tdat,width,/nan,/edge_mirror)
        endfor
        var = bg_vars[id]
        store_data, var, common_times, bg
    endforeach


    ; Cross-over ~3 Re, 2-4 Re.
    dis = snorm(get_var_data(r_gsm_var, times=orbit_times))
    dis = interpol(dis, orbit_times, common_times, /quadratic)
    store_data, prefix+'dis', common_times, dis
    weight_apogee = (tanh((dis-3)/(1*0.5))+1)*0.5
    weight_perigee = 1-weight_apogee
    bg_perigee_var = bg_vars[where(window_sizes eq min(window_sizes))]
    bg_perigee = get_var_data(bg_perigee_var)
    bg_apogee_var = bg_vars[where(window_sizes eq max(window_sizes))]
    bg_apogee = get_var_data(bg_apogee_var)
    b0_gsm = get_var_data(prefix+'bmod_gsm')
    for ii=0,ndim-1 do b0_gsm[*,ii] += bg_apogee[*,ii]*weight_apogee+bg_perigee[*,ii]*weight_perigee
    store_data, b0_gsm_var, common_times, b0_gsm
    add_setting, b0_gsm_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'B0', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz )


    b1_gsm_var = prefix+'b1_gsm'
    b_gsm = get_var_data(b_gsm_var)
    b0_gsm = get_var_data(b0_gsm_var)
    b1_gsm = b_gsm-b0_gsm

    weight_noise = (tanh((dis-1.5)/(1.0*0.5))+1)*0.5
    for ii=0,ndim-1 do b1_gsm[*,ii] *= weight_noise

    ; Remove spikes around perigee.
    index = where(snorm(b1_gsm) ge 30 and dis le 2.2, count)
    if count ne 0 then b1_gsm[index,*] = fillval


    store_data, b1_gsm_var, common_times, b1_gsm
    add_setting, b1_gsm_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'B1', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz )


;---Trim to day_time_range.
    index = where_pro(common_times, '[)', day_time_range)
    times = common_times[index]
    b0_gsm = get_var_data(b0_gsm_var)
    b0_gsm = float(b0_gsm[index,*])
    b1_gsm = get_var_data(b1_gsm_var)
    b1_gsm = float(b1_gsm[index,*])


;---Save data to file.
    cdf_id = (file_test(file))? cdf_open(file): cdf_create(file)

    time_var = 'ut_sec'
    time_var_settings = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix')
    cdf_save_var, time_var, value=times, filename=cdf_id
    cdf_save_setting, time_var_settings, filename=cdf_id, varname=time_var

    ; B1.
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'B1', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz)
    cdf_save_var, b1_gsm_var, value=b1_gsm, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=b1_gsm_var

    ; B0.
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'B0', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz)
    cdf_save_var, b0_gsm_var, value=b0_gsm, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=b0_gsm_var

    cdf_close, cdf_id

end

probe = 'b'
time = time_double('2012-12-04')    ; gap.
;time = time_double('2012-10-02')    ; spikes.
;time = time_double('2014-08-28')    ; spin tone.

file = join_path([homedir(),'test.cdf'])
pflux_grant_read_bfield_gen_file, time, probe=probe, filename=file

end
