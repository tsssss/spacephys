;+
; Test to preprocess RBSP B field for a given time range and probe.
;-

pro test_pflux_grant_preprocess_ebfield, probe=probe, project=project

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project.pflux_calc_setting

    if n_elements(probe) eq 0 then return
    if (where(probe eq ['a','b']))[0] eq -1 then return
    rbspx = 'rbsp'+probe
    secofday = constant('secofday')
    mission_time_range = pflux_calc_setting[rbspx].time_range+[-1,1]*secofday
    nday = total(mission_time_range*[-1,1])/secofday
    dates = mission_time_range[0]+findgen(nday)*secofday


;---Some settings.
    common_time_step = 1d/8
    orbit_period = 9.*constant('secofhour')
    bfield_resolution = 'hires'
    bfield_min_data_ratio = 0.8
    bfield_min_nrec = bfield_min_data_ratio*orbit_period/common_time_step
    orbit_time_step = 60.
    ndim = 3
    model = 't89'
    par = 2.
    xyz = constant('xyz')
    perigee_shell = 4.  ; Use data beyond 4 Re.
    window_large = 3600d    ; in sec, to get B0.
    spin_period = 10.   ; indeed 11-12 sec.
    prefix = rbspx+'_'

    local_root = join_path([default_local_root(),'sdata','rbsp'])
    paths = [rbspx,'ebfield2']


    ;time_range = time_double(['2012-09-29/21:12:00','2012-10-01/00:09:00'])
    time_range = time_double(['2013-03-01/06:30','2013-03-01/09:46'])
    ;time_range = time_double(['2016-06-01/08:42','2016-06-01/11:59'])
    ;time_range = time_double(['2017-09-01/06:09','2017-09-01/09:26'])
    data_time_range = time_range+[-1,1]*window_large


    ; Common times.
    common_times = make_bins(data_time_range, common_time_step)
    ncommon_time = n_elements(common_times)

    ; Load B field. GSM contain weird spin modulation, use UVW is the same.
    rbsp_read_bfield, data_time_range, probe=probe, resolution=bfield_resolution
    b_gsm_var = prefix+'b_gsm'
    interp_time, b_gsm_var, common_times


    ; Load R data.
    r_gsm_var = prefix+'r_gsm'
    rbsp_read_orbit, data_time_range, probe=probe
    interp_time, r_gsm_var, common_times


    ; Get B model.
    r_gsm = get_var_data(r_gsm_var)
    bmod_gsm = fltarr(ncommon_time,ndim)
    foreach time, common_times, ii do begin
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
    store_data, bmod_gsm_var, common_times, bmod_gsm
    add_setting, bmod_gsm_var, /smart, {$
        display_type: 'vector', $
        unit: 'nT', $
        short_name: 'B0!S!U'+strupcase(model)+'!N!R', $
        coord: 'GSM', $
        coord_labels: xyz}


    bmod_gsm = get_var_data(bmod_gsm_var)
    b_gsm = get_var_data(b_gsm_var)
    db_gsm = b_gsm-bmod_gsm


stop

    width1 = 3600d/common_time_step   ; sec.
    width2 = 30d/common_time_step
    db1_gsm = db_gsm
    for jj=0, ndim-1 do db1_gsm[*,jj] = smooth(db1_gsm[*,jj], width1, /nan)
    db2_gsm = db_gsm
    for jj=0, ndim-1 do db2_gsm[*,jj] = smooth(db2_gsm[*,jj], width2, /nan)

    dis = snorm(get_var_data(r_gsm_var, at=common_times))
    weight = fltarr(ncommon_time)
    center_time = mean(time_range)
    max_width = max([width1,width2])
    weight_window = total(time_range*[-1,1])*0.5-max_width*common_time_step
    z = (common_times-center_time)/(weight_window*0.5)
    weight = exp(-z^2*0.5)
    stop

    db_gsm_fg = db_gsm
    for ii=0,ndim-1 do db_gsm_fg[*,ii] = db_gsm[*,ii]-db_gsm_bg[*,ii]
    stop
    db_gsm_bg_var = prefix+'db_gsm_bg'
    store_data, db_gsm_bg_var, common_times, db_gsm_bg
    add_setting, db_gsm_bg_var, /smart, {$
        display_type: 'vector', $
        unit: 'nT', $
        short_name: 'dB', $
        coord: 'GSM', $
        coord_labels: xyz}

    stop





end

; Low frequ fluctations around perigee.
time_range = time_double(['2014-02-20','2014-02-23'])
probe = 'a'
time_range = time_double(['2013-06-06','2013-06-09'])
probe = 'a'
test_pflux_grant_preprocess_ebfield, probe=probe
end
