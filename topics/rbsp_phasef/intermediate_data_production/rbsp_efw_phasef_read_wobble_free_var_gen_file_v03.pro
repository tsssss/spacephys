;+
; Save the following data to CDF:
;   q_uvw2gse. Fixed for artificial wobble at the spin period, Saved at 1 sec cadence.
;   r_mgse. In Re, 10 sec.
;   v_mgse. In km/s, 10 sec.
;   e_mgse. In mV/m, 10 sec, include Ex.
;   b_mgse. In nT, 10 sec.
;
; b_mgse in v02 still have some DC offset (~50 nT) relative to the emfisis data around perigee.
; Here, in v03, b_mgse have 0 DC offset on average.
;-

pro rbsp_efw_phasef_read_wobble_free_var_gen_file_v03, time, probe=probe, filename=file, errmsg=errmsg

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


;---Constants and settings.
    secofday = 86400d
    deg = 180d/!dpi
    rad = !dpi/180
    errmsg = ''
    common_time_step = 10d      ; sec.
    xyz = constant('xyz')


    ; Derived settings.
    date = time[0]-(time[0] mod secofday)
    time_range = date+[0,secofday]
    long_time_range = time_range+[-1,1]*3600d
    prefix = 'rbsp'+probe+'_'
    common_times = make_bins(time_range+[1,0]*common_time_step, common_time_step)
    ncommon_time = n_elements(common_times)


;---Init file.
    out_dir = fgetpath(file)
    if file_test(out_dir,/directory) eq 0 then file_mkdir, out_dir
    data_file = file
    if file_test(data_file) eq 1 then file_delete, data_file  ; overwrite old files.
    cdf_id = cdf_create(data_file)
    ginfo = dictionary($
        'TITLE', 'RBSP EFW prepare quantities related to residue removal around perigee', $
        'TEXT', 'Generated by Sheng Tian at the University of Minnesota' )
    cdf_save_setting, ginfo, filename=cdf_id



;---q_uvw2gse.
    rbsp_read_q_uvw2gse, long_time_range, probe=probe


;---Common times.
    utname = 'time'
    data = common_times
    settings = dictionary($
        'fieldnam', 'time', $
        'unit', 'sec', $
        'var_type', 'support_data' )
    cdf_save_var, utname, value=data, filename=cdf_id, cdf_type='CDF_DOUBLE'
    cdf_save_setting, settings, var=utname, filename=cdf_id


;---r_mgse and v_mgse.
    rbsp_efw_phasef_read_r_mgse, time_range, probe=probe
    rbsp_efw_phasef_read_v_mgse, time_range, probe=probe
    vnames = prefix+['r_mgse','v_mgse']
    foreach vname, vnames do begin
        data = get_var_data(vname, in=time_range, times=orbit_times)
        data = sinterpol(data, orbit_times, common_times, /quadratic)
        store_data, vname, common_times, data
    endforeach

    vname = prefix+'r_mgse'
    data = get_var_data(vname)
    settings = dictionary($
        'display_type', 'vector', $
        'short_name', 'R', $
        'unit', 'Re', $
        'coord', 'MGSE', $
        'coord_labels', xyz, $
        'var_type', 'data', $
        'depend_0', utname )
    cdf_save_var, vname, value=data, filename=cdf_id
    cdf_save_setting, settings, var=vname, filename=cdf_id

    vname = prefix+'v_mgse'
    data = get_var_data(vname)
    settings = dictionary($
        'display_type', 'vector', $
        'short_name', 'V', $
        'unit', 'km/s', $
        'coord', 'MGSE', $
        'coord_labels', xyz, $
        'var_type', 'data', $
        'depend_0', utname )
    cdf_save_var, vname, value=data, filename=cdf_id
    cdf_save_setting, settings, var=vname, filename=cdf_id


;---b_mgse.
    rbsp_efw_phasef_read_b_mgse, time_range+[-1,1]*common_time_step, probe=probe
    get_data, prefix+'b_mgse', times, b_mgse
    time_step = sdatarate(times)
    b_time_step = 1d/16
    min_count = common_time_step*0.5/b_time_step
    ndim = 3
    b_mgse_bg = fltarr(ncommon_time,ndim)+!values.f_nan
    for time_id=0,ncommon_time-1 do begin
        sec_time = common_times[time_id]+[-1,1]*common_time_step*0.5
        time_index = where_pro(times,'[]',sec_time, count=count)
        if count lt min_count then continue
        for dim_id=0,ndim-1 do begin
            dat = b_mgse[time_index,dim_id]
            index = where(finite(dat), count)
            if count lt min_count then continue
            b_mgse_bg[time_id,dim_id] = median(dat)
        endfor
    endfor
    store_data, prefix+'b_mgse', common_times, b_mgse_bg


;    ; Test.
;    b_mgse = sinterpol(b_mgse, times, common_times)
;    for ii=0,ndim-1 do begin
;        store_data, prefix+'db'+xyz[ii], common_times, b_mgse[*,ii]-b_mgse_bg[*,ii]
;    endfor
;    store_data, prefix+'dbmag', common_times, snorm(b_mgse)-snorm(b_mgse_bg)
;    stop

    vname = prefix+'b_mgse'
    data = get_var_data(vname)
    settings = dictionary($
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', 'nT', $
        'coord', 'MGSE', $
        'coord_labels', xyz, $
        'var_type', 'data', $
        'depend_0', utname )
    cdf_save_var, vname, value=data, filename=cdf_id
    cdf_save_setting, settings, var=vname, filename=cdf_id


;---e_mgse.
    rbsp_efw_read_e_mgse, time_range+[-1,1]*common_time_step, probe=probe, keep_spin_axis=1
    get_data, prefix+'e_mgse', times, e_mgse
    time_step = sdatarate(times)
    e_time_step = 1d/32
    min_count = common_time_step*0.5/e_time_step
    ndim = 3
    e_mgse_bg = fltarr(ncommon_time,ndim)+!values.f_nan
    for time_id=0,ncommon_time-1 do begin
        sec_time = common_times[time_id]+[-1,1]*common_time_step*0.5
        time_index = where_pro(times,'[]',sec_time, count=count)
        if count lt min_count then continue
        for dim_id=0,ndim-1 do begin
            dat = e_mgse[time_index,dim_id]
            index = where(finite(dat), count)
            if count lt min_count then continue
            e_mgse_bg[time_id,dim_id] = median(dat)
        endfor
    endfor
    store_data, prefix+'e_mgse', common_times, e_mgse_bg

    vname = prefix+'e_mgse'
    data = get_var_data(vname)
    settings = dictionary($
        'display_type', 'vector', $
        'short_name', 'E', $
        'unit', 'mV/m', $
        'coord', 'MGSE', $
        'coord_labels', xyz, $
        'var_type', 'data', $
        'depend_0', utname )
    cdf_save_var, vname, value=data, filename=cdf_id
    cdf_save_setting, settings, var=vname, filename=cdf_id


    cdf_close, cdf_id

end


;time = time_double('2012-12-31')
;probe = 'a'
;time = time_double('2018-09-27')
;time = time_double('2012-12-25')
;probe = 'b'
;
;time = time_double('2012-09-30')
;probe = 'a'
;
;file = join_path([homedir(),'test.cdf'])
;if file_test(file) then file_delete, file
;rbsp_efw_phasef_read_wobble_free_var_gen_file_v03, time, probe=probe, filename=file
;

probes = ['b']
root_dir = join_path([rbsp_efw_phasef_local_root()])
secofday = constant('secofday')
foreach probe, probes do begin
    prefix = 'rbsp'+probe+'_'
    prefix = 'rbsp'+probe+'_'
    time_range = rbsp_efw_phasef_get_valid_range('wobble_free_var', probe=probe)
    days = make_bins(time_range+[0,-1]*secofday, secofday)
    foreach day, days do begin
        time = day

        if time lt time_double('2012') then continue
        if time ge time_double('2013') then continue
        file = join_path([root_dir,'efw_phasef','wobble_free_var_v03','rbsp'+probe,time_string(time[0],tformat='YYYY'),$
            prefix+'efw_wobble_free_var_'+time_string(time,tformat='YYYY_MMDD')+'_v03.cdf'])
        if file_test(file) eq 1 then continue
        rbsp_efw_phasef_read_wobble_free_var_gen_file_v03, time, probe=probe, filename=file
    endforeach
endforeach
stop

probes = ['b']
root_dir = join_path([rbsp_efw_phasef_local_root()])
secofday = constant('secofday')
foreach probe, probes do begin
    prefix = 'rbsp'+probe+'_'
    prefix = 'rbsp'+probe+'_'
    time_range = rbsp_efw_phasef_get_valid_range('wobble_free_var', probe=probe)
    days = make_bins(time_range+[0,-1]*secofday, secofday)
    foreach day, days do begin
        time = day

if time lt time_double('2012') then continue
if time ge time_double('2013') then continue
        file = join_path([root_dir,'efw_phasef','wobble_free_var_v03','rbsp'+probe,time_string(time[0],tformat='YYYY'),$
            prefix+'efw_wobble_free_var_'+time_string(time,tformat='YYYY_MMDD')+'_v03.cdf'])
if file_test(file) eq 0 then continue
        print, file

        common_time_step = 10.
        day_time_range = time+[0,86400d]
        common_times = make_bins(day_time_range+[1,0]*common_time_step, common_time_step)
        ncommon_time = n_elements(common_times)

        xyz = constant('xyz')
        utname = 'time'

        rbsp_efw_read_e_mgse, day_time_range+[-1,1]*common_time_step, probe=probe, keep_spin_axis=1
        get_data, prefix+'e_mgse', times, e_mgse
        time_step = sdatarate(times)
        e_time_step = 1d/32
        min_count = common_time_step*0.5/e_time_step
        ndim = 3
        e_mgse_bg = fltarr(ncommon_time,ndim)+!values.f_nan
        for time_id=0,ncommon_time-1 do begin
            sec_time = common_times[time_id]+[-1,1]*common_time_step*0.5
            time_index = where_pro(times,'[]',sec_time, count=count)
            if count lt min_count then continue
            for dim_id=0,ndim-1 do begin
                dat = e_mgse[time_index,dim_id]
                index = where(finite(dat), count)
                if count lt min_count then continue
                e_mgse_bg[time_id,dim_id] = median(dat)
            endfor
        endfor
        store_data, prefix+'e_mgse', common_times, e_mgse_bg

        vname = prefix+'e_mgse'
        data = get_var_data(vname)
        settings = dictionary($
            'display_type', 'vector', $
            'short_name', 'E', $
            'unit', 'mV/m', $
            'coord', 'MGSE', $
            'coord_labels', xyz, $
            'var_type', 'data', $
            'depend_0', utname )
        cdf_save_var, vname, value=data, filename=file
        cdf_save_setting, settings, var=vname, filename=file
    endforeach
endforeach

end
