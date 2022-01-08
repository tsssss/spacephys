;+
; Proprocess E field for a given time range.
; Load data by day by checking buffer.
; Process data by orbit.
; Save data by day.
;-

pro pflux_grant_preprocess_efield_per_orbit, orbit_time_range, probe=probe, project=project

    prefix = 'rbsp'+probe+'_'
    errmsg = ''
    fillval = !values.f_nan

    vsvy_var = prefix+'vsvy'
    vsvy = get_var_data(vsvy_var, in=orbit_time_range, times=efield_times)
    nefield_time = n_elements(efield_times)
    efield_time_step = efield_times[1]-efield_times[0]

    e_uv_var = prefix+'e_uv'
    e_colors = sgcolor(['red','blue'])
    e_labels = 'E'+['u','v']


;---Skip orbit during maneuver.
    maneuver_times = rbsp_read_maneuver_time(orbit_time_range, probe=probe)
    has_maneuver = 0
    foreach maneuver_time, maneuver_times do begin
        if product(maneuver_time-orbit_time_range) le 0 then has_maneuver = 1
    endforeach
    if has_maneuver then begin
        e_uv = fltarr(nefield_time,2)
        append_data, e_uv_var, efield_times, e_uv, limits={$
            ytitle: '(mV)', $
            colors: e_colors, $
            labels: e_labels}
        return
    endif

;---Calculate Eu and Ev.
    boom_lengths = [100d,100,12]
    spin_period = 12d
    smooth_width = spin_period/efield_time_step
    effective_boom_lengths = boom_lengths
    eu = (vsvy[*,0]-vsvy[*,1])/effective_boom_lengths[0]*1e3   ; V -> V/m -> mV/m.
    ev = (vsvy[*,2]-vsvy[*,3])/effective_boom_lengths[1]*1e3
    e_uv = [[eu],[ev]]
    for ii=0,1 do e_uv[*,ii] -= smooth(e_uv[*,ii], smooth_width, /nan, /edge_truncate)


;---Load boom_flags, which has incorpreted SDT and eclipse.
    rbsp_efw_read_boom_flag, orbit_time_range, probe=probe, errmsg=errmsg
    if errmsg ne '' then return
    flag_time_step = 10.
    nboom = 4
    flag_var = prefix+'boom_flag'
    boom_flags = get_var_data(flag_var, in=orbit_time_range, times=flag_times)
    index = where(total(boom_flags,2) lt nboom, count)
    if count eq 0 then begin
        bad_time_ranges = !null
        nbad_time_range = 0.
    endif else begin
        bad_time_ranges = time_to_range(flag_times[index], time_step=flag_time_step)
        nbad_time_range = n_elements(bad_time_ranges)/2
    endelse


;---Remove bad data according to the flag.
    pad_time = 60.
    for ii=0, nbad_time_range-1 do begin
        current_time_range = reform(bad_time_ranges[ii,*])+[-1,1]*pad_time
        lprmsg, time_string(current_time_range)
        index = lazy_where(efield_times, '[]', current_time_range, count=count)
        if count eq 0 then continue
        e_uv[index,*] = fillval
    endfor


;---Remove background dE.
    emag = snorm(e_uv)
    width2 = 30d/efield_time_step
    emag0 = smooth(emag, width2, /nan, /edge_truncate)
;    nrec = nefield_time
;    emag0[0:width2] = emag[0:width2]
;    emag0[nrec-1-width2:nrec-1] = emag[nrec-1-width2:nrec-1]

    ; Weight to get the background dE.
    perigee_shell = 4.  ; Use data beyond 4 Re.
    shell_width = 600.  ; sec
    r_var = prefix+'r_gsm'
    dis = snorm(get_var_data(r_var, at=efield_times))
    perigee_shell_time = efield_times[minmax(where(dis ge perigee_shell))]
    shell_times = perigee_shell_time+[-1,1]*shell_width
    weight = fltarr(nefield_time)
    middle_index = nefield_time/2
    weight[0:middle_index-1] = (tanh((efield_times[0:middle_index-1]-shell_times[0])/shell_width)+1)*0.5
    weight[middle_index:-1] = (1-tanh((efield_times[middle_index:-1]-shell_times[1])/shell_width))*0.5
    ratio = (1-weight)*emag0/emag
    for ii=0,1 do begin
        data = e_uv[*,ii]
        e_uv[*,ii] = data-data*ratio
    endfor


;---Weight down data within perigee shell.
    perigee_shell_time = efield_times[minmax(where(dis ge perigee_shell))]
    middle_index = nefield_time*0.5
    weight2 = fltarr(nefield_time)
    decay_times = (orbit_time_range+perigee_shell_time)*0.5
    decay_width = min(abs(orbit_time_range-decay_times))*0.8
    weight2[0:middle_index-1] = (tanh((efield_times[0:middle_index-1]-decay_times[0])/decay_width)+1)*0.5
    weight2[middle_index:-1] = (1-tanh((efield_times[middle_index:-1]-decay_times[1])/decay_width))*0.5
    for ii=0,1 do e_uv[*,ii] *= weight2


;---Remove spikes within perigee shell.
    perigee_shell = 2.
    max_perigee_emag = 5.
    index = where(emag ge max_perigee_emag and dis lt perigee_shell, count)
    if count ne 0 then e_uv[index,*] = !values.f_nan

;---Add to buffer.
    append_data, e_uv_var, efield_times, e_uv, limits={$
        ytitle: '(mV)', $
        colors: e_colors, $
        labels: e_labels}

end



pro pflux_grant_preprocess_efield_save_data, date, probe=probe, project=project

    secofday = constant('secofday')
    if date mod secofday ne 0 then message, 'Invalid date ...'

    rbspx = 'rbsp'+probe
    prefix = rbspx+'_'
    local_root = join_path([default_local_root(),'sdata','rbsp'])
    paths = [rbspx,'ebfield']
    time_range = date+[0,secofday]
    base_name = rbspx+'_ebfield_'+time_string(date,tformat='YYYY_MMDD')+'_v01.cdf'
    year = time_string(date,tformat='YYYY')
    data_file = join_path([local_root,paths,year,base_name])

    lprmsg, 'Save data to '+data_file+' ...'
    if file_test(data_file) eq 0 then return

    time_var = 'ut_sec'
    common_times = cdf_read_var(time_var, filename=data_file)

    e_uv_var = prefix+'e_uv'
    e_uv = get_var_data(e_uv_var, times=efield_times)
    index = lazy_where(efield_times, '[)', time_range)
    e_uv = float(e_uv[index,*])
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'E', $
        'unit', 'mV/m', $
        'coord', 'UVW', $
        'coord_labels', ['u','v'])
    ;if cdf_has_var(e_uv_var) then cdf_del_var, e_uv_var, filename=data_file
    cdf_save_var, e_uv_var, value=e_uv, filename=data_file
    cdf_save_setting, settings, filename=data_file, varname=e_uv_var


    ; Clear buffer.
    lprmsg, 'Clear buffer before '+time_string(time_range[1])+' ...'
    buffer_vars = prefix+['vsvy','r_gsm','e_uv']
    buffer_var = prefix+'buffer_time_range'
    foreach var, buffer_vars do begin
        data = get_var_data(var, times=times)
        index = where(times ge time_range[1], count)
        if count eq 0 then begin
            del_data, var
        endif else begin
            store_data, var, times[index], data[index,*,*,*,*,*,*,*]
        endelse
    endforeach

    var = buffer_vars[0]
    if tnames(var) eq '' then begin
        buffer_time_range = date+secofday+[0,0]
    endif else begin
        buffer_time_range = get_var_data(buffer_var)
        buffer_time_range[0] -= secofday
    endelse
    store_data, buffer_var, 0, buffer_time_range

end

pro pflux_grant_preprocess_efield_load_data, date, probe=probe, project=project

    secofday = constant('secofday')
    ndim = 3
    model = 't89'
    par = 2.
    xyz = constant('xyz')
    rgb = constant('rgb')
    uvw = ['u','v','w']
    uv = uvw[0:1]
    prefix = 'rbsp'+probe+'_'

    efield_time_step = 1d/16
    orbit_time_step = 60.
    spin_period = 12d
    dc_offset_nspin = 1
    max_valid_e = 300.
    max_valid_v = 200.
    omega = (2*!dpi)/86400d  ;Earth's rotation angular frequency
    nboom = 4
    v_colors = sgcolor(['red','green','blue','black'])
    v_labels = 'V'+string(findgen(nboom)+1,format='(I0)')
    buffer_vars = prefix+['vsvy','r_gsm']
    buffer_var = prefix+'buffer_time_range'

    date_time_range = date+[0,secofday]

    ; Load vsvy.
    datatype = (date ge time_double('2016-02-28'))? 'l2%vsvy-highres2': 'l2%vsvy-highres'
    rbsp_read_efw, date_time_range, id=datatype, probe=probe
    original_vsvy_var = 'vsvy'
    get_data, original_vsvy_var, times, vsvy
    if n_elements(vsvy) lt 6 then begin
        msg = strupcase(probe)+ ': no V data on '+time_string(date,tformat='YYYY/MM-DD')
        errmsg = handle_error(msg)
        del_data, buffer_vars
        store_data, buffer_var, 0, date_time_range[1]+[0,secofday]
        return
    endif

    ; Update buffer time range.
    buffer_time_range = get_var_data(buffer_var)
    buffer_time_range = minmax([buffer_time_range,date_time_range])
    store_data, buffer_var, 0, buffer_time_range


    ; Preprocess vsvy for the current day.
    uniform_time, original_vsvy_var, efield_time_step
    get_data, original_vsvy_var, times, vsvy
    vsvy = vsvy[*,0:3]
    index = where(abs(vsvy) gt max_valid_v, count)
    if count ne 0 then vsvy[index] = !values.f_nan
    ntime = n_elements(times)
    ntime0 = secofday/efield_time_step
    if ntime ne ntime0 then begin
        new_times = make_bins(date_time_range, efield_time_step)
        new_vsvy = fltarr(ntime0,nboom)+!values.f_nan
        index = lazy_where(new_times, '[]', minmax(times), count=count)
        if count ne 0 then new_vsvy[index,*] = sinterpol(vsvy, times, new_times[index], /nan)
        times = new_times
        vsvy = new_vsvy
    endif

    ; Add to buffer.
    vsvy_var = prefix+'vsvy'
    append_data, vsvy_var, times, vsvy, limits={$
        ytitle: '(V)', $
        colors: v_colors, $
        labels: v_labels }


    ; Load R.
    r_var = prefix+'r_gsm'
    if tnames(r_var) ne '' then get_data, r_var, old_times, old_data else begin
        old_times = []
        old_data = []
    endelse
    rbsp_read_orbit, date_time_range, probe=probe
    get_data, r_var, times, data
    store_data, r_var, [old_times,times], [old_data,data]
    uniform_time, r_var, orbit_time_step

end

pro pflux_grant_preprocess_efield, probe=probe, project=project

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project.pflux_calc_setting

    if n_elements(probe) eq 0 then return
    if (where(probe eq ['a','b']))[0] eq -1 then return

    rbspx = 'rbsp'+probe
    prefix = 'rbsp'+probe+'_'
    secofday = constant('secofday')

    mission_time_range = pflux_calc_setting[rbspx].time_range
;; test.
mission_time_range = time_double(['2018-04-01','2018-04-02'])
    data_file_dates = make_bins(mission_time_range-[0,secofday], secofday)
    data_file_id = 0


    ; Decompose into orbits.
    orbit_var = prefix+'orbit_time_range'
    if tnames(orbit_var) eq '' then begin
        orbit_time_ranges = pflux_grant_calculate_orbit_time_range(mission_time_range+[-1,1]*secofday, probe=probe)
        index = where(orbit_time_ranges[*,0] le mission_time_range[0])
        orbit_time_ranges = orbit_time_ranges[index[-1]:*,*]
        index = where(orbit_time_ranges[*,1] ge mission_time_range[1])
        orbit_time_ranges = orbit_time_ranges[0:index[0],*]
        store_data, orbit_var, 0, orbit_time_ranges
    endif
    orbit_time_ranges = get_var_data(orbit_var)
    norbit = n_elements(orbit_time_ranges)/2
    orbit_period = mean(orbit_time_ranges[*,1]-orbit_time_ranges[*,0])



    buffer_var = prefix+'buffer_time_range'
    if tnames(buffer_var) eq '' then begin
        start_date = orbit_time_ranges[0,0]-(orbit_time_ranges[0,0] mod secofday)
        store_data, buffer_var, 0, start_date+[0,0]
    endif

    for orbit_id=0, norbit-1 do begin
    ;---See if need to load another day of data.
        buffer_time_range = get_var_data(buffer_var)
        buffer_size = total(buffer_time_range*[-1,1])
        orbit_time_range = reform(orbit_time_ranges[orbit_id,*])

        ; Load a day of data if buffer_size = 0.
        if buffer_size eq 0 then begin
            pflux_grant_preprocess_efield_load_data, buffer_time_range[0], probe=probe
        endif


        ; No data on that day, move to the next orbit.
        buffer_time_range = get_var_data(buffer_var)
        buffer_size = total(buffer_time_range*[-1,1])
        if buffer_size eq 0 then continue

        ; Current orbit is passed because no data is available.
        if orbit_time_range[0] lt buffer_time_range[0] then continue

        ; Current orbit spans 2 days.
        if orbit_time_range[1] ge buffer_time_range[1] then begin
            pflux_grant_preprocess_efield_load_data, buffer_time_range[1], probe=probe, project=project
        endif

        ; No data on the 2nd day, move to the next orbit.
        buffer_time_range = get_var_data(buffer_var)
        buffer_size = total(buffer_time_range*[-1,1])
        if buffer_size eq 0 then continue


    ;---Process data within the current orbit.
        pflux_grant_preprocess_efield_per_orbit, orbit_time_range, probe=probe, project=project


    ;---Check if we need to save data.
        while product(data_file_dates[data_file_id]-buffer_time_range) gt 0 do data_file_id += 1
        data_file_date = data_file_dates[data_file_id]
        if (orbit_time_range[1]-data_file_date) ge secofday then begin
            pflux_grant_preprocess_efield_save_data, data_file_date, probe=probe, project=project
            data_file_id += 1
        endif
    endfor

end


foreach probe, ['b'] do pflux_grant_preprocess_efield, probe=probe, project=project
end
