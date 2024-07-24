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

    e_uv_var = prefix+'e_uv'
    vsvy_var = prefix+'vsvy'
    vsvy = get_var_data(vsvy_var, in=orbit_time_range, times=efield_times)
    nefield_time = n_elements(efield_times)


;---Skip orbit during maneuver.
    maneuver_times = rbsp_read_maneuver_time(orbit_time_range, probe=probe)
    has_maneuver = 0
    foreach maneuver_time, maneuver_times do begin
        if product(maneuver_time-orbit_time_range) le 0 then has_maneuver = 1
    endforeach
    if has_maneuver then begin
        e_uv = fltarr(nefield_time,2)
        append_data, e_uv_var, efield_times, e_uv
        return
    endif

;---Calculate Eu and Ev.
    boom_lengths = [100d,100,12]
    e_colors = sgcolor(['red','blue'])
    e_labels = 'E'+['u','v']
    effective_boom_lengths = boom_lengths
    eu = (vsvy[*,0]-vsvy[*,1])/effective_boom_lengths[0]*1e3   ; V -> V/m -> mV/m.
    ev = (vsvy[*,2]-vsvy[*,3])/effective_boom_lengths[1]*1e3
    e_uv = [[eu],[ev]]
    e_uv_var = prefix+'e_uv'
    append_data, e_uv_var, efield_times, e_uv, limits={$
        ytitle: '(mV)', $
        colors: e_colors, $
        labels: e_labels}


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
    e_uv = get_var_data(e_uv_var, times=efield_times)
    for ii=0, nbad_time_range-1 do begin
        current_time_range = reform(bad_time_ranges[ii,*])+[-1,1]*pad_time
        lprmsg, time_string(current_time_range)
        index = where_pro(efield_times, '[]', current_time_range, count=count)
        if count eq 0 then continue
        e_uv[index,*] = fillval
    endfor
    store_data, e_uv_var, efield_times, e_uv


;---The times within the orbit.
    time_index = where_pro(efield_times, '[]', orbit_time_range, count=ntime)
    times = efield_times[time_index]


;;---Subtract model E field.
;    r_var = prefix+'r_gsm'
;    ndim = 3
;    orbit_time_step = 60.
;    xyz = constant('xyz')
;    ; Get B using T89 model. The measured B have data gap and irregularities.
;    r_gsm = get_var_data(r_var, in=orbit_time_range, times=orbit_times)
;    norbit_time = n_elements(orbit_times)
;    v_gsm = fltarr(norbit_time,ndim)
;    for ii=0, ndim-1 do v_gsm[*,ii] = deriv(r_gsm[*,ii])
;    v_gsm *= (constant('re')/orbit_time_step)
;    v_var = prefix+'v_gsm'
;    store_data, v_var, orbit_times, v_gsm
;    add_setting, v_var, /smart, {$
;        display_type: 'vector', $
;        unit: 'km/s', $
;        short_name: 'V', $
;        coord: 'GSM', $
;        coord_labels: xyz }
;
;    bmod_gsm = fltarr(norbit_time,ndim)
;    par = 2.
;    model = 't89'
;    foreach time, orbit_times, ii do begin
;        tilt = geopack_recalc(time)
;        rx = r_gsm[ii,0]
;        ry = r_gsm[ii,1]
;        rz = r_gsm[ii,2]
;        ; in-situ B field.
;        geopack_igrf_gsm, rx,ry,rz, bx,by,bz
;        geopack_t89, par, rx,ry,rz, dbx,dby,dbz
;        bmod_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
;    endforeach
;    bmod_gsm_var = prefix+'bmod_gsm'
;    store_data, bmod_gsm_var, orbit_times, bmod_gsm
;    add_setting, bmod_gsm_var, /smart, {$
;        display_type: 'vector', $
;        unit: 'nT', $
;        short_name: 'B0!S!U'+strupcase(model)+'!N!R', $
;        coord: 'GSM', $
;        coord_labels: xyz}
;
;    ; E = -vxB.
;    evxb_gsm = scross(v_gsm,bmod_gsm)*1e-3   ; convert to mV/m.
;    evxb_var = prefix+'evxb_gsm'
;    store_data, evxb_var, orbit_times, evxb_gsm
;    add_setting, evxb_var, /smart, {$
;        display_type: 'vector', $
;        unit: 'mV/m', $
;        short_name: 'VxB E', $
;        coord: 'GSM', $
;        coord_labels: xyz }
;
;    ; E co-rotation
;    omega = (2*!dpi)/86400d  ;Earth's rotation angular frequency
;    r_gei = cotran(r_gsm, orbit_times, 'gsm2gei')
;    vcoro_gei = fltarr(norbit_time,ndim)
;    vcoro_gei[*,0] = -r_gei[*,1]*omega
;    vcoro_gei[*,1] =  r_gei[*,0]*omega
;    vcoro_gei[*,2] = 0.0
;    vcoro_gsm = cotran(vcoro_gei, orbit_times, 'gei2gsm')
;    ecoro_gsm = scross(vcoro_gsm, bmod_gsm)
;    ecoro_var = prefix+'ecoro_gsm'
;    store_data, ecoro_var, orbit_times, ecoro_gsm
;    add_setting, ecoro_var, /smart, {$
;        display_type: 'vector', $
;        unit: 'mV/m', $
;        short_name: 'Coro E', $
;        coord: 'GSM', $
;        coord_labels: xyz }
;
;    ; E model = E_vxb + E_corr
;    emod_var = prefix+'e0_gsm'
;    emod_gsm = evxb_gsm+ecoro_gsm
;    store_data, emod_var, orbit_times, emod_gsm
;    add_setting, emod_var, /smart, {$
;        display_type: 'vector', $
;        unit: 'mV/m', $
;        short_name: 'Model E', $
;        coord: 'GSM', $
;        coord_labels: xyz }
;    interp_time, emod_var, times
;
;    ; Get model E in UV.
;    e0_uv_var = prefix+'e0_uv'
;    the_var = prefix+'e0'
;    rbsp_gsm2uvw, the_var+'_gsm', the_var+'_uvw', quaternion=q_var, probe=probe
;    e0_uv = get_var_data(the_var+'_uvw', at=times)
;    e0_uv = e0_uv[*,0:1]
;    store_data, e0_uv_var, times, e0_uv
;    add_setting, e0_uv_var, /smart, {$
;        display_type: 'vector', $
;        unit: 'mV/m', $
;        short_name: 'E0', $
;        coord: 'UVW', $
;        coord_labels: e_labels, $
;        colors: e_colors }
;
;    e_uv[time_index,*] -= e0_uv

    spin_period = 12d
    efield_time_step = efield_times[1]-efield_times[0]
    smooth_width = spin_period/efield_time_step
    for ii=0,1 do begin
        dat = e_uv[*,ii]
        dat = smooth(dat, smooth_width, /nan, /edge_truncate)
        e_uv[*,ii] -= dat
    endfor

    emag = snorm(e_uv[time_index,*])
    nrec = n_elements(emag)
    width2 = 30d/efield_time_step
    emag0 = smooth(emag, width2, /nan, /edge_truncate)
    emag0[0:width2] = emag[0:width2]
    emag0[nrec-1-width2:nrec-1] = emag[nrec-1-width2:nrec-1]

    ; Weight to get the background dE.
    perigee_shell = 4.  ; Use data beyond 4 Re.
    shell_width = 600.  ; sec
    r_var = prefix+'r_gsm'
    dis = snorm(get_var_data(r_var, at=times))
    perigee_shell_time = times[minmax(where(dis ge perigee_shell))]
    shell_times = perigee_shell_time+[-1,1]*shell_width
    ntime = n_elements(times)
    weight = fltarr(ntime)
    middle_index = ntime/2
    weight[0:middle_index-1] = (tanh((times[0:middle_index-1]-shell_times[0])/shell_width)+1)*0.5
    weight[middle_index:-1] = (1-tanh((times[middle_index:-1]-shell_times[1])/shell_width))*0.5
    ratio = (1-weight)*emag0/emag
    for ii=0,1 do begin
        data = e_uv[time_index,ii]
        e_uv[time_index,ii] = data-data*ratio
    endfor
    
    
    ; Remove spikes within perigee shell.
    the_e_uv = e_uv[time_index,*]
    emag = snorm(the_e_uv)
    perigee_index = where(dis lt perigee_shell)
    perigee_emag = emag[perigee_index,*]
    perigee_e_uv = the_e_uv[perigee_index,*]
    max_perigee_emag = 10.
    index = where(perigee_emag ge max_perigee_emag, count, complement=good_index)
    if count ne 0 then begin
        perigee_e_uv[index,*] = !values.f_nan
        the_e_uv[perigee_index,*] = perigee_e_uv
        e_uv[time_index,*] = the_e_uv
    endif


;---Weight down data within perigee shell.
    dis = snorm(get_var_data(r_var, at=times))
    perigee_shell_time = times[minmax(where(dis ge perigee_shell))]
    middle_index = ntime*0.5
    weight2 = fltarr(ntime)
    decay_times = (orbit_time_range+perigee_shell_time)*0.5
    decay_width = min(abs(orbit_time_range-decay_times))*0.8
    weight2[0:middle_index-1] = (tanh((times[0:middle_index-1]-decay_times[0])/decay_width)+1)*0.5
    weight2[middle_index:-1] = (1-tanh((times[middle_index:-1]-decay_times[1])/decay_width))*0.5

    for ii=0,1 do e_uv[time_index,ii] *= weight2
    store_data, e_uv_var, efield_times, e_uv

end



pro pflux_grant_preprocess_efield_save_data, date, probe=probe, project=project

    secofday = constant('secofday')
    if date mod secofday ne 0 then message, 'Invalid date ...'

    rbspx = 'rbsp'+probe
    prefix = rbspx+'_'
    local_root = join_path([default_local_root(),'sdata','rbsp'])
    paths = [rbspx,'ebfield']
    paths = [rbspx,'efield']
    time_range = date+[0,secofday]
    base_name = rbspx+'_ebfield_'+time_string(date,tformat='YYYY_MMDD')+'_v01.cdf'
    year = time_string(date,tformat='YYYY')
    data_file = join_path([local_root,paths,year,base_name])
    
    
    ; Remove spikes within perigee shell.
    e_uv_var = prefix+'e_uv'
    e_uv = get_var_data(e_uv_var, times=efield_times)
    time_index = where_pro(efield_times, '[]', time_range)
    the_e_uv = e_uv[time_index,*]
    times = efield_times[time_index]
    emag = snorm(the_e_uv)

    r_var = prefix+'r_gsm'
    dis = snorm(get_var_data(r_var, at=times))

    perigee_shell = 2.    
    max_perigee_emag = 5.
    index = where(emag ge max_perigee_emag and dis lt perigee_shell, count)
    if count ne 0 then begin
        the_e_uv[index,*] = !values.f_nan
        e_uv[time_index,*] = the_e_uv
        store_data, e_uv_var, efield_times, e_uv
    endif

    lprmsg, 'Save data to '+data_file+' ...'

    ; For survey purpose.
    test = 0
    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    base_name = rbspx+'_test_preprocess_efield_'+time_string(date,tformat='YYYY_MMDD')+'_v01.pdf'
    plot_file = join_path([project.plot_dir,'test_preprocess_efield',rbspx,year,base_name])
    if keyword_set(test) then plot_file = test
    tplot_options, 'labflag', -1
    sgopen, plot_file, xsize=8, ysize=5
    
    rbsp_efw_read_boom_flag, time_range, probe=probe
    vsc_var = prefix+'vsc_median'
    get_data, vsc_var, times, vsc
    options, vsc_var, 'ytitle', 'Vsc (V)'


    r_var = prefix+'r_gsm'
    dis = snorm(get_var_data(r_var, at=times))
    index = where(dis le 2)
    vsc[index] = !values.f_nan
    store_data, vsc_var, times, vsc
    
    lshell = 4.
    orbit_time_step = 60
    dis = snorm(get_var_data(r_var, times=orbit_times))
    index = where(dis le lshell)
    perigee_times = time_to_range(orbit_times[index], time_step=orbit_time_step)
    
    vars = prefix+['e_uv','vsc_median']
    tplot, vars, trange=time_range
    timebar, perigee_times, color=sgcolor('black')
    
    maneuver_times = rbsp_read_maneuver_time(time_range, probe=probe)
    if n_elements(maneuver_times) ne 0 then timebar, maneuver_times, color=sgcolor('red'), linestyle=1
    
    if keyword_set(test) then stop
    sgclose



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
    quaternion_time_step = 1d/4
    spin_period = 12d
    dc_offset_nspin = 1
    max_valid_e = 300.
    max_valid_v = 200.
    omega = (2*!dpi)/86400d  ;Earth's rotation angular frequency
    nboom = 4
    v_colors = sgcolor(['red','green','blue','black'])
    v_labels = 'V'+string(findgen(nboom)+1,format='(I0)')
    buffer_vars = prefix+['vsvy','r_gsm','q_uvw2gsm']
    buffer_var = prefix+'buffer_time_range'

    date_time_range = date+[0,secofday]

    ; Load vsvy.
    vsvy_var = prefix+'vsvy'
    if tnames(vsvy_var) ne '' then get_data, vsvy_var, old_times, old_data else begin
        old_times = []
        old_data = []
    endelse
    datatype = (date ge time_double('2016-02-28'))? 'l2%vsvy-highres2': 'l2%vsvy-highres'
    rbsp_read_efw, date_time_range, id=datatype, probe=probe
    vsvy_var = rename_var('vsvy', output=vsvy_var)
    get_data, vsvy_var, times, vsvy
    if n_elements(vsvy) lt 6 then begin
        msg = strupcase(probe)+ ': no V data on '+time_string(date,tformat='YYYY/MM-DD')
        errmsg = handle_error(msg)
        del_data, buffer_vars
        store_data, buffer_var, 0, date_time_range[1]+[0,secofday]
        return
    endif else begin
        buffer_time_range = get_var_data(buffer_var)
        buffer_time_range = minmax([buffer_time_range,date_time_range])
        store_data, buffer_var, 0, buffer_time_range
    endelse

    vsvy = vsvy[*,0:3]
    index = where(abs(vsvy) gt max_valid_v, count)
    if count ne 0 then vsvy[index] = !values.f_nan
    store_data, vsvy_var, [old_times,times], [old_data,vsvy], limits={$
        ytitle: '(V)', $
        colors: v_colors, $
        labels: v_labels, $
        yrange: [0,1]+[-1,1]*0.2, $
        ytickv: [0,1], $
        yticks: 1, $
        yminor: 0}
    uniform_time, vsvy_var, efield_time_step

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

;    ; Load quaternion.
;    q_var = prefix+'q_uvw2gsm'
;    if tnames(q_var) ne '' then get_data, q_var, old_times, old_data else begin
;        old_times = []
;        old_data = []
;    endelse
;    rbsp_read_quaternion, date_time_range, probe=probe
;    get_data, q_var, times, data
;    store_data, q_var, [old_times,times], [old_data,data]
;    uniform_time, q_var, quaternion_time_step



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

probe = 'a'
pflux_grant_preprocess_efield, probe=probe, project=project
end