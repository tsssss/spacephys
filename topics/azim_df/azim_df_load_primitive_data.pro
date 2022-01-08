;+
; Load orbit data xxx_r_gsm, and magnetic field data xxx_b_gsm
; Save them to cdf, to avoid reloading again.
; Load data for long times: [2007-03-01,2009-09-01], and [2012-10-01,2017-10-01].
;
; time_range. The time range, does not count pad_time.
; probe. A string for mission_probe.
; time_step. A number of the common time step.
; data_types. A list of the requeste quantities.
;-

function azim_df_load_primitive_data_orbit, probe_info=probe_info, times=common_times
    ; Internal.

    routine_name = probe_info.routine_name+'_read_orbit'
    prefix = probe_info.prefix
    time_step = common_times[1]-common_times[0]
    full_time_range = minmax(common_times)
    r_gsm_var = prefix+'r_gsm'
    ndim = 3
    ncommon_time = n_elements(common_times)
    r_gsm = fltarr(ncommon_time,ndim)+!values.f_nan

    easy_missions = list(['themis','rbsp','goes'],/extract)
    if easy_missions.where(probe_info.routine_name) ne !null then begin
        call_procedure, routine_name, full_time_range, probe=probe_info.probe, coord='gsm'
        r_gsm = get_var_data(r_gsm_var, at=common_times)
    endif else begin
        secofday = constant('secofday')
        ndate = total(full_time_range*[-1,1])/secofday
        dates = smkarthm(full_time_range[0], secofday, ndate, 'x0')
        foreach date, dates do begin
            lprmsg, 'Processing '+time_string(date,tformat='YYYY-MM-DD')+' ...'
            time_range = date+[0,secofday]
            times = make_bins(time_range, time_step)
            del_data, r_gsm_var
            call_procedure, routine_name, time_range, probe=probe_info.probe, errmsg=errmsg
            if tnames(r_gsm_var) eq '' then begin
                lprmsg, 'No data, skip ...'
                continue
            endif
            index = (times-full_time_range[0])/time_step
            r_gsm[index,*] = get_var_data(r_gsm_var, at=times)
        endforeach
    endelse

    store_data, r_gsm_var, common_times, r_gsm
    return, r_gsm

end

function azim_df_load_primitive_data_bfield, probe_info=probe_info, times=common_times
    ; Internal.

    routine_name = probe_info.routine_name+'_read_bfield'
    prefix = probe_info.prefix
    time_step = common_times[1]-common_times[0]
    full_time_range = minmax(common_times)
    b_gsm_var = prefix+'b_gsm'
    ndim = 3
    ncommon_time = n_elements(common_times)
    b_gsm = fltarr(ncommon_time,ndim)+!values.f_nan

    secofday = constant('secofday')
    ndate = total(full_time_range*[-1,1])/secofday
    dates = smkarthm(full_time_range[0], secofday, ndate, 'x0')
    foreach date, dates do begin
        lprmsg, 'Processing '+time_string(date,tformat='YYYY-MM-DD')+' ...'
        time_range = date+[0,secofday]
        times = make_bins(time_range, time_step)
        del_data, b_gsm_var
        call_procedure, routine_name, time_range, probe=probe_info.probe, errmsg=errmsg
        if tnames(b_gsm_var) eq '' then begin
            lprmsg, 'No data, skip ...'
            continue
        endif
        index = (times-full_time_range[0])/time_step
        b_gsm[index,*] = get_var_data(b_gsm_var, at=times)
    endforeach

    b_sm_var = prefix+'b_sm'
    b_sm = cotran(b_gsm, common_times, 'gsm2sm')
    store_data, b_sm_var, common_times, b_sm
    return, b_sm

end


pro azim_df_load_primitive_data, time_range=time_range, probe=probe, project=project, data_file=data_file, reset=reset

    if n_elements(time_range) ne 2 then begin
        errmsg = handle_error('No input time range ...')
        return
    endif
    if n_elements(probe) eq 0 then begin
        errmsg = handle_error('No input probe ...')
        return
    endif
    probe_info = resolve_probe(probe)
    prefix = probe_info.prefix
    the_probe = probe_info.probe
    xyz = constant('xyz')

    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(data_file) eq 0 then data_file = join_path([project.data_dir,'azim_df_primitive_data_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hh'),'_to_')+'_v01.cdf'])
    if keyword_set(reset) then file_delete, data_file, /allow_nonexistent


;---Prepare event_info.
    info_var = 'event_info'
    if file_test(data_file) eq 0 then begin
        event_info = dictionary()
        event_info['time_range'] = time_range
        event_info['file_suffix'] = fgetbase(data_file)
        event_info['orbit_time_step'] = project.orbit_time_step
        event_info['bfield_time_step'] = project.time_step
        cdf_save_var, info_var, value=0, filename=data_file
        cdf_save_setting, event_info, varname=info_var, filename=data_file
    endif
    event_info = cdf_read_setting(info_var, filename=data_file)
    event_info_updated = 0

;---Orbit data.
    orbit_time_var = 'orbit_ut'
    orbit_time_step = event_info.orbit_time_step
    if ~cdf_has_var(orbit_time_var, filename=data_file) then begin
        orbit_times = make_bins(time_range, orbit_time_step)
        cdf_save_var, orbit_time_var, value=orbit_times, filename=data_file
        setting = dictionary($
            'description', 'ut sec of orbit data', $
            'unit', 'sec')
        cdf_save_setting, setting, varname=orbit_time_var, filename=data_file
    endif

    r_gsm_var = prefix+'r_gsm'
    if ~cdf_has_var(r_gsm_var, filename=data_file) then begin
        orbit_times = cdf_read_var(orbit_time_var, filename=data_file)
        r_gsm = azim_df_load_primitive_data_orbit(probe_info=probe_info, times=orbit_times)
        cdf_save_var, r_gsm_var, value=r_gsm, filename=data_file
        setting = dictionary($
            'display_type', 'vector', $
            'short_name', 'R', $
            'unit', 'Re', $
            'coord', 'GSM', $
            'coord_labels', xyz, $
            'depend_0', orbit_time_var)
        cdf_save_setting, setting, varname=r_gsm_var, filename=data_file
    endif

;---B field data.
    bfield_time_var = 'bfield_ut'
    bfield_time_step = event_info.bfield_time_step
    if ~cdf_has_var(bfield_time_var, filename=data_file) then begin
        bfield_times = make_bins(time_range, bfield_time_step)
        cdf_save_var, bfield_time_var, value=bfield_times, filename=data_file
        setting = dictionary($
            'description', 'ut sec of bfield data', $
            'unit', 'sec')
        cdf_save_setting, setting, varname=bfield_time_var, filename=data_file
    endif

    b_sm_setting = dictionary($
            'display_type', 'vector', $
            'short_name', 'B', $
            'unit', 'nT', $
            'coord', 'SM', $
            'coord_labels', xyz, $
            'depend_0', bfield_time_var)
    b_sm_var = prefix+'b_sm'
    b_gsm_var = prefix+'b_gsm'
    
    ; See if B GSM exist, if so then convert it to SM.
    if cdf_has_var(b_gsm_var, filename=data_file) then begin
        cdf_load_var, b_gsm_var, filename=data_file
        get_data, b_gsm_var, times, b_gsm
        b_sm = float(cotran(b_gsm, times, 'gsm2sm'))
        store_data, b_sm_var, times, b_sm
        cdf_save_var, b_sm_var, value=b_sm, filename=data_file
        setting = dictionary($
            'display_type', 'vector', $
            'short_name', 'B', $
            'unit', 'nT', $
            'coord', 'SM', $
            'coord_labels', xyz, $
            'depend_0', bfield_time_var)
        cdf_save_setting, b_sm_setting, varname=b_sm_var, filename=data_file
        cdf_del_var, b_gsm_var, filename=data_file
    endif
    
    if ~cdf_has_var(b_sm_var, filename=data_file) then begin
        bfield_times = cdf_read_var(bfield_time_var, filename=data_file)
        b_sm = azim_df_load_primitive_data_bfield(probe_info=probe_info, times=bfield_times)
        cdf_save_var, b_sm_var, value=b_sm, filename=data_file
        cdf_save_setting, b_sm_setting, varname=b_gsm_var, filename=data_file
    endif


    if event_info_updated then cdf_save_setting, event_info, varname=info_var, filename=data_file

end


    project = azim_df_load_project()
    search_settings = project.search_settings
    foreach search, search_settings do begin
        time_range = search.time_range
        probes = search.probes
        data_file = join_path([project.data_dir,search.data_file_suffix])
        foreach probe, probes do azim_df_load_primitive_data, $
            time_range=time_range, probe=probe, project=project, data_file=data_file
    endforeach

end
