;+
; Load data from low-level to high:
;   1. orbit data: r_gsm, r_sm, mlat, mlt, bmod_gsm
;   2. magnetic field data: b_gsm
;   3. ...
;
; time_range. The time range, does not count pad_time.
; probe. A string for mission_probe.
; time_step. A number of the common time step.
; data_types. A list of the requeste quantities.
;-

pro azim_df_calc_tilt, b_var, b_tilt_var

    deg = constant('deg')
    coord = get_setting(b_var, 'coord', exist)
    if ~exist then coord = ''

    get_data, b_var, times, bvec
    tilt = asin(bvec[*,2]/snorm(bvec))*deg
    the_var = b_tilt_var
    store_data, the_var, times, tilt
    add_setting, the_var, /smart, {$
        display_type: 'scalar', $
        unit: 'deg', $
        coord: coord, $
        short_name: 'Tilt'}
end


pro azim_df_load_data_bfield, time_range, probe, time_step, data_types=data_types
    ; Internal.

    probe_info = resolve_probe(probe)
    prefix = probe_info.prefix
    common_times = make_bins(time_range, time_step)
    if n_elements(data_types) eq 0 then data_types = list()

    ; xxx_b_gsm.
    bgsm_var = prefix+'b_gsm'
    if check_if_update(b_gsm_var, time_range) then begin
        call_procedure, probe_info.routine_name+'_read_bfield', time_range, probe=probe_info.probe
        interp_time, bgsm_var, common_times
    endif
    bgsm = get_var_data(bgsm_var)

    ; xxx_b_sm.
    xyz = constant('xyz')
    bsm_var = prefix+'b_sm'
    if data_types.where('b_sm') ne !null then begin
        bsm = cotran(bgsm, common_times, 'gsm2sm')
        store_data, bsm_var, common_times, bsm
        add_setting, bsm_var, /smart, {$
            display_type: 'vector', $
            short_name: 'B', $
            unit: 'nT', $
            coord: 'SM', $
            coord_labels: xyz}
    endif

    ; xxx_b_tilt.
    b_tilt_var = prefix+'b_sm_tilt'
    if data_types.where('b_sm_tilt') ne !null then begin
        if check_if_update(bsm_var) then azim_df_load_data_bfield, time_range, probe, time_step, data_types='b_sm'
        azim_df_calc_tilt, bsm_var, btilt_var
    endif

end

pro azim_df_load_data_orbit, time_range, probe, time_step, data_types=data_types
    ; Internal.

    probe_info = resolve_probe(probe)
    prefix = probe_info.prefix
    common_times = make_bins(time_range, time_step)
    if n_elements(data_types) eq 0 then data_types = list()

    ; xxx_r_gsm.
    rgsm_var = prefix+'r_gsm'
    if check_if_update(r_gsm_var, time_range) then begin
        call_procedure, probe_info.routine_name+'_read_orbit', time_range, probe=probe_info.probe
        interp_time, rgsm_var, common_times
    endif
    rgsm = get_var_data(rgsm_var)

    ; xxx_r_sm.
    xyz = constant('xyz')
    rsm_var = prefix+'r_sm'
    if data_types.where('r_sm') ne !null then begin
        rsm = cotran(rgsm, common_times, 'gsm2sm')
        store_data, rsm_var, common_times, rsm
        add_setting, rsm_var, /smart, {$
            display_type: 'vector', $
            short_name: 'R', $
            unit: 'Re', $
            coord: 'SM', $
            coord_labels: xyz}
    endif

    ; xxx_mlt.
    deg = constant('deg')
    mlt_var = prefix+'mlt'
    if data_types.where('mlt') ne !null then begin
        rmag = cotran(rgsm, common_times, 'gsm2mag')
        mlat = asin(rmag[*,2]/snorm(rmag))*deg
        mlon = atan(rmag[*,1],rmag[*,0])*deg
        mlt = mlon2mlt(mlon, times)

        store_data, mlt_var, times, mlt
        add_setting, mlt_var, /smart, {$
            display_type: 'scalar', $
            unit: 'hr', $
            range: [-12.,12], $
            short_name: 'MLT'}
    endif

    ; xxx_mlat.
    deg = constant('deg')
    mlat_var = prefix+'mlat'
    if data_types.where('mlat') ne !null then begin
        rmag = cotran(rgsm, common_times, 'gsm2mag')
        mlat = asin(rmag[*,2]/snorm(rmag))*deg
        mlon = atan(rmag[*,1],rmag[*,0])*deg
        mlt = mlon2mlt(mlon, times)

        store_data, mlat_var, times, mlat
        add_setting, mlat_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'MLat'}
    endif

    ; xxx_bmod_gsm.
    bmod_gsm_var = prefix+'bmod_gsm'
    if data_types.where('bmod_gsm') ne !null then begin
        par = 2.
        ndim = 3
        ncommon_time = n_elements(common_times)
        bmod_gsm = fltarr(ncommon_time,ndim)
        foreach time, common_times, ii do begin
            tilt = geopack_recalc(times[ii])
            ; in-situ position
            rx = rgsm[ii,0]
            ry = rgsm[ii,1]
            rz = rgsm[ii,2]
            ; in-situ B field.
            geopack_igrf_gsm, rx,ry,rz, bx,by,bz
            geopack_t89, par, rx,ry,rz, dbx,dby,dbz
            bmod_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endforeach
        store_data, bmod_gsm_var, times, bmod_gsm
        add_setting, bmod_gsm_var, /smart, {$
            display_type: 'vector', $
            unit: 'nT', $
            short_name: 'B!S!U'+strupcase(model)+'!N!S', $
            model: strupcase(model), $
            model_par: par, $
            coord: 'GSM', $
            coord_labels: xyz}
    endif

    ; xxx_bmod_sm.
    bmod_sm_var = prefix+'bmod_sm'
    if data_types.where('bmod_sm') ne !null then begin
        if check_if_update(bmod_gsm_var) then azim_df_load_data_orbit, time_range, probe, time_step, data_types='bmod_gsm'
        bmod_gsm = get_var_data(bmod_gsm_var)
        bmod_sm = cotran(bmod_gsm, common_times, 'gsm2sm')

        store_data, bmod_sm_var, times, bmod_gsm
        add_setting, bmod_sm_var, /smart, {$
            display_type: 'vector', $
            unit: 'nT', $
            short_name: 'B!S!U'+strupcase(model)+'!N!S', $
            model: strupcase(model), $
            model_par: par, $
            coord: 'SM', $
            coord_labels: xyz}
    endif

    ; xxx_bmod_tilt.
    bmod_tilt_var = prefix+'bmod_sm_tilt'
    if data_types.where('bmod_sm_tilt') ne !null then begin
        if check_if_update(bmod_sm_var) then azim_df_load_data_orbit, time_range, probe, time_step, data_types='bmod_sm'
        azim_df_calc_tilt, bmod_sm_var, bmod_tilt_var
    endif

end


pro azim_df_load_data, time_range, probes=probes, project=project, data_file=data_file, reset=reset

    if n_elements(time_range) ne 2 then begin
        errmsg = handle_error('No input time range ...')
        return
    endif

    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(probes) eq 0 then probes = project.all_probes

    if n_elements(data_file) eq 0 then data_file = join_path([project.data_dir,'azim_df_data_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hh'),'_to_')+'_v01.tplot'])

    if keyword_set(reset) then file_delete, data_file, /allow_nonexistent
    del_data, '*'


;---Init file with the most basic event_info.
    if file_test(data_file) eq 0 then begin
        event_info = dictionary()
        event_info['all_probes'] = probes
        event_info['time_range'] = time_range
        event_info['file_suffix'] = fgetbase(data_file)
        event_info['time_step'] = project.time_step
        store_data, 'event_info', 0, event_info
        tplot_save, 'event_info', filename=data_file
    endif else tplot_restore, filename=data_file


;---Prepare event_info with more settings.
    event_info = get_var_data('event_info')
    time_step = event_info.time_step

    the_key = 'smooth_window'
    if ~event_info.haskey(the_key) then event_info[the_key] = constant('secofhour')
    smooth_window = event_info[the_key]

    the_key = 'pad_time'
    if ~event_info.haskey(the_key) then event_info[the_key] = smooth_window
    pad_time = event_info[the_key]

    the_key = 'full_time_range'
    if ~event_info.haskey(the_key) then event_info[the_key] = time_range+[-1,1]*pad_time
    full_time_range = event_info[the_key]


;---Orbit data.
    foreach probe, probes do begin
        prefix = probe+'_'
        if check_if_update(prefix+'r_gsm', full_time_range) then azim_df_load_data_orbit, full_time_range, probes, time_step
    endforeach

end
