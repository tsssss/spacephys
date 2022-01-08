;+
; Load data of r_gsm, b_gsm. Calculate dis, MLT, MLat, the detrended tilt angle.
; Downsample to the common data rate.
;
; Should be used internally, or for a general event.
;-

pro azim_df_load_data_per_event, time_range=time_range, reload=reload, filename=file, mission_probes=mission_probes, time_step=time_step, $
    available_probes=available_probes

    if n_elements(time_range) ne 2 then begin
        errmsg = handle_error('Invalid input time range ...')
        return
    endif
    if n_elements(file) eq 0 then begin
        event_id = time_string(mean(time_range),tformat='YYYY_MMDD_hh')
        file = join_path([homedir(),event_id+'_data.tplot'])
    endif
    ; Delete file if reload is wanted.
    if keyword_set(reload) then file_delete, file, /allow_nonexistent
    ; Do not reload data if it exists.
    if file_test(file) then begin
        tplot_restore, filename=file
        return
    endif


;---Load data.
    deg = 180d/!dpi
    ;if n_elements(mission_probes) eq 0 then mission_probes = ['rbsp'+letters('b'),'th'+letters('e'),'g'+['13','14','15'],'mms1','c1']
    if n_elements(mission_probes) eq 0 then mission_probes = ['rbsp'+letters('b'),'th'+letters('e'),'g'+['13','14','15'],'mms1']
    if n_elements(time_step) eq 0 then time_step = 10.
    del_data, '*'
    flags = list()  ; 1 for have all data.
    common_times = make_bins(time_range, time_step)
    foreach mission_probe, mission_probes do begin
        the_flag = 1
        mission_info = resolve_probe(mission_probe)
        routine_name = mission_info.routine_name
        probe = mission_info.probe
        prefix = mission_info.prefix

    ;---Load xxx_r_gsm.
        bad_data = 0
        r_gsm_var = prefix+'r_gsm'
        call_procedure, routine_name+'_read_orbit', time_range, probe=probe, errmsg=errmsg
        if tnames(r_gsm_var) eq '' then bad_data = 1 else begin
            get_data, r_gsm_var, times
            if n_elements(times) eq 1 then bad_data = 1
        endelse
        if errmsg ne '' then bad_data = 1
        if bad_data then begin
            the_flag = 0
            flags.add, the_flag
            continue
        endif

    ;---Load xxx_b_gsm.
        bad_data = 0
        b_gsm_var = prefix+'b_gsm'
        call_procedure, routine_name+'_read_bfield', time_range, probe=probe, errmsg=errmsg
        if tnames(b_gsm_var) eq '' then bad_data = 1 else begin
            get_data, b_gsm_var, times
            if n_elements(times) eq 1 then bad_data = 1
        endelse
        if errmsg ne '' then bad_data = 1
        if bad_data then begin
            the_flag = 0
            flags.add, the_flag
            continue
        endif
        
    ;---Interpolate to the common_times.
        flags.add, the_flag
        foreach var, [r_gsm_var,b_gsm_var] do interp_time, var, common_times

    ;---Calculate xxx_r_sm.
        r_sm_var = prefix+'r_sm'
        get_data, r_gsm_var, times, r_gsm
        r_sm = cotran(r_gsm, times, 'gsm2sm')
        store_data, r_sm_var, times, r_sm
        add_setting, r_sm_var, /smart, {$
            display_type: 'vector', $
            unit: 'Re', $
            short_name: 'R', $
            coord: 'SM', $
            coord_labels: ['x','y','z']}
            
    ;---Calculate xxx_b_sm.
        b_sm_var = prefix+'b_sm'
        get_data, b_gsm_var, times, b_gsm
        b_sm = cotran(b_gsm, times, 'gsm2sm')
        store_data, b_sm_var, times, b_sm
        add_setting, b_sm_var, /smart, {$
            display_type: 'vector', $
            unit: 'nT', $
            short_name: 'B', $
            coord: 'SM', $
            coord_labels: ['x','y','z']}

    ;---Calculate xxx_mlat, xxx_mlt.
        r_mag = cotran(r_gsm, times, 'gsm2mag')
        mlat = asin(r_mag[*,2]/snorm(r_mag))*deg
        mlon = atan(r_mag[*,1],r_mag[*,0])*deg
        mlt = mlon2mlt(mlon, times)

        mlat_var = prefix+'mlat'
        ;vars.add, mlat_var
        store_data, mlat_var, times, mlat
        add_setting, mlat_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'MLat'}

        mlt_var = prefix+'mlt'
        store_data, mlt_var, times, mlt
        add_setting, mlt_var, /smart, {$
            display_type: 'scalar', $
            unit: 'hr', $
            range: [-12,12], $
            short_name: 'MLT'}

    ;---Calculate the model bfield.
        par = 2.
        ndim = 3
        ntime = n_elements(times)
        b0gsm = fltarr(ntime, ndim)
        for ii=0, ntime-1 do begin
            tilt = geopack_recalc(times[ii])
            ; in-situ position
            rx = r_gsm[ii,0]
            ry = r_gsm[ii,1]
            rz = r_gsm[ii,2]
            ; in-situ B field.
            geopack_igrf_gsm, rx,ry,rz, bx,by,bz
            geopack_t89, par, rx,ry,rz, dbx,dby,dbz
            b0gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endfor
        bmod_var = prefix+'bmod_gsm'
        store_data, bmod_var, times, b0gsm
        add_setting, bmod_var, /smart, {$
            display_type: 'vector', $
            unit: 'nT', $
            short_name: 'B!S!UT89!N!S', $
            coord: 'GSM', $
            coord_labels: ['x','y','z']}

    ;---Calculate tilt angle.
        vars = prefix+['b','bmod']
        foreach var, vars do begin
            get_data, var+'_gsm', times, bgsm
            bsm = cotran(bgsm, times, 'gsm2sm')
            tilt = atan(bsm[*,2],snorm(bsm[*,0:1]))*deg
            the_var = var+'_tilt'
            store_data, the_var, times, tilt
            add_setting, the_var, /smart, {$
                display_type: 'scalar', $
                unit: 'deg', $
                short_name: 'Tilt'}
        endforeach

    ;---Calculate the tilt angle with model subtracted.
        tilt_var = prefix+'db_tilt'
        sys_subtract, prefix+'b_tilt', prefix+'bmod_tilt', to=tilt_var
        add_setting, tilt_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'Tilt'}
    endforeach


;---Make sure every probe have all data.
    flags = flags.toarray()
    index = where(flags eq 1, count)
    if count eq 0 then begin
        errmsg = handle_error('No data to save ...')
        return
    endif
    
    ; Gather all prefixes.
    var_suffixes = ['r_sm','mlt','mlat','db_tilt','b_sm']
    prefixes = list()
    foreach var_suffix, var_suffixes do begin
        vars = tnames('*'+var_suffix)
        foreach var, vars do begin
            the_prefix = get_prefix(var)
            if prefixes.where(the_prefix) ne !null then continue
            prefixes.add, the_prefix
        endforeach
    endforeach

    ; For each prefix, check all data exist, otherwise delete all data for the current prefix.
    foreach prefix, prefixes do begin
        vars = tnames(prefix+'*')
        foreach var, vars do begin
            del_var = 0
            if tnames(var) eq '' then del_var = 1 else begin
                get_data, var, times, data
                ntime = n_elements(times)
                ndata = (size(data,/dimensions))[0]
                del_var = (ntime eq ndata)? 0: 1
            endelse
            if del_var then begin
                del_data, vars
                break
            endif
        endforeach
    endforeach

    prefixes = tnames('*_'+var_suffixes[0])
    available_probes = prefixes
    foreach prefix, prefixes, ii do begin
        prefixes[ii] = get_prefix(prefixes[ii])
        available_probes[ii] = strmid(prefixes[ii],0,strpos(prefixes[ii],'_'))
    endforeach
    all_vars = ['dst','ae']
    omni_read_index, time_range
    foreach prefix, prefixes do all_vars = [all_vars, prefix+var_suffixes]
    
    common_times = make_bins(time_range, time_step)
    foreach var, vars do interp_time, var, common_times
    tplot_save, all_vars, filename=file

end
