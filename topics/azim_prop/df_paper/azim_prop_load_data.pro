;+
; Load data of r_gsm, b_gsm for each spacecraft. Calculate in-situ MLT, MLon/MLat.
; Calculate the tilt angle, bmod_gsm_t89, |B|-|B| T89.
;-

pro azim_prop_load_data_per_event, project, event_id=event_id

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    trace_dir = -1
    model = 't89'
    string_theta = '!9'+string(113b)+'!X'

    constant = project.constant
    deg = constant.deg
    rgb = constant.rgb

    event = events[event_id]
    pad_time = 30.  ; in min.
    time_range = event.time_range+[0,1]*pad_time*60
    probes = event.probes


    vars = list()
    foreach probe, probes do begin
        ;if probe eq 'mms1' then continue
        pre0 = probe+'_'
        label = strupcase(project[probe].short_name)

    ;---Load xxx_r_gsm.
        call_procedure, project[probe].routine_name+'_read_orbit', time_range, probe=project[probe].probe
        rgsm_var = pre0+'r_gsm'
        if tnames(rgsm_var) eq '' then stop
        vars.add, rgsm_var
        options, rgsm_var, 'colors', rgb

        ; Convert xxx_r_gsm to xxx_r_sm.
        get_data, rgsm_var, times, rgsm
        rsm = cotran(rgsm, times, 'gsm2sm')
        rsm_var = pre0+'r_sm'
        ;vars.add, rsm_var
        store_data, rsm_var, times, rsm
        add_setting, rsm_var, /smart, {$
            display_type: 'vector', $
            unit: 'Re', $
            short_name: 'R', $
            coord: 'SM', $
            coord_labels: ['x','y','z'], $
            colors: rgb}

        ; Convert xxx_r_gsm to xxx_r_mag and get MLT, MLon/MLat.
        get_data, rgsm_var, times, rgsm
        rmag = cotran(rgsm, times, 'gsm2mag')
        mlat = asin(rmag[*,2]/snorm(rmag))*deg
        mlon = atan(rmag[*,1],rmag[*,0])*deg
        mlt = mlon2mlt(mlon, times)

        mlat_var = pre0+'mlat'
        ;vars.add, mlat_var
        store_data, mlat_var, times, mlat
        add_setting, mlat_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'MLat'}

        mlon_var = pre0+'mlon'
        ;vars.add, mlon_var
        store_data, mlon_var, times, mlon
        add_setting, mlon_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'MLon'}

        mlt_var = pre0+'mlt'
        vars.add, mlt_var
        store_data, mlt_var, times, mlt
        add_setting, mlt_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'MLT'}

        ; Get model field.
        read_geopack_info, pre0+'r_gsm', model=model, direction=trace_dir
        bmod_var = pre0+'bmod_gsm_'+model
        vars.add, bmod_var
        options, bmod_var, 'colors', rgb


    ;---Load xxx_bgsm.
        call_procedure, project[probe].routine_name+'_read_bfield', time_range, probe=project[probe].probe
        bgsm_var = pre0+'b_gsm'
        if tnames(bgsm_var) eq '' then stop
        vars.add, bgsm_var
        options, bgsm_var, 'colors', rgb

        ; Calculate |B|-|B| model.
        get_data, bgsm_var, times, bgsm
        bmag = snorm(bgsm)
        bmodgsm = get_var_data(bmod_var, at=times)
        bmodmag = snorm(bmodgsm)
        dbmag = bmag-bmodmag
        dbmag_var = pre0+'db_mag'
        vars.add, dbmag_var
        store_data, dbmag_var, times, dbmag
        add_setting, dbmag_var, /smart, {$
            display_type: 'scalar', $
            unit: 'nT', $
            short_name: '|B|-|B!D'+strupcase(model)+'!N|'}

        ; Calculate B tilt in SM.
        get_data, bgsm_var, times, bgsm
        bsm = cotran(bgsm, times, 'gsm2sm')
        btilt = atan(bsm[*,2],sqrt(bsm[*,0]^2+bsm[*,1]^2))*deg
        bmodgsm = get_var_data(bmod_var, at=times)
        bmodsm = cotran(bmodgsm, times, 'gsm2sm')
        bmodtilt = atan(bmodsm[*,2],sqrt(bmodsm[*,0]^2+bmodsm[*,1]^2))*deg
        dbtilt = btilt-bmodtilt
        dtilt_var = pre0+'db_tilt'
        vars.add, dtilt_var
        store_data, dtilt_var, times, dbtilt
        add_setting, dtilt_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: string_theta}
    endforeach
    vars = vars.toarray()
    
    print, 'Saved variables are: ', vars
    tplot_save, vars, filename=event.file
    print, 'Saved data to '+event.file+' ...'
    
    project.events[event_id] = event
    azim_prop_update_project, project
end


;+
; There is an overall hash to keep track of directory, event_ids.
; For each event, there is a key
;-
pro azim_prop_load_data, reload=reload, project=project, event_id=id

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    if n_elements(id) eq 0 then id = events.keys()

    foreach event_id, id do begin
        event = events[event_id]
        load = 0
        if ~event.haskey('file') then event.file = join_path([project.data_dir,event_id+'_basic_data.tplot'])
        if file_test(event.file) eq 0 then load = 1
        if keyword_set(reload) then load = 1
        if load eq 0 then begin
            tplot_restore, filename=event.file
        endif else begin
            file_delete, event.file, /allow_nonexistent
            azim_prop_load_data_per_event, project, event_id=event_id
        endelse
    endforeach

end


;ids = ['2014_0828_10']
;foreach id, ids do azim_prop_load_data, event_id=id, /reload
azim_prop_load_data, /reload
end
