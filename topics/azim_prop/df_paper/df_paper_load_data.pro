;+
; Load data of r_gsm, b_gsm for each spacecraft. Calculate in-situ MLT, MLon/MLat.
; Calculate the tilt angle, bmod_gsm_t89, |B|-|B| T89.
;-



pro df_paper_load_data_per_event, events, id=id, reload=reload

    trace_dir = events['trace_dir']
    model = events['model']

    deg = events['constant'].deg
    rgb = events['constant'].rgb
    string_theta = events['constant'].string_theta

    event = events[id]
    time_range = event['time_range']
    probes = event['probes']
    event_id = event['event_id']


    if ~event.haskey('file') then event['file'] = join_path([events['data_dir'],events['name']+'_'+event_id+'_data.tplot'])

    ; Force to reload data.
    if keyword_set(reload) then begin
        file_delete, event.file, /allow_nonexistent
    endif


    if file_test(event.file) then begin
        ; Because every event has the same vars, we have to load file to overwrite the vars of the old event.
        tplot_restore, filename=event.file
    endif else begin
        vars = list()
        foreach probe, probes do begin
            mission = resolve_probe(probe)
            probe = mission.name+mission.probe
            pre0 = probe+'_'
            label = strupcase(mission.short_name+mission.probe)

        ;---Load xxx_r_gsm.
            call_procedure, mission.routine_name+'_read_orbit', time_range, probe=mission.probe
            rgsm_var = pre0+'r_gsm'
            if tnames(rgsm_var) eq '' then stop
            vars.add, rgsm_var

            ; Convert xxx_r_gsm to xxx_r_sm.
            get_data, rgsm_var, times, rgsm
            rsm = cotran(rgsm, times, 'gsm2sm')
            rsm_var = pre0+'r_sm'
            vars.add, rsm_var
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
            vars.add, mlat_var
            store_data, mlat_var, times, mlat
            add_setting, mlat_var, /smart, {$
                display_type: 'scalar', $
                unit: 'deg', $
                short_name: 'MLat'}

            mlon_var = pre0+'mlon'
            vars.add, mlon_var
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


        ;---Load xxx_bgsm.
            call_procedure, mission.routine_name+'_read_bfield', time_range, probe=mission.probe
            bgsm_var = pre0+'b_gsm'
            if tnames(bgsm_var) eq '' then stop
            vars.add, bgsm_var

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
            tilt = atan(bsm[*,2],sqrt(bsm[*,0]^2+bsm[*,1]^2))*deg
            tilt_var = pre0+'b_tilt'
            vars.add, tilt_var
            store_data, tilt_var, times, tilt
            add_setting, tilt_var, /smart, {$
                display_type: 'scalar', $
                unit: 'deg', $
                short_name: string_theta}
        endforeach
        event.vars = vars
        events[id] = event
        store_data, events['var'], 0, events
        tplot_save, events['var'], filename=events['file']

        ; Tune options.
        vars = list()
        foreach probe, probes do vars.add, probe+'_'+['b_gsm','bmod_gsm','r_gsm','r_sm'], /extract
        vars = vars.toarray()
        options, vars, 'colors', rgb

        tplot_save, event.vars.toarray(), filename=event.file
    endelse
end


function df_paper_default_event_info, event

    event_info = hash()
    id = '2014_0828_10'
    event_info[id] = dictionary($
        'event_id', id, $
        'time_range', time_double(['2014-08-28/09:50','2014-08-28/11:00']), $
        'probes', ['rbsp'+['b'],'th'+['a','d','e'],'g'+['13','15']])

    id = '2012_1001_02'
    event_info[id] = dictionary($
        'event_id', id, $
        'time_range', time_double(['2012-10-01/01:30','2012-10-01/03:30']), $
        'probes', ['th'+['a','d','e'],'g'+['13','15']])

    id = '2013_0607_05'
    event_info[id] = dictionary($
        'event_id', id, $
        'time_range', time_double(['2013-06-07/04:00','2013-06-07/07:00']), $
        'probes', ['rbsp'+['a','b'],'th'+['d','e'],'g'+['13','15']])

    return, event_info
end

;+
; There is an overall hash to keep track of directory, event_ids.
; For each event, there is a key
;-
pro df_paper_load_data, reload=reload, event=event, id=id

    event = hash()
    event['name'] = 'df_paper'
    event['var'] = event['name']+'_event_info'
    event['root_dir'] = sparentdir(srootdir())
    event['data_dir'] = join_path([event['root_dir'],'data'])
    event['plot_dir'] = join_path([event['root_dir'],'plot'])
    event['file'] = join_path([event['data_dir'],event['name']+'_event_info.tplot'])
    event['event_info'] = df_paper_default_event_info(event)
    event['model'] = 't89'
    event['trace_dir'] = -1

    constant = dictionary()
    constant.string_theta = '!9'+string(113b)+'!X'
    constant.rgb = sgcolor(['red','green','blue'])
    constant.re = 6378d
    constant.deg = 180d/!dpi
    constant.rad = !dpi/180d
    event['constant'] = constant


    ; Update event info.
    if file_test(event['file']) then begin
        tplot_restore, filename=event['file']
        event = get_var_data(event['var'])
    endif else begin
        store_data, event['var'], 0, event
        tplot_save, event['var'], filename=event['file']
    endelse


    ; Force to reload event.
    if keyword_set(reload) then begin
        if event.haskey(id) then begin
            file_delete, event[id].file, /allow_nonexistent
            event.remove, id
            (event['ids']).remove, (event['ids']).where(id)
        endif else if ~event.haskey('ids') then event['ids'] = list()
    endif

    ; Update variable for the wanted event.
    if n_elements(id) eq 0 then return
    if event.haskey(id) then begin
        tplot_restore, filename=event[id].file
    endif else begin
        event['ids'].add, id
        event_info = (event['event_info'])[id]
        event[id] = dictionary('event_id', id)
        event[id].time_range = event_info.time_range
        event[id].probes = event_info.probes

        df_paper_load_data_per_event, event, id=id, reload=reload

        store_data, event['var'], 0, event
        tplot_save, event['var'], filename=event['file']
    endelse

end


ids = ['2012_1001_02','2013_0607_05','2014_0828_10']
;ids = ['2013_0828_01', '2013_1002_20', '2013_1009_08', '2014_0828_10', '2014_0912_19', '2016_1013_12', '2016_1025_08', '2016_1028_23', '2017_0301_15', '2017_0328_10']
foreach id, ids do df_paper_load_data, event=event, id=id, /reload
end
