;+
; Load Eu and Ev, compare the waveforms of the two, to identify bad efield.
;   1. Wake effect.
;   2. Spikes.
;-


pro test_analyze_euv, date, probe=mission_probe

    test = 1
    secofday = 86400d

;---Prepare to load data.
    if n_elements(date) eq 0 then message, 'No input date ...'
    if n_elements(mission_probe) eq 0 then message, 'No input probe ...'
    case strmid(mission_probe,0,2) of
        'po': begin
            mission = 'polar'
            probe = ''
            pre0 = 'po_'
            end
        'rb': begin
            mission = 'rbsp'
            probe = strmid(mission_probe,0,1,/reverse_offset)
            pre0 = 'rbsp'+probe+'_'
            end
        'th': begin
            mission = 'themis'
            probe = strmid(mission_probe,0,1,/reverse_offset)
            pre0 = 'th'+probe+'_'
            end
    endcase
    spin_rate = call_function(mission+'_info', 'spin_rate')
    duration = 10*spin_rate
    full_time_range = time_double(date)+[0,secofday]
    timespan, full_time_range, secofday, /seconds

    nsection = round((full_time_range[1]-full_time_range[0])/duration)
    section_times = smkarthm(full_time_range[0], duration, nsection, 'x0')+duration*0.5

;---Check if we need to load data.
    data_var = pre0+'e_uv'
    if tnames(data_var) ne '' then begin
        get_data, data_var, times
        index = where_pro(times,full_time_range, count=count)
        if count le 1 then store_data, data_var, /delete
    endif
    if tnames(data_var) eq '' then load_var = 1 else load_var = 0
    if load_var then begin
        call_procedure, mission+'_read_euv', full_time_range, probe=probe
    endif

    ; Split into components.
    components = ['u','v']
    ncomponent = n_elements(components)
    comp_vars = pre0+'e_'+components
    load_var = 0
    foreach comp_var, comp_vars do if tnames(comp_var) eq '' then load_var = 1
    if load_var then begin
        stplot_split, data_var, newnames=comp_vars
        get_data, data_var, times, data
        yrange = [-1,1]*max(abs(sg_autolim(data)))
        for ii=0, ncomponent-1 do begin
            comp_var = comp_vars[ii]
            add_setting, comp_var, /smart, {$
                yrange: yrange, $
                short_name: 'UVW E!D'+components[ii]+'!N', $
                display_type: 'scalar'}
            ; Remove cwt to let it refresh.
            cwt_var = comp_var+'_cwt'
            if tnames(cwt_var) ne '' then store_data, cwt_var, /delete
        endfor
    endif


;---Work on each component.
    target = wake_effect_init_target_scales(spin_rate)
    target_periods = target['periods']
    target_scales = target['scales']
    ntarget_scale = target['nscale']
    const_periods = spin_rate/[1.,2,3,4]

    f1_index = [1,2,3,4,5,6]

    red = 6
    black = 0
    white = 255

    period_range = [1,20]
    period_tickv = [1.,10]
    amp_log = 1
    unit = 'mV/m'

    foreach comp_var, comp_vars, jj do begin
        get_data, comp_var, times, data, limits=lim

    ;---Calculate CWT.
        load_cwt = 0
        cwt_var = comp_var+'_cwt'
        amp_var = comp_var+'_amp'
        if tnames(cwt_var) ne '' then begin
            get_data, cwt_var, 0, cwt
            if n_elements(cwt.s_j) ne ntarget_scale then load_cwt = 1
        endif else load_cwt = 1
        if load_cwt then begin
            calc_psd, comp_var, scales=target_scales

            ; Calculate the amplitude spectrogram.
            amps = dblarr(nsection,ntarget_scale)
            get_data, comp_var, times, edata
            get_data, cwt_var, 0, cwt
            ww = cwt.w_nj
            s2t = cwt.s2t
            cdelta = cwt.cdelta
            dt = cwt.dt

            foreach time, section_times, ii do begin
                section_time_range = time+[-0.5,0.5]*duration
                index = where_pro(times, section_time_range, count=N)
                wps = abs(ww[index,*])^2
                gws = total(wps,1)/N
                psd = gws*2*s2t*dt/cdelta
                amps[ii,*] = sqrt(psd)
            endforeach
            store_data, amp_var, section_times, amps, target_periods
            add_setting, amp_var, /smart, {$
                ytitle: 'Period (sec)', $
                yrange: period_range, $
                zrange: [0.05,5], $
                zlog: amp_log, $
                display_type: 'spec', $
                unit: '['+unit+'] /Hz!U1/2!N', $
                constant: const_periods, $
                const_color: white, $
                short_name: '|E|'}
            options, amp_var, 'no_interp', 0
            options, amp_var, 'x_no_interp', 1
            options, amp_var, 'y_no_interp', 0
        endif

    ;---Calculate the corrected data.
        fix_var = comp_var+'_fix'
        load_var = check_if_update(fix_var, full_time_range)
        if load_var then begin
            get_data, comp_var, times, data, limits=lim
            fix_data = wavelet_reconstruct(cwt, index=f1_index)
            store_data, fix_var, times, fix_data, limits=lim
            combo_var = comp_var+'_combo'
            store_data, combo_var, times, [[data],[fix_data]], limits=lim
            options, combo_var, 'labels', 'E!D'+components[jj]+'!N '+['orig', 'fixed']
            options, combo_var, 'colors', [red,black]
        endif

    ;---Calculate the spectrogram for the corrected data.
        fix_amp_var = fix_var+'_amp'
        load_var = check_if_update(fix_amp_var)
        if load_var then begin
            fix_cwt_var = fix_var+'_cwt'
            amps = dblarr(nsection,ntarget_scale)
            get_data, fix_var, times, edata
            calc_psd, fix_var, scales=target_scales
            get_data, fix_cwt_var, 0, cwt
            ww = cwt.w_nj
            s2t = cwt.s2t
            cdelta = cwt.cdelta
            dt = cwt.dt

            foreach time, section_times, ii do begin
                section_time_range = time+[-0.5,0.5]*duration
                index = where_pro(times, section_time_range, count=N)
                wps = abs(ww[index,*])^2
                gws = total(wps,1)/N
                psd = gws*2*s2t*dt/cdelta
                amps[ii,*] = sqrt(psd)
            endforeach
            get_data, amp_var, limits=lim
            store_data, fix_amp_var, section_times, amps, target_periods, limits=lim
        endif
    endforeach


;---Combine U and V, convert to GSM.
    fix_var = pre0+'e_fix'
    load_var = check_if_update(fix_var, full_time_range)
    if load_var then begin
        rgb = sgcolor(['red','green','blue'])
        stplot_merge, comp_vars+'_fix', newname=fix_var
        get_data, fix_var, times, data
        nrec = n_elements(times)
        data = [[data],[fltarr(nrec)]]
        store_data, fix_var, times, data
        add_setting, fix_var, /smart, {$
            display_type: 'vector', $
            short_name: 'E', $
            unit: unit, $
            colors: rgb, $
            coord: 'UVW', $
            coord_labels: ['u','v','w']}
    endif
    
    
;    cotrans_set_coord, dl, 'uvw'
;    store_data, fix_var, dlimit=dl
;    rbsp_spinfit, fix_var

    egsm_var = pre0+'e_gsm'
;    rbsp_uvw2gsm, fix_var, egsm_var
;    uniform_time, egsm_var, spin_rate

    emgse_var = pre0+'e_mgse'
    get_data, egsm_var, times, egsm
    emgse = cotran(egsm, times, 'gsm2mgse')
    store_data, emgse_var, times, emgse
    add_setting, emgse_var, /smart, {$
        display_type: 'vector', $
        short_name: 'E', $
        unit: unit, $
        colors: rgb, $
        coord: 'MGSE', $
        coord_labels: ['x','y','z']}
    stop
    
    rbsp_read_efield, full_time_range, probe=probe, resolution='survey'
end

    events = list()

; Kazue example 1.
    probe = 'rbspa'
    date = '2012-10-06'
    events.add, list(date,probe)

    foreach event, events do begin
        date = event[0]
        probe = event[1]
        test_analyze_euv, date, probe=probe
    endforeach

end
