;+
; Calculate the frequency spectrum at full frequency and selected frequency.
;-


pro test_identify_wake, date, probe=mission_probe

    test = 1
    secofday = 86400d

;---Load data.
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

;---Check if we need to load data.
    data_var = pre0+'e_uv'
    if tnames(data_var) ne '' then begin
        get_data, data_var, times
        index = lazy_where(times,full_time_range, count=count)
        if count le 1 then store_data, data_var, /delete
    endif
    if tnames(data_var) eq '' then load_var = 1 else load_var = 0
    if load_var then begin
        call_procedure, mission+'_read_euv', full_time_range, probe=probe
        get_data, data_var, times, data
    endif


    components = ['u','v']
    comp_index = 0
    comp_string = components[comp_index]
    comp_var = pre0+'e_'+comp_string
    if tnames(comp_var) ne '' then begin
        get_data, comp_var, times
        index = lazy_where(times,full_time_range, count=count)
        if count le 1 then store_data, comp_var, /delete
    endif
    if tnames(comp_var) eq '' then load_var = 1 else load_var = 0
    if load_var then begin
        stplot_index, data_var, comp_index, newnames=comp_var
        get_data, comp_var, times, data
        yrange = [-1,1]*max(abs(sg_autolim(data)))
        add_setting, comp_var, /smart, {$
            yrange: yrange, $
            short_name: 'UVW E!D'+comp_string+'!N', $
            display_type: 'scalar'}
    endif


;---Do wavelet analysis at the target freqs.
    cwt_var = comp_var+'_cwt'
    target = wake_effect_init_target_scales(spin_rate)
    target_periods = target['periods']
    target_scales = target['scales']
    ntarget_scale = target['nscale']
    load_cwt = load_var
    if load_cwt eq 0 then begin
        get_data, cwt_var, 0, cwt
        if n_elements(cwt.s_j) ne ntarget_scale then load_cwt = 1
    endif
    if load_cwt then begin
        calc_psd, comp_var, scales=target_scales
    endif


;---Plot settings.
    white = 255
    red = 6
    const_periods = spin_rate/[1.,2,3,4]
    xticklen = -0.06
    yticklen = -0.01
    unit = 'mV/m'
    xtitle = 'Period (sec)'
    xticklen = xticklen
    xlog = 1
    ytitle = '|'+unit+'|/Hz!U1/2!N'
    yticklen = yticklen
    ylog = 1
    yrange = [-3,3]
    ytickv = make_bins(yrange,1)
    yrange = 10.^minmax(ytickv)
    yticks = n_elements(ytickv)-1
    ytickn = strarr(yticks+1)
    for ii=0, yticks do begin
        case ytickv[ii] of
            0: ytickn[ii] = '1'
            else: ytickn[ii] = '10!U'+string(ytickv[ii],format='(I0)')
        endcase
    endfor
    if yticks ge 6 then begin
        target = 3
        for ii=0, yticks do begin
            if ii mod target eq 0 then continue
            ytickn[ii] = ' '
        endfor
    endif
    ytickv = 10.^ytickv
    yminor = 10

    xrange = [0.5,20]
    xtickv = [1.,10]

    label_size = 0.7
    psym = 1    ; dot.
    symsize = 0.5


;---Calculate the amplitude.
    amp_var = comp_var+'_amp'
    if load_cwt then begin
        nsection = round((full_time_range[1]-full_time_range[0])/duration)
        section_times = smkarthm(full_time_range[0], duration, nsection, 'x0')+duration*0.5
        amps = dblarr(nsection,ntarget_scale)

        get_data, comp_var, times, edata
        get_data, cwt_var, 0, cwt
        ww = cwt.w_nj
        s2t = cwt.s2t
        cdelta = cwt.cdelta
        dt = cwt.dt

        foreach time, section_times, ii do begin
            section_time_range = time+[-0.5,0.5]*duration
            index = lazy_where(times, section_time_range, count=N)
            uts = times[index]
            ees = edata[index]
            wps = abs(ww[index,*])^2
            gws = total(wps,1)/N
            psd = gws*2*s2t*dt/cdelta
            amps[ii,*] = sqrt(psd)
        endforeach
        store_data, amp_var, section_times, amps, target_periods
        add_setting, amp_var, /smart, {$
            ytitle: 'Period (sec)', $
            yrange: xrange, $
            zrange: [0.05,5], $
            zlog: ylog, $
            display_type: 'spec', $
            unit: '['+unit+'] /Hz!U1/2!N', $
            constant: const_periods, $
            const_color: white, $
            short_name: '|E|'}
        options, amp_var, 'no_interp', 0
        options, amp_var, 'x_no_interp', 1
        options, amp_var, 'y_no_interp', 0
    end


;---Loop through sections to tell data are good or bad.
    get_data, amp_var, section_times, amps
    nsection = n_elements(section_times)
    bad_data_flags = bytarr(nsection)   ; 1 for bad data.
    bg_index = [2,7,9,12]
    xf_index = [4,8,10]
    foreach time, section_times, ii do begin
        target_amp = reform(amps[ii,*])

    ;---Calculate spectral peak.
        xf_bg = interpol(target_amp[bg_index], target_periods[bg_index], target_periods[xf_index])
        xf_height = target_amp[xf_index]-xf_bg
        xf_ratio = target_amp[xf_index]/xf_bg

    ;---Identify peaks.
        min_peak_ratio = [2,1.3,1.3]
        xf_have_peak = xf_ratio ge min_peak_ratio
        have_1f_peak = xf_have_peak[0]
        have_2f_peak = xf_have_peak[1]
        have_3f_peak = xf_have_peak[2]

    ;---Tell good or bad.
        data_type = 'good data'
        if have_3f_peak or have_2f_peak then begin
            data_type = 'bad data'
        endif
        bad_data_flags[ii] = data_type ne 'good data'
    endforeach

    ; Flip isolated flags.
    ; Change 1-0-1 to 1-1-1.
    iso_index = list()
    for ii=1, nsection-2 do begin
        if bad_data_flags[ii] eq 1 then continue
        if bad_data_flags[ii-1]+bad_data_flags[ii+1] eq 2 then iso_index.add, ii
    endfor
    if n_elements(iso_index) gt 0 then begin
        iso_index = iso_index.toarray()
        bad_data_flags[iso_index] = 1
    endif
    ; Change 0-1-0 to 0-0-0.
    iso_index = list()
    for ii=1, nsection-2 do begin
        if bad_data_flags[ii] eq 0 then continue
        if bad_data_flags[ii-1]+bad_data_flags[ii+1] eq 0 then iso_index.add, ii
    endfor
    if n_elements(iso_index) gt 0 then begin
        iso_index = iso_index.toarray()
        bad_data_flags[iso_index] = 0
    endif

    ; Extend by 1 duration.
    index = where(bad_data_flags eq 1, count)
    for ii=0, count-1 do begin
        i0 = index[ii]-1 > 0
        i1 = index[ii]+1 < nsection-1
        bad_data_flags[i0:i1] = 1
    endfor
    ; Change 1-0-1 to 1-1-1.
    iso_index = list()
    for ii=1, nsection-2 do begin
        if bad_data_flags[ii] eq 1 then continue
        if bad_data_flags[ii-1]+bad_data_flags[ii+1] eq 2 then iso_index.add, ii
    endfor
    if n_elements(iso_index) gt 0 then begin
        iso_index = iso_index.toarray()
        bad_data_flags[iso_index] = 1
    endif

    ; Treat edges.
    bad_data_flags[0] = bad_data_flags[1]
    bad_data_flags[nsection-1] = bad_data_flags[nsection-2]

    flag_var = comp_var+'_wake_flag'
    store_data, flag_var, section_times, bad_data_flags
    add_setting, flag_var, /smart, {$
        display_type: 'scalar', $
        yrange: [-0.2,1.2], $
        yticks: 1, $
        ytickv: [0,1], $
        yminor: 0, $
        ytitle: '', $
        short_name: 'Bad E'}

    ; Convert to full time resolution.
    get_data, comp_var, times
    bad_data_flags = interpol(bad_data_flags, section_times, times) ne 0
    store_data, flag_var, times, bad_data_flags


;---Remove 2f and 3f.
    f1_index = [3,4,5]
    badf_index = [6,7,8,9,10,11,12]
    get_data, cwt_var, 0, cwt
    get_data, comp_var, times, edata, limit=lim

    f1_e = wavelet_reconstruct(cwt, index=f1_index)
    edata2 = f1_e
;    f23_e = wavelet_reconstruct(cwt, index=badf_index)
;    edata2 = edata-f23_e

    comp_var_new = comp_var+'_fixed'
    store_data, comp_var_new, times, edata2, limit=lim
    options, comp_var_new, 'labels', 'E fixed'
    options, comp_var_new, 'colors', 0

    comp_var_combo = comp_var+'_combo'
    store_data, comp_var_combo, times, [[edata],[edata2]], limit=lim
    options, comp_var_combo, 'labels', ['E orig','E fixed']
    options, comp_var_combo, 'colors', [6,0]

    stop

end


    events = list()

;    probe = 'polar'
;    date = '1998-09-25'
;    events.add, list(date,probe)

    probe = 'rbspb'
    date = '2013-06-07'
    events.add, list(date,probe)

    foreach event, events do begin
        date = event[0]
        probe = event[1]
        test_identify_wake, date, probe=probe
    endforeach

end
