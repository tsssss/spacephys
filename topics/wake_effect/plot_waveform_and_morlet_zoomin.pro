;+
; Calculate the frequency spectrum at full frequency and selected frequency.
;-


pro plot_wavelet_and_morlet_zoomin, zoomin_times, probe=mission_probe

    test = 1
    secofday = 86400d

;---Load data.
    if n_elements(zoomin_times) eq 0 then message, 'No input times ...'
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
    date = time_string(zoomin_times[0],tformat='YYYY-MM-DD')
    full_time_range = time_double(date)+[0,secofday]

    components = ['u','v']
    comp_index = 0
    comp_string = components[comp_index]

;---Check if we need to load data.
    data_var = pre0+'e_uv'
    comp_var = pre0+'e_'+comp_string
    if tnames(data_var) ne '' then begin
        get_data, data_var, times
        index = lazy_where(times,full_time_range, count=count)
        if count le 1 then store_data, data_var, /delete
    endif
    if tnames(data_var) eq '' then load_var = 1 else load_var = 0

    if load_var then begin
        call_procedure, mission+'_read_euv', full_time_range, probe=probe
        get_data, data_var, times, data
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
            zrange: [0.5,50], $
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


;---Plot.
    get_data, amp_var, section_times, amps
    foreach time, zoomin_times do begin
        tmp = min(section_times-(time+0.5*duration), /absolute, index)
        target_amp = reform(amps[index,*])
        time = section_times[index]
        time_range = time+[-0.5,0.5]*duration

    ;---Calculate the full psd.
        get_data, comp_var, times, edata
        index = lazy_where(times, time_range)
        tmp_var = comp_var+'_tmp'
        store_data, tmp_var, times[index], edata[index]
        yrange = sg_autolim(edata[index])
        yrange = [-1,1]*max(abs(yrange))
        tmp_yminor = 2
        options, tmp_var, 'constant', 0
        options, tmp_var, 'labels', ''
        options, tmp_var, 'xticklen', xticklen
        options, tmp_var, 'yticklen', yticklen
        options, tmp_var, 'yticks', 4
        options, tmp_var, 'yminor', tmp_yminor
        options, tmp_var, 'ystyle', 1
        options, tmp_var, 'ytitle', '('+unit+')'
        ylim, tmp_var, yrange
        calc_psd, tmp_var


    ;---Plot.
        file = join_path([srootdir(),'plot','fig_e'+comp_string+'_wavelet_'+mission+probe+'_'+time_string(time,tformat='YYYY_MMDD_hhmm')+'.pdf'])
        if keyword_set(test) then file = test
        sgopen, file, xsize=8, ysize=2

        poss = sgcalcpos(1,2,xpans=[3,1.2], xpad=10, $
            tmargin=2, bmargin=4, lmargin=11, rmargin=2, xchsz=xchsz, ychsz=ychsz)

        ; Panel 1.
        tpos = poss[*,0]
        tplot, tmp_var, position=tpos, trange=time_range, vlab_margin=10
        store_data, tmp_var, /delete
        tx = tpos[0]
        ty = tpos[3]+ychsz*0.35
        xyouts, tx,ty,/normal, 'a. Waveform E!D'+comp_string+'!N '+strupcase(mission)

        ; Panel 2.
        tpos = poss[*,1]
        psd_var = tmp_var+'_psd_cwt'
        get_data, psd_var, xs, ys
        ys = sqrt(ys)
        xs = 1d/xs
        xtitle = 'Period (sec)'
        xticklen = xticklen
        xlog = 1
        ytitle = '|'+unit+'|/Hz!U1/2!N'
        yticklen = yticklen
        ylog = 1
        yrange = minmax(alog10(ys))
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

        plot, xs, ys, /nodata, $
            xlog=xlog, xstyle=1, xticklen=xticklen, xtitle=xtitle, $
            ylog=ylog, ystyle=1, yticklen=yticklen, ytitle=ytitle, yrange=yrange, $
            yticks=yticks, ytickv=ytickv, ytickname=ytickn, yminor=yminor, $
            /noerase, position=tpos
        oplot, xs,ys, color=sgcolor('silver')

        tx = tpos[0]
        ty = tpos[3]+ychsz*0.35
        xyouts, tx,ty,/normal, 'b. Amplitude'

        ; Add the psd at target freqs.
        xs = target_periods
        ys = target_amp
        oplot, xs,ys, color=sgcolor('blue'), psym=psym, symsize=symsize

        ; Add constant freqs.
        foreach number, [1d,2,3,4] do begin
            period = spin_rate/number
            plots, period+[0,0], yrange, linestyle=1
            msg = (number eq 1)? 'f  ': string(number,format='(I0)')+'f'
            tmp = convert_coord(period,yrange[1], /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*1
            xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size
            if number eq 1 then xyouts, tx,ty,/normal, alignment=0.1, 'spin', charsize=label_size*0.8
        endforeach

        if keyword_set(test) then stop
    endforeach

    if keyword_set(test) then stop
    sgclose
end


    events = list()

    probe = 'polar'
    date = '1998-09-14'
    ;times = ['21:10','21:15',21:20']   ; 2f.
    ;times = ['14:00']   ; 3f.
    date = '1998-09-25'
    times = ['12:40']
    times = time_double(date+'/'+times)
    events.add, list(times,probe)

    foreach event, events do begin
        times = event[0]
        probe = event[1]
        plot_wavelet_and_morlet_zoomin, times, probe=probe
    endforeach

end
