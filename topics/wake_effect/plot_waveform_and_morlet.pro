;+
; Plot waveform and morlet wavelet for a whole day, to show overall features in time-frequency domain.
;-


pro plot_wavelet_and_morlet, date, probe=mission_probe

    test = 0
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
    psym = 3    ; dot.
    symsize = 0.1


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


;---Plot.
    file = join_path([srootdir(),'plot','fig_e'+comp_string+'_wavelet_'+mission+probe+'_'+time_string(full_time_range[0],tformat='YYYY_MMDD')+'.pdf'])
    if keyword_set(test) then file = test
    sgopen, file, xsize=10, ysize=4
    device, decomposed=0
    loadct2, 43

    vars = [comp_var,amp_var]
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, ypans=[2,1.5], rmargin=8, tmargin=1, xchsz=xchsz,ychsz=ychsz)

    tvar = comp_var
    ypans = [1,2,1]
    lin_yrange = [-1,1]*10
    lin_ytickv = [-2,-1,0,1,2]*5
    lin_yticks = n_elements(lin_ytickv)-1
    lin_ytickn = string(lin_ytickv,format='(I0)')
    lin_yminor = 5

    log_yrange = [-1,1]*1000
    log_ytickv = 10.^[1,2,3]
    log_yticks = n_elements(log_ytickv)-1
    log_ytickn = string(log_ytickv,format='(I0)')
    log_yminor = 10

    tpos = poss[*,0]
    tx = xchsz*0.5
    ty = tpos[3]-ychsz*0.5
    xyouts, tx,ty,/normal, 'a. E!D'+comp_string
    tposs = sgcalcpos(3, ypans=ypans, ypad=0, position=tpos)
    get_data, tvar, xx, yy, limits=lim

    ; negative part.
    tpos = tposs[*,2]
    yrange = -[min(log_yrange),min(lin_yrange)]
    ylog = 1
    ytickn = '-'+log_ytickn & ytickn[0] = ' '
    store_data, tvar+'_neg', xx, -yy, limits=lim
    store_data, tvar+'_neg', limits={$
        ytitle: '', $
        labels: '', $
        xticklen: xticklen, $
        ytickv: log_ytickv, $
        yticks: log_yticks, $
        yminor: log_yminor, $
        ytickname: ytickn, $
        yrange: yrange, $
        yticklen: yticklen, $
        ylog: ylog}
    tplot, tvar+'_neg', trange=full_time_range, position=tpos, /noerase, /nouttick
    store_data, tvar+'_neg', /delete

    ; posivie part.
    tpos = tposs[*,0]
    yrange = [max(lin_yrange),max(log_yrange)]
    ylog = 1
    ytickn = log_ytickn & ytickn[0] = ' '
    store_data, tvar, limits={$
        ytitle: '', $
        labels: '', $
        xticklen: xticklen, $
        ytickv: log_ytickv, $
        yticks: log_yticks, $
        yminor: log_yminor, $
        ytickname: ytickn, $
        yrange: yrange, $
        yticklen: yticklen, $
        ylog: ylog}
    tplot, tvar, trange=full_time_range, position=tpos, /noerase, /nouttick

    ; middle part.
    tpos = tposs[*,1]
    yrange = lin_yrange
    ylog = 0
    ytickn = lin_ytickn
    yspace = (yticklen ge 0)? 0: yticklen
    tx = tpos[0]-xchsz*0.15+yspace
    ty = tpos[1]-ychsz*0.25
    xyouts, tx,ty,/normal,alignment=1, lin_ytickn[0]
    ty = tpos[3]-ychsz*0.25
    xyouts, tx,ty,/normal,alignment=1, lin_ytickn[lin_yticks]

    ytickn[[0,lin_yticks]] = ' '
    store_data, tvar, limits=lim
    store_data, tvar, limits={$
        yrange: yrange, $
        ytitle: '('+unit+')', $
        labels: '', $
        xticklen: xticklen, $
        ytickv: lin_ytickv, $
        yticks: lin_yticks, $
        yminor: lin_yminor, $
        ytickname: ytickn, $
        yticklen: yticklen, $
        ylog: ylog}
    tplot, tvar, trange=full_time_range, position=tpos, /noerase, /nouttick
    tx = tpos[2]+xchsz*0.5
    ty = (tpos[1]+tpos[3])*0.5
    msg = strupcase(mission)
    if mission ne 'polar' then msg += '-'+strupcase(probe)
    msg += '!CUVW E!D'+comp_string
    xyouts, tx,ty,/normal, msg
    plots, tpos[[0,2]],(tpos[1]+tpos[3])*0.5+[0,0],/normal, linestyle=1, color=white


;---Secend panel.
    tpos = poss[*,1]
    tx = xchsz*0.5
    ty = tpos[3]-ychsz*0.5
    xyouts, tx,ty,/normal, 'b. Morlet'
    options, amp_var, 'yticklen', yticklen
    tplot, amp_var, trange=full_time_range, position=tpos, /noerase
    get_data, amp_var, limits=lim
    plot, full_time_range, lim.yrange, /nodata, /noerase, position=tpos, $
        xstyle=5, xlog=0, ystyle=5, ylog=1
    foreach period, const_periods do begin
        tx = tpos[0]-xchsz*0.5
        tmp = convert_coord(full_time_range[0],period, /data, /to_normal)
        ty = tmp[1]-ychsz*0.2
        msg = string(round(spin_rate/period),format='(I0)')
        if msg eq '1' then msg = 'f' else msg += 'f'
        xyouts, tx,ty,/normal, alignment=1, msg, charsize=label_size;, color=white
    endforeach


    if keyword_set(test) then stop
    sgclose

end

    events = list()

;    probe = 'polar'
;    dates = ['1998-09-25','1998-09-26','1998-09-14']
;    foreach date, dates do events.add, [date,probe]

;    probe = 'rbspb'
;    dates = ['2013-05-01','2013-06-07']
;    foreach date, dates do events.add, [date,probe]

    probe = 'rbspa'
    dates = ['2012-11-14','2013-06-07']
    dates = ['2012-10-06']
    foreach date, dates do events.add, [date,probe]

    foreach event, events do begin
        date = event[0]
        probe = event[1]
        plot_wavelet_and_morlet, date, probe=probe
    endforeach

end
