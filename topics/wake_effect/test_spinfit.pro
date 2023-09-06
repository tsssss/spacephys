;+
; A sophisticated spinfit, adjusting the phase and amplitude b/w Eu and Ev.
;-


pro test_spinfit, date, probe=mission_probe

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
    spin_period = call_function(mission+'_info', 'spin_period')
    duration = 10*spin_period
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
    target = wake_effect_init_target_scales(spin_period)
    target_periods = target['periods']
    target_scales = target['scales']
    ntarget_scale = target['nscale']
    const_periods = spin_period/[1.,2,3,4]

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


;---Combine U and V.
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


;---Do spinfit.
    deg = 180d/!dpi
    rad = !dpi/180d

    get_data, fix_var, times, efix
    dt = sdatarate(times)
    nrec_per_spin = round(spin_period/dt)
    wt = 2*!dpi*findgen(nrec_per_spin)/nrec_per_spin
    coswt = cos(wt)
    sinwt = sin(wt)
    xx = [transpose(coswt),transpose(sinwt)]

    ntime = n_elements(times)
    nspinfit_section = floor(ntime/nrec_per_spin)
    spinfit_times = dblarr(nspinfit_section)
    spinfit_emags = fltarr(nspinfit_section)
    spinfit_ramps = fltarr(nspinfit_section)
    spinfit_dphis = fltarr(nspinfit_section)

    for ii=0, nspinfit_section-1 do begin
        i0 = nrec_per_spin*ii
        i1 = i0+nrec_per_spin-1
        tt = times[i0:i1]
        eu = efix[i0:i1,0]
        ev = efix[i0:i1,1]
        eu = eu-mean(eu)
        ev = ev-mean(ev)

        ;fitu = regress(xx, eu, const=const)
        fitu = la_least_squares(xx, eu)
        Ecosp = fitu[0]
        Esinp = -fitu[1]
        E = sqrt(Ecosp^2+Esinp^2)
        p = atan(Esinp, Ecosp)*deg

        fitv = la_least_squares(xx, ev)
        aEsinq = fitv[0]
        aEcosq = fitv[1]
        aE = sqrt(aEcosq^2+aEsinq^2)
        q = atan(aEsinq, aEcosq)*deg
        a = aE/E

        dphi = q-p
        if dphi gt  90 then dphi -= 180
        if dphi lt -90 then dphi += 180

        spinfit_times[ii] = mean(times[[i0,i1]])
        spinfit_emags[ii] = E
        spinfit_ramps[ii] = a
        spinfit_dphis[ii] = dphi
    endfor

    sf_emag_var = pre0+'spinfit_emag'
    store_data, sf_emag_var, spinfit_times, spinfit_emags
    sf_ramp_var = pre0+'spinfit_ramp'
    store_data, sf_ramp_var, spinfit_times, spinfit_ramps
    sf_dphi_var = pre0+'spinfit_dphi'
    store_data, sf_dphi_var, spinfit_times, spinfit_dphis

;---Generate plots.
    plot_dt = 3600*2
    nplot_time = floor(secofday/plot_dt)
    plot_times = time_double(date)+smkarthm(1800d,plot_dt,nplot_time,'x0')
    ticklen = -0.01
    xtitle = 'E!Du!N (mV/m)'
    ytitle = 'E!Dv!N (mV/m)'
    str_dphi = '!9'+string(68b)+'!X'+'!9'+string(102b)+'!X'
    tpos = [0.15,0.15,0.9,0.9]

    for ii=0, nplot_time-1 do begin
        the_time = plot_times[ii]
        plot_time_range = the_time+[0,spin_period]
        efix = get_var_data(fix_var, in=plot_time_range)
        eu = efix[*,0]
        ev = efix[*,1]
        nrec = n_elements(eu)
        dphi = get_var_data(sf_dphi_var, at=the_time)
        ramp = get_var_data(sf_ramp_var, at=the_time)
        emag = get_var_data(sf_emag_var, at=the_time)
        erange = sg_autolim([eu,ev,ev/ramp,emag])
        erange = abs(max(erange))*[-1,1]

        file = join_path([srootdir(),'plot','fig_e_spinfit_'+mission+probe+'_'+time_string(the_time,tformat='YYYY_MMDD_hhmm')+'.pdf'])
        if keyword_set(test) then file = test
        sgopen, file, xsize=5, ysize=5
        tmp = sgcalcpos(xchsz=xchsz,ychsz=ychsz)

        plot, erange, erange, /nodata, /isotropic, $
            xstyle=1, xrange=erange, xtitle=xtitle, xticklen=ticklen, $
            ystyle=1, yrange=erange, ytitle=ytitle, yticklen=ticklen, $
            position=tpos

        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*3
        color = sgcolor('red')
        plots, emag*coswt, emag*sinwt, color=color
        xyouts, tx,ty,/normal, 'Theoretical', color=color

        ty = tpos[3]-ychsz*2
        color = sgcolor('black')
        plots, shift(eu,-dphi/360*nrec), ev/ramp, psym=1, symsize=0.5, color=color
        xyouts, tx,ty,/normal, 'Corrected', color=color

        ty = tpos[3]-ychsz*1
        color = sgcolor('silver')
        plots, eu, ev, psym=1, symsize=0.5, color=color
        xyouts, tx,ty,/normal, 'Raw Data', color=color

        tx = tpos[0]
        msg = strupcase(strmid(fix_var,0,5))+'  '+time_string(plot_time_range[0])+' to '+$
            time_string(plot_time_range[1],tformat='hh:mm:ss')
        ty = tpos[3]+ychsz*1.3
        xyouts, tx,ty,/normal, msg

        msg = 'E = '+strtrim(string(emag,format='(F7.2)'),2)+' mV/m  '+$
        'E!Du!N/E!Dv!N = '+strtrim(string(ramp,format='(F7.2)'),2)+'  '+$
        str_dphi+' = '+strtrim(string(dphi,format='(F7.1)'),2)+' deg'
        ty = tpos[3]+ychsz*0.3
        xyouts, tx,ty,/normal, msg

        if keyword_set(test) then stop
        sgclose
    endfor


end


    events = list()

; Kazue example 1.
    probe = 'rbspa'
    date = '2012-10-06'
    probe = 'rbspb'
    date = '2013-05-01'
    probe = 'polar'
    date = '2003-05-01'
    events.add, list(date,probe)

    foreach event, events do begin
        date = event[0]
        probe = event[1]
        test_spinfit, date, probe=probe
    endforeach

end
