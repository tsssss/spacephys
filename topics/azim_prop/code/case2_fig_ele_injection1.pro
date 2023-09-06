;+
; Study the 1st electron injection.
;-


;---Load electron flux.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 0


;---Settings.
    einfo = get_var_data(info_var)
    time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
    plot_time_range = time_double(['2014-08-28/10:00','2014-08-28/10:40'])
    plot_time = time_double('2014-08-28/10:10')
    species = 'e'
    ; This is the target time range for xcor.
    test_time_range = time_double(['2014-08-28/10:05','2014-08-28/10:15'])
    left_label = 'b'
    right_label = 'c'
    
    ; Probe color
    all_probes = ['thd','g15','1991-080','g13','1994-084','LANL-01A','LANL-02A','LANL-04A']
    ct = 40
    top_color = 250
    bottom_color = 100
    nprobe = n_elements(all_probes)
    index_colors = bytscl(findgen(nprobe), top=top_color-bottom_color)+bottom_color
    probe_colors = fltarr(nprobe)
    foreach color, index_colors, color_id do $
        probe_colors[color_id] = sgcolor(color, ct=ct)


;---Get the flux in the wanted energy range, and the drift periods.
    probes = ['thd','g15','1991-080','g13','1994-084']
    nprobe = n_elements(probes)
    index = intarr(nprobe)
    foreach probe, probes, probe_id do index[probe_id] = where(all_probes eq probe)
    probe_colors = probe_colors[index]
    
    foreach probe, probes do begin
        prefix = probe+'_'
        the_var = prefix+'kev_'+species+'_flux'
        get_data, the_var, times, fluxes, energies

        ; Filter energy.
        prefix2 = prefix+species+'_'
        drift_period_var = prefix2+'drift_period'
        get_data, drift_period_var, tmp, tmp, target_energies
        ntarget_energy = n_elements(target_energies)
        index = fltarr(n_elements(target_energies))
        for ii=0,ntarget_energy-1 do begin
            index[ii] = where(energies eq target_energies[ii])
        endfor
        energies = energies[index]
        fluxes = fluxes[*,index]
        ; Filter time.
        index = where_pro(times, '[]', time_range)
        fluxes = fluxes[index,*]
        times = times[index]

        ; Store new flux.
        yrange = alog10(minmax(fluxes))
        yrange = 10.^[floor(yrange[0]),ceil(yrange[1])]
        log_ytickv = make_bins(alog10(yrange),1)
        ytickv = 10.^log_ytickv
        yticks = n_elements(ytickv)-1
        yminor = 10
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        index = where(log_ytickv eq 0, count)
        if count ne 0 then ytickn[index] = '1'
        index = where(log_ytickv eq 1, count)
        if count ne 0 then ytickn[index] = '10'
        index = where(log_ytickv eq -1, count)
        if count ne 0 then ytickn[index] = '0.1'
        if yticks ge 5 then ytickn[0:*:2] = ' '

        var = prefix+'kev_'+species
        store_data, var, times, fluxes, energies
        add_setting, var, /smart, {$
            display_type: 'list', $
            ylog: 1, $
            yrange: yrange, $
            ytickv: ytickv, $
            yticks: yticks, $
            ytickname: ytickn, $
            yminor: yminor, $
            color_table: 52, $
            unit: '#/cm!U2!N-s-sr-keV', $
            value_unit: 'keV', $
            short_name: 'e!U-!N flux'}
    endforeach


;---For each probe, find the optimal dMLT to maximize xcor among fluxes.
    mlts = fltarr(nprobe)
    dmlt = 0.1
    test_mlts = make_bins([-0.5,4.5],dmlt)
    test_xcors = dblarr(n_elements(test_mlts))
    foreach probe, probes, probe_id do begin
        prefix = probe+'_'
        prefix2 = prefix+species+'_'
        drift_period_var = prefix2+'drift_period'

        new_var = prefix+'kev_'+species+'_new_flux'
        if tnames(new_var) ne '' then begin
            the_test_time = get_setting(new_var, 'test_time')
;if ~keyword_set(test) then if the_test_time eq plot_time then continue
            if total(the_test_time-test_time_range) eq 0 then continue
        endif
        get_data, prefix+'kev_'+species, limits=lim


        foreach test_mlt, test_mlts, jj do begin
            ; Get the MLT.
            r_gsm = get_var_data(prefix+'r_gsm', at=test_time_range)
            r_sm = cotran(r_gsm, test_time_range, 'gsm2sm')
            mlts[probe_id] = mean(pseudo_mlt(r_sm))

            the_var = prefix+'kev_'+species
            get_data, the_var, times, fluxes, energies, limits=lim

            drift_periods = get_var_data(drift_period_var, at=times)
            new_fluxes = fluxes
            foreach energy, energies, ii do begin
                dtimes = test_mlt/24*drift_periods[*,ii]
                new_times = times-dtimes
                new_fluxes[*,ii] = interpol(fluxes[*,ii], new_times, times)
            endforeach
            store_data, new_var, times, new_fluxes, energies, limits=lim
            tplot, [the_var, new_var]

            new_fluxes = get_var_data(new_var, in=test_time_range)
            index = where(finite(new_fluxes,/nan), count)
            if count ne 0 then new_fluxes[index] = 0
            xcor = 0
            foreach energy, energies, ii do begin
                if ii eq 0 then continue
                xcor += c_correlate(new_fluxes[*,ii-1],new_fluxes[*,ii],0)
            endforeach
            test_xcors[jj] = xcor/(n_elements(energies)-1)
        endforeach
        max_xcor = max(test_xcors, index, /nan)
        test_mlt = test_mlts[index]>0


        prefix = probe+'_'
        prefix2 = prefix+species+'_'
        drift_period_var = prefix2+'drift_period'

        the_var = prefix+'kev_'+species
        get_data, the_var, times, fluxes, energies, limits=lim

        drift_periods = get_var_data(drift_period_var, at=times)
        new_fluxes = fluxes
        foreach energy, energies, ii do begin
            dtimes = test_mlt/24*drift_periods[*,ii]
            new_times = times-dtimes
            new_fluxes[*,ii] = interpol(fluxes[*,ii], new_times, times)
        endforeach
        store_data, new_var, times, new_fluxes, energies
        options, new_var, 'test_time', test_time_range
        options, new_var, 'dmlt', test_mlt

        tplot, [the_var, new_var]
;if keyword_set(test) then stop
    endforeach


    dmlts = fltarr(nprobe)
    foreach probe, probes, probe_id do begin
        prefix = probe+'_'
        new_var = prefix+'kev_'+species+'_new_flux'
        dmlts[probe_id] = get_setting(new_var, 'dmlt')
    endforeach




;---Caclulate the MLT of dispersionless injection.
    new_mlts = fltarr(nprobe)
    foreach probe, probes, probe_id do begin
        prefix = probe+'_'
        prefix2 = prefix+species+'_'

        new_var = prefix+'kev_'+species+'_new_flux'
        get_data, new_var, new_times, new_fluxes, energies
        ; Get drift period.
        drift_period_var = prefix2+'drift_period'
        drift_periods = get_var_data(drift_period_var, at=new_times)
        tmp = min(new_times-plot_time, time_index, /absolute)
        new_time = new_times[time_index]

        nenergy = n_elements(energies)
        the_mlts = fltarr(nenergy)
        pos_var = prefix+'r_gsm'
        dmlt = dmlts[probe_id]
        foreach energy, energies, energy_id do begin
            dtime = dmlt/24*drift_periods[time_index,energy_id]
            the_time = new_time+dtime
            r_gsm = get_var_data(pos_var, at=the_time)
            r_sm = cotran(r_gsm, the_time, 'gsm2sm')
            the_mlts[energy_id] = pseudo_mlt(r_sm)
        endforeach
        new_mlts[probe_id] = mean(the_mlts)-dmlt
    endforeach

    ; The orig MLT.
    old_mlts = fltarr(nprobe)
    foreach probe, probes, probe_id do begin
        prefix = probe+'_'
        pos_var = prefix+'r_gsm'
        r_gsm = get_var_data(pos_var, at=plot_time)
        r_sm = cotran(r_gsm, plot_time, 'gsm2sm')
        old_mlts[probe_id] = pseudo_mlt(r_sm)
    endforeach




    ; Prepare to plot.
;    foreach probe, probes do begin
;        get_data, probe+'_kev_'+species, limits=lim
;        store_data, probe+'_kev_'+species+'_new_flux', limits=lim
;    endforeach


    root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
    ofn = join_path([root_dir,'case2_fig_ele_injection1.pdf'])
    if keyword_set(test) then ofn = 0

    pan_xsize = 3
    pan_ysize = 1
    xpad = 4
    ypad = 0.4
    margins = [5,4,8,2]
    sgopen, 0, xsize=1, ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    nypan = nprobe
    nxpan = 2
    fig_ysize = pan_ysize*nypan+abs_ychsz*(total(margins[[1,3]])+ypad*(nypan-1))
    fig_xsize = pan_xsize*nxpan+abs_xchsz*(total(margins[[0,2]])+xpad*(nxpan-1))

    tplot_options, 'version', 2

    sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize
    poss = sgcalcpos(nypan,nxpan, margins=margins, xchsz=xchsz, ychsz=ychsz, xpad=xpad, ypad=ypad)

    left_poss = reform(poss[*,0,*])
    vars = probes+'_kev_'+species
    options, vars, 'labels', ' '
    options, vars, 'ytitle', ' '
    tplot, vars, trange=plot_time_range, position=left_poss, /noerase, /novtitle
    timebar, plot_time, linestyle=1

    ; Add figure label.
    for ii=0, nprobe-1 do begin
        tpos = left_poss[*,ii]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,/normal, left_label+'-'+string(ii+1,format='(I0)')+'.'
    endfor

    ; Add title.
    tpos = left_poss[*,0]
    tx = total(tpos[[0,2]])*0.5
    ty = tpos[3]+ychsz*0.5
    xyouts, tx,ty,/normal,alignment=0.5, 'Electron flux (#/cm!U2!n-s-sr-keV) in real time'

    ; Add new MLT.
    for ii=0, nprobe-1 do begin
        tpos = left_poss[*,ii]
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,/normal, alignment=1, string(old_mlts[ii],format='(F5.1)')+' MLT'
    endfor


    right_poss = reform(poss[*,1,*])
    vars = probes+'_kev_'+species+'_new_flux'
    options, vars, 'ytitle', ' '
    options, vars, 'ytickformat', '(A1)'
    tplot, vars, trange=plot_time_range, position=right_poss, /noerase, /novtitle
    timebar, plot_time, linestyle=1

    ; Add figure label.
    for ii=0, nprobe-1 do begin
        tpos = right_poss[*,ii]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,/normal, right_label+'-'+string(ii+1,format='(I0)')+'.'
    endfor

    ; Add title.
    tpos = right_poss[*,0]
    tx = total(tpos[[0,2]])*0.5
    ty = tpos[3]+ychsz*0.5
    xyouts, tx,ty,/normal,alignment=0.5, 'Traced back to MLT of dispersionless injection'

    ; Add new MLT.
    for ii=0, nprobe-1 do begin
        tpos = right_poss[*,ii]
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,/normal, alignment=1, string(new_mlts[ii],format='(F5.1)')+' MLT'
    endfor

    ; Add Probe.
    for ii=0, nprobe-1 do begin
        tpos = left_poss[*,ii]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.2
        xyouts, tx,ty,/normal,alignment=0, strupcase(probes[ii]), color=probe_colors[ii]

        tpos = right_poss[*,ii]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.2
        xyouts, tx,ty,/normal,alignment=0, strupcase(probes[ii]), color=probe_colors[ii]
    endfor

    if keyword_set(test) then stop
    sgclose




end
