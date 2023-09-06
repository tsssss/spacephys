;+
; Study the local time of electron injection as a function of time.
;-


;---Load electron flux.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 0

;---Settings.
    einfo = get_var_data(info_var)
    time_range_plot = time_double(['2014-08-28/09:50','2014-08-28/11:00'])

    probes = ['LANL-02A','LANL-04A','LANL-97A','1994-084','rbspb','g13','1991-080','tha','the','g15','thd','LANL-01A']
    nprobe = n_elements(probes)
    color_start = 50
    color_end = 250
    color_table = 40
    colors = round(smkarthm(color_start,color_end, nprobe, 'n'))
    for ii=0, nprobe-1 do colors[ii] = sgcolor(colors[ii], ct=color_table)
    colors[where(probes eq 'the')] = sgcolor('gold')    ; the yellow is too bright.


    the_probes = ['thd','g15','1991-080','g13','rbspb','1994-084']
    plot_info = dictionary()
    plot_info['probes'] = the_probes
    plot_info['colors'] = []
    plot_info['energy_range'] = [50,500]
    foreach probe, the_probes do begin
        plot_info.colors = [plot_info.colors,colors[where(probes eq probe)]]
    endforeach

    foreach probe, plot_info.probes do begin
        get_data, probe+'_kev_e_flux', times, flux, energy_bins
        index = where_pro(energy_bins, plot_info.energy_range, count=nenergy_bin)
        flux = flux[*,index]
        energy_bins = energy_bins[index]

;        if nenergy_bin/2 ge 4 then begin
;            flux = flux[*,0:*:2]
;            energy_bins = energy_bins[0:*:2]
;        endif

        yrange = alog10(minmax(flux))
        yrange = 10.^[floor(yrange[0]),ceil(yrange[1])]

        var = probe+'_kev_e'
        store_data, var, times, flux, energy_bins
        add_setting, var, /smart, {$
            display_type: 'list', $
            ylog: 1, $
            yrange: yrange, $
            color_table: 52, $
            unit: '#/cm!U2!N-s-sr-keV', $
            value_unit: 'keV', $
            short_name: 'e!U-!N flux'}
    endforeach
    
    stop
    
    
    probes = ['thd','g15','1991-080','g13','rbspb']
    nprobe = n_elements(probes)
    mlts = list()
    times = list()
    foreach probe, probes do begin
        if probe eq 'g13' then begin
            mlts.add, [4.2,4.4]
            times.add, time_double('2014-08-28/'+['10:15:10','10:15:30'])
        endif else begin
            case probe of
                'thd': the_time = '10:08:30'
                'g15': the_time = '10:10:30'
                '1991-080': the_time = '10:13:20'
                'rbspb': the_time = '10:10:00'
            endcase
            the_time = time_double('2014-08-28/'+the_time)
            times.add, the_time
            r_gsm = get_var_data(probe+'_r_gsm', at=the_time)
            r_sm = cotran(r_gsm, the_time, 'gsm2sm')
            the_mlt = pseudo_mlt(r_sm)
            mlts.add, the_mlt
        endelse
    endforeach
stop    
    mlt_range = [0,6]
    time_range = time_double('2014-08-28/'+['10:00','10:30'])
    xrange = time_range-time_range[0]
    yrange = [0,6]
    plot, xrange, yrange, /nodata, xstyle=1, ystyle=1
    foreach probe, probes, ii do plots, times[ii]-time_range[0], mlts[ii], psym=-1
    stop
    
    
;---G13.
    info = dictionary($
        'probe', 'g13', $
        'energy_range', [50.,500], $
        ;'time_range', time_double('2014-08-28/'+['10:15','10:23']))
        'time_range', time_double('2014-08-28/'+['10:15','10:20']))
    pos_var = info.probe+'_r_gsm'
    eflux_var = info.probe+'_kev_e'
    get_data, eflux_var, times, efluxes, energies
    index = where_pro(energies,'[]', info.energy_range)
    efluxes = efluxes[*,index]
    energies = energies[index]
    index = where_pro(times,'[]', info.time_range)
    times = times[index]
    efluxes = efluxes[index,*]
    
    nenergy = n_elements(energies)-1
    the_energies = sqrt(energies[1:nenergy]*energies[0:nenergy-1])
    the_times = dblarr(nenergy)
    the_mlts = fltarr(nenergy)
    drift_periods = fltarr(nenergy)
    foreach the_energy, the_energies, ii do begin
        the_eflux = efluxes[*,where(energies eq max(energies[ii:ii+1]))]
        max_eflux = min(the_eflux, index)
        the_time = times[index]
        the_times[ii] = the_time
        r_gsm = get_var_data(pos_var, at=the_time)
        r_sm = cotran(r_gsm, the_time, 'gsm2sm')
        the_mlts[ii] = pseudo_mlt(r_sm)
        drift_periods[ii] = calc_drift_period(rgsm=r_gsm, time=the_time, energy=the_energy, species='e')
    endforeach
    info['drift_periods'] = drift_periods
    info['times'] = the_times
    info['energies'] = the_energies
    info['mlts'] = the_mlts

    ;info['times'] = time_double('2014-08-28/'+['10:21:52','10:19:29','10:18:09'])
    tpos = sgcalcpos(1)
    yrange = [0,6]
    store_data, 'test', info.times, info.mlts, limits={psym:1, yrange:yrange}
    tplot, 'test', trange=info.time_range, position=tpos
    plot, info.time_range, yrange, /nodata, /noerase, xstyle=5, ystyle=5, position=tpos
    foreach drift_period, drift_periods, ii do begin
        xxs = [info.time_range[0],info.times[ii]]
        slope = (24d/drift_period)
        yys = (xxs-info.times[ii])*slope+info.mlts[ii]
        oplot, xxs, yys
    endforeach
    
    ; The injection originated from 10:16:09 to 10:16:34, in 2.9 to 3.1 MLT.
    
    
    stop



;---Plot.
    vars = plot_info.probes+'_kev_e'
    nvar = n_elements(vars)
    fig_xsize = 6
    aspect_ratio = 0.25
    fig_ysize = nvar*fig_xsize*aspect_ratio

    test = 1
    if keyword_set(test) then begin
        file = test
        magnify = 1.2
    endif else begin
        file = join_path([sparentdir(srootdir()),'plot','fig_case2_kev_e_flux.pdf'])
        magnify = 1
    endelse

    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify

    poss = sgcalcpos(nvar, margins=[9,4,8,2], xchsz=xchsz, ychsz=ychsz)
    tplot, vars, position=poss, /novtitle, trange=time_range_plot

    fig_labels = ['a','b','c','d','e','f','g']
    full_ysz = 0.8
    foreach probe, plot_info.probes, ii do begin
        tpos = poss[*,ii]
        tx = tpos[0]+xchsz*1
        ty = tpos[3]-ychsz*1
        label = fig_labels[ii]+'.  '
        xyouts, tx,ty, /normal, label+strupcase(probe), color=(plot_info.colors)[ii]
        xyouts, tx,ty, /normal, label
    endforeach


    if keyword_set(test) then stop
    sgclose

end