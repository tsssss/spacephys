;+
; Test to trace keV electron backward for selected energy bins.
;-

;---Load electron flux.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 0

;---Settings.
    einfo = get_var_data(info_var)
    time_range_plot = einfo.injection['time_range']
    model = 'igrf'
    species = 'e'


    pre0 = '1991-080_'
    energy_range = [50.,300]    ; keV.
    ref_times = time_double([$
        '2014-08-28/10:09:08', $
        '2014-08-28/10:07:48', $
        '2014-08-28/10:07:18', $
        '2014-08-28/10:06:58', $
        '2014-08-28/10:06:58'])
    
;    pre0 = 'g15_'
;    energy_range = [40,300]     ; keV.
;    ref_times = time_double([$
;        '2014-08-28/10:08:02', $
;        '2014-08-28/10:07:39', $
;        '2014-08-28/10:07:27', $
;        '2014-08-28/10:07:23'])

    pre0 = 'g13_'
    energy_range = [30,300]     ; keV.
    ref_times = time_double([$
        '2014-08-28/10:16:04', $
        '2014-08-28/10:13:00', $
        '2014-08-28/10:09:55', $
        '2014-08-28/10:08:40'])
    ref_times = time_double([$
        '2014-08-28/10:18:34', $
        '2014-08-28/10:13:46', $
        '2014-08-28/10:10:10', $
        '2014-08-28/10:08:44'])
        
    pre0 = 'rbspb_'
    energy_range = [40.,200]    ; keV.
    ref_times = time_double([$
        '2014-08-28/10:13:38', $
        '2014-08-28/10:11:49', $
        '2014-08-28/10:10:23', $
        '2014-08-28/10:09:28', $
        '2014-08-28/10:08:56'])


    fvar0 = pre0+'kev_'+species+'_flux'    ; electron flux.
    rvar = pre0+'r_gsm'


;---Filter using energy.
    get_data, fvar0, times, flux, energy_bins, limits=lim
    index = lazy_where(energy_bins, energy_range)

    fvar = pre0+'_flux'
    store_data, fvar, times, flux[*,index], energy_bins[index], limits=lim
    options, fvar, 'labels', lim.labels[index]
    options, fvar, 'colors', lim.colors[index]
    yrange = alog10(minmax(flux[*,index]))
    yrange = 10.^[floor(yrange[0]),ceil(yrange[1])]
    options, fvar, 'yrange', yrange

    tplot, fvar, trange=time_range_plot
    stop


;---Get the velocity for each energy bin.
    get_data, fvar, times, flux, energy_bins
    nenergy_bin = n_elements(energy_bins)
    drift_vels = dblarr(nenergy_bin)
    foreach energy_bin, energy_bins, ii do begin
        time = ref_times[ii]
        rgsm0 = get_var_data(rvar, at=time)
        drift_vels[ii] = calc_drift_velocity(rgsm=rgsm0, time=time, $
            energy=energy_bin, pitch_angle=90, species=species, $
            model=model, bounce_period=bounce_period)
        print, 'Energy (keV): ', energy_bin
        print, 'Drift velocity (km/s): ', drift_vels[ii]
    endforeach
    stop



;---Get the flux peak for given time period.
    ut_start = time_double('2014-08-28/09:30')
    ut_ref = time_double('2014-08-28/10:00')
    colors = get_setting(fvar, 'colors')
    xrange = [-2000,1000]
    yrange = [-6,0]
    plot, xrange, yrange, /nodata, $
        xstyle=1, xtitle='UT (sec)', $
        ystyle=1, ytitle='MLT (hr)'
    for i=0, nenergy_bin-1 do begin
        txs = [ut_start,ref_times[i]]-ut_ref
        tys = (txs-txs[1])*drift_vels[i]/15
        oplot, txs, tys, color=colors[i]
    endfor
    stop

end
