;+
; Test to trace the injection location in Turner+2017.
;-


    ; Load Themis electron data.
    time_range = time_double(['2016-04-07/01:15','2016-04-07/02:15'])
    probe = 'd'
    energy_range = [90,300]
    themis_read_kev_electron, time_range, probe=probe, energy=energy_range
    themis_read_orbit, time_range, probe=probe
    
    stop
    
    prefix = 'th'+probe+'_'
    the_var = prefix+'kev_e_flux'
    get_data, the_var, times, fluxes, energies, limits=lim

    
    r_var = prefix+'r_gsm'
    model = 'igrf'
    species = 'e'
    pitch_angle0 = 90.
    time0 = mean(time_range)
    r_gsm0 = get_var_data(r_var, at=time0)
    
    
;---Calculate drift period.
    v_drifts = list()
    foreach energy, energies do begin
        v_drift = calc_drift_velocity(rgsm=r_gsm0, time=time0, $
            energy=energy, pitch_angle=pitch_angle0, species=species, $
            model=model, bounce_period=bounce_period)
        print, energy, v_drift
        v_drifts.add, v_drift
    endforeach
    v_drifts = v_drifts.toarray()
    
    drift_periods = 360/v_drifts/60 ; in min.
    xxs = make_bins(minmax(energies),10)
    yys = interpol(drift_periods,energies, xxs, spline=1)
    
    plot_file = 0
    sgopen, plot_file, xsize=6, ysize=6, /inch
    tpos = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz)
    plot, xxs, yys, position=tpos, $
        xtitle='Ele energy (keV)', ytitle='Drift period (min)'
    plots, energies, drift_periods, psym=1, color=sgcolor('blue')
    xyouts, tpos[2]-xchsz*1, tpos[3]-ychsz*1, /normal, alignment=1, 'Sheng', color=sgcolor('blue')
    
    
    test_energies = [241,165,113]
    test_drift_periods = [19.7,28.9,41.3]
    test_times = time_double('2016-04-07/'+['01:35:20','01:44:15','01:50:25'])
    plots, test_energies, test_drift_periods, psym=1, color=sgcolor('red')
    xyouts, tpos[2]-xchsz*1, tpos[3]-ychsz*2, /normal, alignment=1, 'Turner', color=sgcolor('red')

    

end