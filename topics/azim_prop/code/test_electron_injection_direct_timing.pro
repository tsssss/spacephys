;+
; Plot keV electron.
;-

;pro test_electron_injection_driect_timing

;---Load electron flux.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 0

;---Settings.
    einfo = get_var_data(info_var)
    time_range = time_double(['2014-08-28/10:00','2014-08-28/10:20'])
    time_range_plot = time_double(['2014-08-28/09:50','2014-08-28/11:00'])


    probes = ['thd','g15']
    nprobe = n_elements(probes)
    the_energy = 40
    data_rate = 3
    times = make_bins(time_range, data_rate)
    ntime = n_elements(times)
    fluxes = dblarr(ntime,nprobe)
    the_energys = intarr(nprobe)
    foreach probe, probes, ii do begin
        get_data, probe+'_kev_e_flux', uts, flux, energy_bins
        tmp = min(energy_bins-the_energy, /absolute, index)
        fluxes[*,ii] = alog10(interpol(flux[*,index], uts, times))
        the_energys[ii] = energy_bins[index]
    endforeach

    time_shift_range = [-1,1]*5*60
    record_shifts = make_bins(time_shift_range/data_rate, 1)
    time_shifts = record_shifts*data_rate    
    corrs = c_correlate(fluxes[*,0],fluxes[*,1], record_shifts)
    max_corr = max(corrs, index)
    time_lag = time_shifts[index]
    record_lag = record_shifts[index]
    
    times = make_bins(time_range_plot, data_rate)
    ntime = n_elements(times)
    fluxes = dblarr(ntime,nprobe)
    foreach probe, probes, ii do begin
        get_data, probe+'_kev_e_flux', uts, flux, energy_bins
        tmp = min(energy_bins-the_energy, /absolute, index)
        fluxes[*,ii] = interpol(flux[*,index], uts, times)
    endforeach
    fluxes[*,1] = shift(fluxes[*,1], -record_lag)
    
    store_data, 'test_direct_shift', times, [[fluxes[*,0]],[shift(fluxes[*,1],-index)]], limits={$
        ytitle:'(#/cm!U2!N-s-sr-keV)', colors:sgcolor(['black','red']), $
        labflag: 1, ylog:1, $
        labels:[strupcase(probes[0])+' '+sgnum2str(the_energys[0])+' keV', strupcase(probes[1])+' '+sgnum2str(the_energys[1])+' keV']}



;---Plot.
    fig_xsize = 4
    fig_ysize = 4

    if keyword_set(test) then begin
        file = test
        magnify = 1.2
    endif else begin
        file = join_path([sparentdir(srootdir()),'plot','test_electron_injection_direct_timing.pdf'])
        magnify = 1
    endelse

    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify

    poss = sgcalcpos(2, margins=[9,4,10,2], xchsz=xchsz, ychsz=ychsz, ypad=5)
    label_size = 0.8
    
    tpos = poss[*,0]
    yrange = [0,1]
    xrange = time_shift_range
    plot, time_shifts, corrs, /noerase, position=tpos, $
        xstyle=1, xrange=xrange, xticks=4, xminor=5, xtitle='Time shift (sec)', xticklen=-0.02, $
        ystyle=1, yrange=yrange, yticks=2, yminor=2, ytitle='Cross Correlation (#)', yticklen=-0.01
    plots, time_lag+[0,0], yrange, linestyle=1
    tx = tpos[0]+xchsz*0
    ty = tpos[3]+ychsz*0.5
    xyouts, tx,ty,/normal, 'a. Cross correlation'
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    xyouts, tx, ty, /normal, 'Time lag (sec): '+sgnum2str(time_lag,ndec=0), charsize=label_size
    ty = tpos[3]-ychsz*2
    xyouts, tx, ty, /normal, 'Max correlation: '+sgnum2str(max_corr,ndec=1), charsize=label_size
    
    tpos = poss[*,1]
    tplot, 'test_direct_shift', trange=time_range_plot+[1,-1]*max(time_shifts), /novtitle, position=tpos, /noerase
    tx = tpos[0]+xchsz*0
    ty = tpos[3]+ychsz*0.5
    xyouts, tx,ty,/normal, 'b. After time shift'

    if keyword_set(test) then stop
    sgclose

end
