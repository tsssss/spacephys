;+
; Test L3 PA.
;-

function test_tracers_l3_pa

    compile_opt idl2

    ; Sergei's event.
    tr = ['2026-04-01/03:08','2026-04-01/03:10']
    tr = ['2026-04-01/03:08:50','2026-04-01/03:09:20']
    tr = time_double(tr)
    probe = '2'
    energy_range = [1e2,1e4]
    prefix = 'ts'+probe+'_'

;---L2.
    files = tracers_load_ace(tr, probe=probe, id='cdaweb%l2')
    prefix2 = prefix+'l2_ace_'
    anode_var = prefix2+'TSCS_anode_angle'
    anode_angles = cdf_read_var(anode_var, filename=files[0])
    flux_var = prefix2+'def'
    fluxs = cdf_read_var(flux_var, filename=files)
    energy_var = prefix2+'energy'
    energy_bins = cdf_read_var(energy_var, filename=files[0])
    epoch_var = 'Epoch'
    epochs = cdf_read_var(epoch_var, filename=files)
    times = convert_time(epochs, from='tt2000', to='unix')

    pa_spec_var = prefix2+'pa_spec_l2'
    en_index = where_pro(energy_bins,'[]',energy_range)
    pa_spec = mean(fluxs[*,en_index,*], dimension=2, nan=1)
    anode_angles = 180-anode_angles
    unit = 'eV/cm!U2!N-s-sr-eV'
    settings = dictionary($
        'display_type', 'spec', $
        'ylog', 0, $
        'yrange', [-10,190], $
        'ytickv', [0,45,90,135,180], $
        'yticks', 4, $
        'yminor', 5, $
        'ytitle', 'PA!C(deg)', $
        'zlog', 1, $
        'zrange', [1e6,1e9], $
        'ztitle', '('+unit+')' )
    pa_l2_var = var_store(pa_spec_var, pa_spec, times, anode_angles, settings=settings)


;---My version.
    ; Load the files.
    files = tracers_load_mag(tr, probe=probe, id='l1b%bdc-roi_x232')
    prefix2 = prefix+'l1b_'
    b_var = prefix2+'bdc_roi'
    epoch_var = 'Epoch'
    epochs = cdf_read_var(epoch_var, filename=files)
    nepoch = n_elements(epochs)
    del_epochs = cdf_read_var('EpochOffset', filename=files[0])
    ndel_epoch = n_elements(del_epochs)
    b_ts_mag = cdf_read_var(b_var, filename=files)
    ntime = nepoch*ndel_epoch
    ndim = 3
    b_ts_mags = dblarr(ntime,ndim)
    raw_times = convert_time(epochs, from='tt2000', to='unix')
    del_times = del_epochs*1e-9   ; ns to sec.
    times = dblarr(ntime)
    for ii=0,nepoch-1 do begin
        i0 = ii*ndel_epoch
        i1 = (ii+1)*ndel_epoch
        b_ts_mags[i0:i1-1,*] = transpose(reform(b_ts_mag[ii,*,*]))
        times[i0:i1-1] = raw_times[ii]+del_times
    endfor

    b_tscs = cotran_pro(b_ts_mags, probe=probe, coord_msg='ts_'+['mag','tscs'])
    b_tscs_var = prefix+'b_tscs'
    settings = dictionary('coord', 'ts_tscs', 'id', 'bfield')
    b_tscs_var = var_store(b_tscs_var, b_tscs, times, settings=settings)

    ; Calculate pitch angles, in [0,360] deg.
    pa_range = [0,180d]
    pa_step = 10d
    pa_bins = smkarthm(pa_range[0],pa_range[1],pa_step,'dx')
    npa_bin = n_elements(pa_bins)


    ; Reload data, just to be sure.
    files = tracers_load_ace(tr, probe=probe, id='cdaweb%l2')
    prefix2 = prefix+'l2_ace_'
    anode_var = prefix2+'TSCS_anode_angle'
    anode_angles = cdf_read_var(anode_var, filename=files[0])
    flux_var = prefix2+'def'
    full_fluxs = cdf_read_var(flux_var, filename=files)
    energy_var = prefix2+'energy'
    energy_bins = cdf_read_var(energy_var, filename=files[0])
    epoch_var = 'Epoch'
    epochs = cdf_read_var(epoch_var, filename=files)
    full_times = convert_time(epochs, from='tt2000', to='unix')
    index = where_pro(full_times, '[]', tr, count=ntime)
    times = full_times[index]
    fluxs = full_fluxs[index,*,*]
    pa_spec = fltarr(ntime, npa_bin)
    

    ; Annode angles. 0 deg is along the +Z axis, 90 deg is along the +Y axis.
    nanode = n_elements(anode_angles)
    anode_dirs = fltarr(nanode,ndim)
    anode_dirs[*,0] = 0
    tts = anode_angles*constant('rad')
    anode_dirs[*,1] = sin(tts)
    anode_dirs[*,2] = cos(tts)

    foreach time, times, tid do begin
        the_b_tscs = var_get_data(b_tscs_var, at=time)
        ; expand to [n,3]
        b_dirs = replicate(1.,nanode) # transpose(the_b_tscs)
        the_pas = sang(anode_dirs, b_dirs, degree=1)
        the_fluxs = mean(reform(fluxs[tid,en_index,*]),nan=1,dimension=1)
        rounded_pas = round(the_pas/pa_step)*pa_step
        uniq_pas = suniq(rounded_pas)
        foreach pa, uniq_pas do begin
            index = where(rounded_pas eq pa)
            pa_index = where(pa_bins eq pa)
            pa_spec[tid,pa_index] = mean(the_fluxs[index],nan=1)
        endforeach
    endforeach


    unit = 'eV/cm!U2!N-s-sr-eV'
    pa_ts1_var = prefix+'pa_spec_ts1'
    settings = dictionary($
        'display_type', 'spec', $
        'ylog', 0, $
        'yrange', [0,180], $
        'ytickv', [0,45,90,135,180], $
        'yticks', 4, $
        'yminor', 5, $
        'ytitle', 'PA!C(deg)', $
        'zlog', 1, $
        'zrange', [1e6,1e9], $
        'ztitle', '('+unit+')' )
    pa_ts1_var = var_store(pa_ts1_var, pa_spec, times, pa_bins, settings=settings)



;---My version2.
    ; Calculate pitch angles, in [0,360] deg.
    pa_range = [0,360d]
    pa_step = 10d
    pa_bins = smkarthm(pa_range[0],pa_range[1],pa_step,'dx')
    npa_bin = n_elements(pa_bins)
    pa_spec = fltarr(ntime, npa_bin)

    spin_period = 3d
    foreach time, times, tid do begin
        the_b_tscs = var_get_data(b_tscs_var, at=time)
        the_b_tscs[[0:1]] = 0
        if the_b_tscs[2] lt 0 then begin
            the_pas = 180-anode_angles
        endif else begin
            the_pas = anode_angles
        endelse
        the_fluxs = mean(reform(fluxs[tid,en_index,*]),nan=1,dimension=1)
        rounded_pas = round(the_pas/pa_step)*pa_step
        uniq_pas = suniq(rounded_pas)
        foreach pa, uniq_pas do begin
            index = where(rounded_pas eq pa)
            pa_index = where(pa_bins eq pa)
            pa_spec[tid,pa_index] = mean(the_fluxs[index],nan=1)
        endforeach

        ; Now find the next half.
        the_time = time+spin_period*0.5
        the_b_tscs = var_get_data(b_tscs_var, at=the_time)
        the_b_tscs[[0:1]] = 0
        if the_b_tscs[2] lt 0 then begin
            the_pas = 180-anode_angles
        endif else begin
            the_pas = anode_angles
        endelse
        tmp = min(full_times-the_time, absolute=1, full_tid)
        the_fluxs = mean(reform(full_fluxs[full_tid,en_index,*]),nan=1,dimension=1)
        rounded_pas = round(the_pas/pa_step)*pa_step
        uniq_pas = suniq(rounded_pas)
        foreach pa, uniq_pas do begin
            index = where(rounded_pas eq pa)
            pa_index = where(pa_bins eq pa)
            pa_spec[tid,pa_index] = mean(the_fluxs[index],nan=1)
        endforeach
    endforeach

    ; Convert from [0,360] to [-90,270] deg.
    yrange = [-90,270]
    ytickv = make_bins(yrange, 90)
    yticks = n_elements(ytickv)-1
    index = where(pa_bins ge max(yrange), count)
    if count ne 0 then begin
        pa_bins[index] = pa_bins[index]-360
    endif
    index = sort(pa_bins)
    pa_bins = pa_bins[index]
    pa_spec = pa_spec[*,index]

    pa_ts2_var = prefix+'pa_spec_ts2'
    settings = dictionary($
        'display_type', 'spec', $
        'ylog', 0, $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', 5, $
        'ytitle', 'PA!C(deg)', $
        'zlog', 1, $
        'zrange', [1e6,1e9], $
        'ztitle', '('+unit+')' )
    pa_ts2_var = var_store(pa_ts2_var, pa_spec, times, pa_bins, settings=settings)





;---L3.
    files = tracers_load_ace(tr, probe=probe, id='iowa%l3')

    prefix = 'ts'+probe+'_'
    prefix2 = prefix+'l3_ace_'
    pa_var = prefix2+'pitch_angle'
    en_var = prefix2+'energy'
    flux_var = prefix2+'pitch_def'
    fluxs = cdf_read_var(flux_var, filename=files)
    energy_bins = cdf_read_var(en_var, filename=files)
    pa_bins = cdf_read_var(pa_var, filename=files)
    epoch_var = 'Epoch'
    epochs = cdf_read_var(epoch_var, filename=files)
    times = convert_time(epochs, from='tt2000', to='unix')

    unit = 'eV/cm!U2!N-s-sr-eV'
    pa_spec_var = prefix+'pa_spec_l3'
    en_index = where_pro(energy_bins,'[]',energy_range)
    pa_spec = mean(fluxs[*,en_index,*], dimension=2, nan=1)
    settings = dictionary($
        'display_type', 'spec', $
        'ylog', 0, $
        'yrange', [0,180], $
        'ytickv', [0,45,90,135,180], $
        'yticks', 4, $
        'yminor', 5, $
        'ytitle', 'PA!C(deg)', $
        'zlog', 1, $
        'zrange', [1e6,1e9], $
        'ztitle', '('+unit+')' )
    pa_l3_var = var_store(pa_spec_var, pa_spec, times, pa_bins, settings=settings)

    en_spec_var = prefix+'en_spec_l3'
    en_spec = mean(fluxs, dimension=3, nan=1)
    settings = dictionary($
        'display_type', 'spec', $
        'ylog', 1, $
        'ytitle', 'Energy!C(eV)', $
        'zlog', 1, $
        'zrange', [1e6,1e9], $
        'ztitle', '('+unit+')' )
    en_spec_var = var_store(en_spec_var, en_spec, times, energy_bins, settings=settings)


    plot_vars = [en_spec_var,pa_l2_var,pa_l3_var,pa_ts1_var,pa_ts2_var]
    nplot_var = n_elements(plot_vars)
    ypans = replicate(1.,nplot_var)
    pid = where(plot_vars eq pa_ts2_var, count)
    if count ne 0 then ypans[pid] = 2
    poss = sgcalcpos(nplot_var, ypans=ypans)

    tplot, plot_vars, trange=tr, position=poss
    stop

end

compile_opt idl2
print, test_tracers_l3_pa()
end