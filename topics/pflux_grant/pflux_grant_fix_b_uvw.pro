;+
; Remove wobble in B.
;-

function pflux_grant_fix_b_uvw_init_target_scales, spin_period

    spin_freq = 1d/spin_period
    base = 2d^(1d/4)
    freq_range = [0.5,4]/spin_period
    selected_freqs = smkgmtrc(min(freq_range),max(freq_range),base,'dx')
    nselected_freq = n_elements(selected_freqs)
    selected_periods = 1d/selected_freqs
    wavelet_info = wavelet_info()
    t2s = wavelet_info[7]
    selected_scales = selected_periods*t2s

    return, dictionary($
        'base', base, $                 ; separation between terms.
        'spin_period', spin_period, $   ; spin period in sec.
        'spin_freq', spin_freq, $       ; spin freq in Hz.
        'nscale', nselected_freq, $     ; # of target scales.
        'scales', selected_scales, $    ; scales in sec.
        'freqs', selected_freqs, $      ; freqs in Hz.
        'periods', selected_periods, $  ; periods in sec.
        't2s', t2s, $                   ; convert period to scale.
        's2t', 1d/t2s)                  ; convert scale to period.
end

pro pflux_grant_fix_b_uvw, time_range, probe=probe, test=test, pad_suffix=pad_suffix

    prefix = 'rbsp'+probe+'_'
    b_uvw_var = prefix+'b_uvw'
    if check_if_update(b_uvw_var) then stop

    ; Mask invalid data with NaN.
    cal_state = get_var_data(prefix+'cal_state', times=uts)
    mag_valid = get_var_data(prefix+'mag_valid')
    bad_index = where(cal_state ne 0 or mag_valid eq 1, count, complement=good_index)
    fillval = !values.f_nan
    pad_time = 5.   ; sec.
    if count ne 0 then begin
        time_ranges = uts[time_to_range(bad_index,time_step=1)]
        ntime_range = n_elements(time_ranges)*0.5
        b_uvw = get_var_data(b_uvw_var, times=uts)
        for ii=0,ntime_range-1 do begin
            index = where_pro(uts, '[]', time_ranges[ii,*]+[-1,1]*pad_time, count=count)
            if count eq 0 then continue
            b_uvw[index,*] = fillval
        endfor
        store_data, b_uvw_var, uts, b_uvw
    endif


;---Read data.
    ndim = 3
    uvw = constant('uvw')
    xyz = constant('xyz')
    fillval = !values.f_nan
    if n_elements(common_time_step) eq 0 then common_time_step = 1d/16
    if common_time_step gt 1 then message, 'Cadence too low ...'
    common_times = make_bins(time_range, common_time_step)
    ncommon_time = n_elements(common_times)
    data_gap_window = 4*common_time_step
    interp_time, b_uvw_var, common_times, data_gap_window=data_gap_window
    add_setting, b_uvw_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'nT', $
        'short_name', 'B', $
        'coord', 'UVW', $
        'coord_labels', uvw )
    b_uvw = get_var_data(b_uvw_var)


;---Convert to DSC.
    rad = constant('rad')
    spin_phase_var = prefix+'spin_phase'
    rbsp_read_spin_phase, time_range, probe=probe, times=common_times
    spin_phase = get_var_data(spin_phase_var)*rad
    cost = cos(spin_phase)
    sint = sin(spin_phase)
    b_dsc = dblarr(ncommon_time,ndim)
    b_dsc_var = prefix+'b_dsc'
    b_dsc[*,0] = b_uvw[*,0]*cost-b_uvw[*,1]*sint
    b_dsc[*,1] = b_uvw[*,0]*sint+b_uvw[*,1]*cost
    b_dsc[*,2] = b_uvw[*,2]
    store_data, b_dsc_var, common_times, b_dsc
    add_setting, b_dsc_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', '#', $
        'coord', 'DSC', $
        'coord_labels', xyz )


;---Correct in DSC.
    ; Get the background field.
    section_window = 60.
    section_times = make_bins(time_range, section_window)
    nsection_time = n_elements(section_times)-1
    b_dsc_bg = fltarr(nsection_time, ndim)
    b_dsc = get_var_data(b_dsc_var)
    for ii=0,nsection_time-1 do begin
        index = where_pro(common_times, '[]', section_times[ii:ii+1])
        for jj=0,ndim-1 do b_dsc_bg[ii,jj] = median(b_dsc[index,jj])
    endfor
    center_times = section_times[0:nsection_time-1]+section_window*0.5
    b_dsc_bg = sinterpol(b_dsc_bg, center_times, common_times, /quadratic)
    if keyword_set(test) then $
        store_data, prefix+'b_dsc_bg', common_times, b_dsc_bg, limits={colors:constant('rgb')}



    ; Remove spikes and bad data around apogee.
    dbmag = abs(snorm(b_uvw)-snorm(b_dsc_bg))
    ; Normalize with R.
    r_var = prefix+'r_gse'
    rbsp_read_orbit, time_range, probe=probe
    dis = snorm(get_var_data(r_var, times=orbit_times))
    dis = interpol(dis, orbit_times, common_times)
    perigee_shell = 2.5
    dbmag *= (dis/perigee_shell)^1.5
    ; Normalize according to range_flag.
    range_flag = get_var_data(prefix+'range_flag', times=uts)
    range_index = where(range_flag ne 1, range_count)
    if range_count ne 0 then begin
        time_ranges = uts[time_to_range(range_index,time_step=1)]
        ntime_range = n_elements(time_ranges)*0.5
        for ii=0,ntime_range-1 do begin
            index = where_pro(common_times, '[]', time_ranges[ii,*], count=count)
            if count eq 0 then continue
            dbmag[index] *= 0.25
        endfor
    endif
    ; Remove mode switch around apogee.
    range_flag = interpol(range_flag, uts, common_times)
    index = where(range_flag ne 1 and dis ge perigee_shell, count)
    if count ne 0 then begin
        time_ranges = common_times[time_to_range(index,time_step=1)]
        ntime_range = n_elements(time_ranges)*0.5
        for ii=0,ntime_range-1 do begin
            index = where_pro(common_times, '[]', time_ranges[ii,*]+[-1,1]*pad_time, count=count)
            if count eq 0 then continue
            dbmag[index] = fillval
        endfor
    endif
    ; Remove spikes around the perigee mode switch.
    index = where(dbmag ge 10 and dis lt perigee_shell, count)
    if count ne 0 then dbmag[index] = fillval
    ; Remove obvious bad data.
    bmag = snorm(b_uvw)
    index = where(bmag ge 3.4e4, count)
    if count ne 0 then begin
        time_ranges = common_times[time_to_range(index,time_step=1)]
        ntime_range = n_elements(time_ranges)*0.5
        for ii=0,ntime_range-1 do begin
            index = where_pro(common_times, '[]', time_ranges[ii,*]+[-1,1]*300, count=count)
            if count eq 0 then continue
            dbmag[index] = fillval
        endfor
    endif


    ; Mask invalid data.
    index = where(finite(dbmag,/nan), count)
    if count ne 0 then begin
        b_uvw[index,*] = fillval
        b_dsc[index,*] = fillval
    endif
    db_dsc = b_dsc-b_dsc_bg


    ; Now dB should be spike-free, but may contain spin tone.
    foreach component, xyz, ii do store_data, prefix+'db_dsc_'+component, common_times, db_dsc[*,ii]

    ; Remove spin tone using wavelet.
    spin_period = 11.
    target = pflux_grant_fix_b_uvw_init_target_scales(spin_period)
    target_periods = target['periods']
    target_scales = target['scales']
    ntarget_scale = target['nscale']
    bad_index = []
    index_range = make_bins([-1,1]*2,1)
    foreach period, spin_period*[1,0.5] do begin
        tmp = min(target_periods-period, index, /abs)
        bad_index = [bad_index,index_range+index]
    endforeach
    bad_index = sort_uniq(bad_index)


    foreach component, xyz, ii do begin
        comp_var = prefix+'db_dsc_'+component
        data = get_var_data(comp_var, limit=lim)
        uniform_time, comp_var, 0.5
        tdat = get_var_data(comp_var, times=cwt_times)
        tmp = where(finite(tdat,/nan), count, complement=index)
        if count ne 0 then tdat = interpol(tdat[index], cwt_times[index], cwt_times)
        store_data, comp_var, cwt_times, tdat

    ;---Calculate cwt.
        cwt_var = comp_var+'_cwt'
        ;if check_if_update(comp_var+'_wps', time_range) then begin
        calc_psd, comp_var, scales=target_scales
        ;endif

    ;---Remove f and 2f.
        get_data, cwt_var, 0, cwt
        bad_data = wavelet_reconstruct(cwt, index=bad_index)
        db_dsc[*,ii] = data-interp(bad_data, cwt_times, common_times)


        if keyword_set(test) then begin
            comp_var_new = comp_var+'_fixed'
            store_data, comp_var_new, common_times, db_dsc[*,ii], limit=lim
            options, comp_var_new, 'labels', 'B fixed'
            options, comp_var_new, 'colors', 0

            comp_var_combo = comp_var+'_combo'
            store_data, comp_var_combo, common_times, [[data],[db_dsc[*,ii]]], limit=lim
            options, comp_var_combo, 'labels', ['B orig','B fixed']
            options, comp_var_combo, 'colors', sgcolor(['blue','red'])
        endif
    endforeach
    b_dsc = b_dsc_bg+db_dsc


    ; Update the data.
    b_dsc_var = prefix+'b_dsc_fix'
    store_data, b_dsc_var, common_times, b_dsc
    add_setting, b_dsc_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', '#', $
        'coord', 'DSC', $
        'coord_labels', xyz )

    ; Convert back to UVW.
    b_dsc = get_var_data(prefix+'b_dsc_fix')
    b_uvw[*,0] = b_dsc[*,0]*cost+b_dsc[*,1]*sint
    b_uvw[*,1] =-b_dsc[*,0]*sint+b_dsc[*,1]*cost
    b_uvw[*,2] = b_dsc[*,2]
    store_data, b_uvw_var, common_times, b_uvw
    add_setting, b_uvw_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', '#', $
        'coord', 'UVW', $
        'coord_labels', uvw )


    if keyword_set(test) then begin
        two_colors = sgcolor(['blue','red'])
        b_dsc_var = prefix+'b_dsc'
        vec_old = get_var_data(b_dsc_var)
        vec_new = get_var_data(b_dsc_var+'_fix')
        for ii=0,ndim-1 do begin
            the_var = prefix+'b'+xyz[ii]+'_dsc'
            store_data, the_var, common_times, [[vec_old[*,ii]],[vec_new[*,ii]]], $
                limits={colors:two_colors, labels:['orig','fixed'], ytitle:'(nT)'}
        endfor
        tplot_options, 'labflag', -1
        tplot_options, 'ynozero', 1
    endif

    if keyword_set(pad_suffix) then begin
        the_var = prefix+'b_dsc_fix'
        copy_data, the_var, the_var+pad_suffix
    endif

end


probe = 'a'

time = time_double('2012-12-04')    ; gap.
time = time_double('2012-10-02')    ; spikes and gap.
;time = time_double('2012-10-11')    ; spikes.
;time = time_double('2012-10-01')    ; storm.
;time = time_double('2013-11-02')    ; spikes.
;time = time_double('2013-07-15')    ; waves.
;time = time_double('2013-06-24')    ; waves.
;time = time_double('2013-06-10')    ; waves.
;time = time_double('2014-08-28')    ; spin tone.


file = join_path([homedir(),'test.cdf'])
time_range = time+[0,86400d]
rbsp_read_emfisis, time_range, probe=probe, id='l2%magnetometer'
pflux_grant_fix_b_uvw, time_range, probe=probe

stop


;day_time_range = time+[0,86400d]
;xyz = constant('xyz')
;prefix = 'rbsp'+probe+'_'
;
;time_range = day_time_range+[-1,1]*3600
;pflux_grant_fix_b_uvw, time_range, probe=probe, test=1, pad_suffix='_3600_sec'
;time_range = day_time_range+[-1,1]*300
;pflux_grant_fix_b_uvw, time_range, probe=probe, test=1, pad_suffix='_300_sec'
;
;foreach comp_var, prefix+['b_dsc','b_uvw','b_dsc_fix'] do begin
;    foreach comp, xyz, ii do begin
;        data1 = get_var_data(comp_var+'_3600_sec', in=day_time_range, times=times)
;        data2 = get_var_data(comp_var+'_300_sec', at=times)
;        data = [[data1[*,ii]],[data2[*,ii]]]
;        store_data, comp_var+'_'+comp+'_diff', times, data[*,1]-data[*,0]
;    endforeach
;endforeach
;comp_var = prefix+'spin_phase'
;data1 = get_var_data(comp_var+'_3600_sec', in=day_time_range, times=times)
;data2 = get_var_data(comp_var+'_300_sec', at=times)
;data = [[data1[*,0]],[data2[*,0]]]
;store_data, comp_var+'_diff', times, data[*,1]-data[*,0]



end
