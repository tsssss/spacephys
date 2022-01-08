;+
; Test wavelet with different padding.
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


;---Input.
    day = time_double('2014-08-28')
    probe = 'b'

;---Settings.
    day_time_range = day+[0,86400d]
    time_range = day_time_range+[-1,1]*86400d
    pad_times = [0,300,3600d]
    prefix = 'rbsp'+probe+'_'
    ndim = 3
    uvw = constant('uvw')
    xyz = constant('xyz')
    fillval = !values.f_nan
    common_time_step = 1d
    common_times = make_bins(time_range, common_time_step)
    ncommon_time = n_elements(common_times)

;---Read data.
    b_uvw_var = prefix+'b_uvw'
    if check_if_update(b_uvw_var, time_range, dtime=60.) then $
        rbsp_read_emfisis, time_range, probe=probe, id='l2%magnetometer'
    interp_time, b_uvw_var, common_times
    add_setting, b_uvw_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'nT', $
        'short_name', 'B', $
        'coord', 'UVW', $
        'coord_labels', uvw )
    b_uvw = get_var_data(b_uvw_var)

    ; Load the fixed q_uvw2gse.
    rbsp_read_q_uvw2gse, time_range, probe=probe

    ; Mask invalid data.
    index = where(b_uvw[*,2] le -99999, count)
    if count ne 0 then b_uvw[index,*] = fillval
    b_gse_var = prefix+'b_gse'
    b_gse = cotran(b_uvw, common_times, 'uvw2gse', probe=probe)
    store_data, b_gse_var, common_times, b_gse
    add_setting, b_gse_var, /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', '#', $
        'coord', 'GSE', $
        'coord_labels', xyz )

    ; Remove spin tone using wavelet.
    spin_period = 12.
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


    b_gse = get_var_data(b_gse_var)
    foreach component, xyz, ii do begin
        comp_var = prefix+'b'+xyz[ii]+'_gse'
        foreach pad_time, pad_times do begin
            index = lazy_where(common_times, '[]', day_time_range+[-1,1]*pad_time)
            store_data, comp_var, common_times[index], b_gse[index,ii]

        ;---Calculate cwt.
            cwt_var = comp_var+'_cwt'
            if check_if_update(comp_var+'_wps', time_range) then begin
                tic
                calc_psd, comp_var, scales=target_scales
                toc
                print, pad_time
            endif

        ;---Trim to day_time_range.
            get_data, cwt_var, 0, cwt
            data = get_var_data(comp_var, limit=lim, times=cwt_times)
            data2 = data-wavelet_reconstruct(cwt, index=bad_index)
            the_var = comp_var+'_'+string(pad_time,format='(I0)')+'_sec'
            store_data, the_var, cwt_times, data2
        endforeach

        vars = comp_var+'_'+string(pad_times,format='(I0)')+'_sec'
        the_time_range = day_time_range+[1,-1]*spin_period*2*3
        times = common_times[lazy_where(common_times, '[]', the_time_range)]
        data0 = get_var_data(vars[0], at=times)
        for jj=1,n_elements(vars)-1 do begin
            ddata = get_var_data(vars[jj], at=times)-data0
            store_data, comp_var+'_diff'+string(jj,format='(I0)'), times, ddata
        endfor
    endforeach

end