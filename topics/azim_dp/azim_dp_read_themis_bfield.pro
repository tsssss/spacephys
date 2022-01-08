;+
; Read THEMIS B field for a given time and probe.
;-

pro azim_dp_read_themis_bfield, time_range, probe=probe, errmsg=errmsg


;---Load 3 sec data.
    common_time_step = 3.0
    themis_read_fgm, time_range, id='l2%fgs', probe=probe, errmsg=errmsg
    if errmsg ne '' then return

    prefix = 'th'+probe+'_'
    b_var = prefix+'b_gsm'
    flag_var = prefix+'fgm_fgs_quality'
    rename_var, prefix+'fgs_gsm', to=b_var
    uniform_time, b_var, common_time_step
    uniform_time, flag_var, common_time_step
    add_setting, b_var, /smart, {$
        display_type: 'vector', $
        unit: 'nT', $
        short_name: 'B', $
        coord: 'GSM', $
        coord_labels: ['x','y','z']}

;---Flag for bad data.
    ; Looks like: 2 for eclipse, 1 for commisional phase.
    pad_time = 120.  ; sec.
    get_data, flag_var, times, flags
    ntime = n_elements(times)
    all_flags = bytarr(ntime)+1
    index = where(flags ge 2, count)
    if count ne 0 then all_flags[index] = 0

    index = where(all_flags eq 0, count)
    if count ne 0 then begin
        get_data, b_var, times, b_gsm, limits=lim
        store_data, b_var+'_before', times, b_gsm, limits=lim
        bad_times = times[time_to_range(index,time_step=1)]
        bad_times[*,0] -= pad_time
        bad_times[*,1] += pad_time
        nbad_time = n_elements(bad_times)/2
        for ii=0, nbad_time-1 do b_gsm[lazy_where(times,'[]',reform(bad_times[ii,*])),*] = !values.f_nan

        ; Check out L1 data with eclipse correction.
        thm_load_fgm, probe=probe, level=1, type='calibrated', use_eclipse_corrections=1, trange=time_range
        dsl2gse, prefix+'fgl', prefix+'state_spinras', prefix+'state_spindec', prefix+'l1_b_gse'
        get_data, prefix+'l1_b_gse', l1_times, b_gse
        l1_b_gsm = cotran(b_gse, l1_times, 'gse2gsm')
        for ii=0, nbad_time-1 do begin
            index = lazy_where(l1_times, '[]', bad_times[ii,*], count=count)
            if count eq 0 then continue
            the_times = l1_times[index]
            the_b_gsm = l1_b_gsm[index,*]
            b1_time_range = the_times[[0,count-1]]
            index = lazy_where(times, '[]', b1_time_range, count=count)
            if count eq 0 then continue
            b_gsm[index,*] = sinterpol(the_b_gsm, the_times, times[index])
        endfor
        
        store_data, b_var, times, b_gsm
    endif

end

time_range = time_double(['2019-08-04/20:00','2019-08-05/08:00'])
probe = 'e'
azim_dp_read_themis_bfield, time_range, probe=probe
end
