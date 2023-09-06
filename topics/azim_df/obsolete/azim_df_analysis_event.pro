;+
; Analysis on event-basis.
;
; The input is the return value of azim_df_find_event.
;-

function azim_df_analysis_event, event, project=project, no_load_data=no_load_data

    if n_elements(project) eq 0 then project = azim_df_load_project()


    event_info_var = 'azim_df_event_info'
    data_file = join_path([project.data_dir,'event_info',$
        'azim_df_event_info_'+strjoin(time_string(event.time_range,tformat='YYYY_MMDD_hhmmss'),'_')+'_v01.tplot'])
    if file_test(data_file) ne 0 then begin
        tplot_restore, filename=data_file
        if tnames(event_info_var) ne '' then begin
            event_info = get_var_data(event_info_var)
            if total(event_info.time_range-event.time_range) eq 0 then return, event_info
        endif
    endif


    ; Load data.
    if ~keyword_set(no_load_data) then azim_df_load_event_data, event.time_range, project=project, event=event

    ; Determine ref_time.
    timing_info = dictionary()
    ; Use DF time.
    ref_time_name = 'df_time'
    ref_times = list()
    foreach df, event.df_list do ref_times.add, df.arrival_time
    timing_info[ref_time_name] = dictionary($
        'ref_times', ref_times.toarray(), $
        'probes', event.probes)
    ; Use xcor time.
    ref_time_name = 'xcor_time'
    ref_times = list()
    foreach df, event.df_list do ref_times.add, df.xcor_time
    timing_info[ref_time_name] = dictionary($
        'ref_times', ref_times.toarray(), $
        'probes', event.probes)
    ; Use combined time.
    ref_time_name = 'combined_time'
    ref_times = list()
    foreach df, event.df_list do ref_times.add, mean([df.arrival_time,df.xcor_time])
    timing_info[ref_time_name] = dictionary($
        'ref_times', ref_times.toarray(), $
        'probes', event.probes)


    ; Determine ref_times, and the associated time_lags.
    ndim = 3
    hr_per_sec_2_deg_per_min = 15.*60
    km_per_s_2_deg_per_min = 1d/constant('re')*60*constant('deg')
    foreach key, timing_info.keys() do begin
        tinfo = timing_info[key]
        probes = tinfo.probes
        nprobe = n_elements(probes)
        ref_times = tinfo.ref_times
        ref_rsms = fltarr(nprobe,ndim)
        foreach probe, probes, ii do begin
            prefix = probe+'_'
            ref_rsms[ii,*] = get_var_data(prefix+'r_sm', at=ref_times[ii])
        endforeach
        ref_mlts = azim_df_calc_pseudo_mlt(ref_rsms)
        tinfo['ref_rsms'] = ref_rsms
        tinfo['ref_mlts'] = ref_mlts
        tinfo['ref_rxys'] = snorm(ref_rsms[*,0:1])

    ;---Sort probes.
        tinfo['mlt_sorted_probes'] = probes[sort(abs(ref_mlts))]
        tinfo['time_sorted_probes'] = probes[sort(ref_times)]

    ;---Linear fit.
        xxs = ref_times
        yys = ref_mlts
        fit_result = linfit(xxs, yys, yfit=yfit)
        ss_res = total((yys-yfit)^2)
        ss_tot = total((yys-mean(yys))^2)
        r_square = 1-ss_res/ss_tot
        omega_fit = fit_result[1]*hr_per_sec_2_deg_per_min
        linear_fit_info = dictionary($
            'ref_times', ref_times, $
            'ref_mlts', ref_mlts, $
            'omega_fit', omega_fit, $
            'fit_result', fit_result, $
            'r_square', r_square)
        tinfo['linear_fit_info'] = linear_fit_info
        tinfo['omega_fit'] = omega_fit

    ;---3-spacecraft timing.
        nvertex = 3
        ref_rsms[*,2] = 0   ; force x-y plane.
        combos = choose_from(probes, nvertex)
        combo_infos = dictionary()
        foreach combo, combos do begin
            sorted_combo = combo[sort(combo)]
            key = strjoin(sorted_combo,'_')
            combo_index = []
            foreach probe, combo do combo_index = [combo_index, where(probes eq probe)]
            vertex_rsms = ref_rsms[combo_index,*]
            vertex_center_rsms = reform(total(vertex_rsms,1)/nvertex)
            triad_angles = triangle_angles(reform(transpose(vertex_rsms), [1,nvertex,ndim]))
            vertex_times = ref_times[combo_index]
            triad_time_diffs = abs(vertex_times-vertex_times[shift(findgen(nvertex),1)])

            ; Solve for vsm and vmag.
            time_index = sort(vertex_times)
            tt = vertex_times[time_index]
            rr = vertex_rsms[time_index,0:1]
            tt = tt[1:*]-tt[0]
            rr = transpose(rr[1:*,*]-(rr[0,*] ## [1,1]))
            vv = la_linear_equation(rr,tt)
            v_hat = sunitvec(vv)
            rr_normal = dblarr(ndim-1)
            for ii=0, ndim-2 do rr_normal[ii] = sdot(rr[*,ii],v_hat)
            fit_result = linfit(tt, rr_normal)
            v_mag = fit_result[1]*constant('re')
            rxy_center = [vertex_center_rsms[0:1],0]
            omega_triad = snorm(vec_cross([v_hat,0],sunitvec(rxy_center)))*v_mag/snorm(rxy_center)*km_per_s_2_deg_per_min

            combo_info = dictionary()
            combo_info['probes'] = combo
            combo_info['times'] = vertex_times
            combo_info['rsms'] = vertex_rsms
            combo_info['center_rsm'] = vertex_center_rsms
            combo_info['angles'] = triad_angles
            combo_info['time_diffs'] = triad_time_diffs
            combo_info['v_hat'] = v_hat
            combo_info['v_mag'] = v_mag
            combo_info['omega_triad'] = omega_triad
            combo_infos[key] = combo_info
        endforeach

        good_triad_flags = intarr(n_elements(combos))+1
        good_triad_angle_range = [15.,165]
        good_triad_angle_range = project.search_candidate.search_triad.triad_angle_range
        good_triad_tdiff_range = [70.,1e10]
        foreach key, combo_infos.keys(), ii do begin
            combo_info = combo_infos[key]
            index = where_pro(combo_info.angles, '[]', good_triad_angle_range, count=count)
            if count ne nvertex then good_triad_flags[ii] = 0
            index = where_pro(combo_info.time_diffs, '[]', good_triad_tdiff_range, count=count)
            if count ne nvertex then good_triad_flags[ii] = 0
        endforeach
        tinfo['good_triad_flags'] = good_triad_flags
        tinfo['triad_timing'] = combo_infos
    endforeach

    timing_info['good_triad_angle_range'] = good_triad_angle_range
    timing_info['good_triad_tdiff_range'] = good_triad_tdiff_range
    event['timing_info'] = timing_info

    the_key = 'max_ae'
    if ~event.haskey(the_key) then begin
        omni_read_index, event.time_range
        event[the_key] = max(get_var_data('ae'))
    endif

    store_data, event_info_var, 0, event
    tplot_save, event_info_var, filename=data_file
    return, event

end
