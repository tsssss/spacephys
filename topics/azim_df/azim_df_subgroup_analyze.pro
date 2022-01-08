;+
;-

function azim_df_subgroup_analyze, subgroup, project=project, log_file=log_file, test_time=test_time, settings=settings

    if n_elements(log_file) eq 0 then log_file = -1

    retval = dictionary()
    lprmsg, '', log_file
    lprmsg, 'Processing candidate '+string(subgroup.id), log_file
    lprmsg, strjoin(time_string(subgroup.time_range),' to ')+' ...', log_file

;test_time = time_double('2014-08-28/10:30')
;test_time = time_double('2016-10-13/12:30')
if keyword_set(test_time) then begin
    if product(subgroup.time_range-test_time) ge 0 then return, retval else stop
endif

    tab = constant('4space')
    common_data_rate = project.time_step

    if n_elements(settings) eq 0 then settings = dictionary()
    obs_time_diff_range = settings.haskey('obs_time_diff_range')? settings.obs_time_diff_range: [-1.,1800]
    obs_time_max_diff = settings.haskey('obs_time_max_diff')? settings.obs_time_max_diff: 900
    obs_time_max_diff_count = settings.haskey('obs_time_max_diff_count')? settings.obs_time_max_diff_count: 2
    xcor_section_ratio = settings.haskey('xcor_section_ratio')? settings.xcor_section_ratio: [-1.5,1.5]
    xcor_min_value = settings.haskey('xcor_min_value')? settings.xcor_min_value: 0.5
    max_xcor_error = settings.haskey('max_xcor_error')? settings.max_xcor_error: 1d/3
    max_xcor_error_abs = settings.haskey('max_xcor_error_abs')? settings.max_xcor_error_abs: 60
    triad_angle_range = settings.haskey('triad_angle_range')? settings.triad_angle_range: [0.,180]+[1,-1]*15


;---Check if need to load data.
    search_name = subgroup.search_name
    the_var = 'current_search_name'
    if get_var_data(the_var) ne search_name then begin
        if search_name eq 'beyond_15Re' then begin
            full_time_range = time_double(['2007-03-01','2009-09-01'])
            all_probes = 'th'+letters('e')
            data_file_suffix = 'azim_df_primitive_data_themis.cdf'
        endif else begin
            full_time_range = time_double(['2012-10-01','2017-10-01'])
            all_probes = ['rbsp'+letters('b'),'th'+['a','d','e'],'mms1','g'+['13','14','15']]   ; turns out no thb/thc candidate.
            data_file_suffix = 'azim_df_primitive_data.cdf'
        endelse

        time_step = project.time_step
        common_times = make_bins(full_time_range, time_step)
        probe_infos = project.probe_infos
        foreach probe, all_probes do begin
            lprmsg, 'Processing '+strupcase(probe)+' ...'
            prefix = probe_infos[probe].prefix
            if check_if_update(prefix+'theta', full_time_range) then begin
                azim_df_read_data, 'theta', time_range=full_time_range, probe=probe, project=project
                interp_time, prefix+'theta', common_times
            endif
            if check_if_update(prefix+'r_sm', full_time_range) then begin
                azim_df_read_data, 'r_sm', time_range=full_time_range, probe=probe, project=project
                interp_time, prefix+'r_sm', common_times
            endif
            if check_if_update(prefix+'pseudo_mlt', full_time_range) then begin
                r_sm = get_var_data(prefix+'r_sm', times=times)
                mlt = azim_df_calc_pseudo_mlt(r_sm)
                store_data, prefix+'pseudo_mlt', times, mlt
            endif
        endforeach
        store_data, the_var, 0, search_name
    endif


;---Gather info of the current subgroup.
    probes = subgroup.probes
    df_list = subgroup.df_list
    nprobe = n_elements(probes)

    ndim = 3
    df_r_sms = fltarr(nprobe,ndim)
    foreach df, df_list, ii do df_r_sms[ii,*] = df.obs_r_sm
    df_obs_times = dblarr(nprobe)
    foreach df, df_list, ii do df_obs_times[ii] = df.obs_time
    df_widths = fltarr(nprobe)
    foreach df, df_list, ii do df_widths[ii] = df.width


;---Filter by dt.
    dtimes = df_obs_times[1:nprobe-1]-df_obs_times[0:nprobe-2]
    msg = tab+'obs_time diff (sec): '
    foreach tmp, dtimes do msg += string(tmp,format='(I0)')+' '
    lprmsg, msg, log_file
    index = lazy_where(dtimes, ')(', obs_time_diff_range, count=count)
    if count ne 0 then begin
        lprmsg, tab+'obs_time diff too large or too small, skip ...', log_file
        return, retval
    endif else lprmsg, tab+'obs_time diff good ...', log_file
    index = where(dtimes ge obs_time_max_diff, count)
    if count ge obs_time_max_diff_count then begin
        lprmsg, tab+'obs_time diff too large for '+string(obs_time_max_diff_count,format='(I0)')+' times, skip ...', log_file
        return, retval
    endif else lprmsg, tab+'obs_time diff good ...', log_file
    index = where(dtimes[[0,nprobe-2]] ge obs_time_max_diff, count)
    if count ne 0 then begin
        lprmsg, tab+'obs_time large at beginning or end ...', log_file
        return, retval
    endif


;---Do xcor for all edges.
    nedge_vertex = 2
    probe_list = list(probes,/extract)
    edge_combos = choose_from(probe_list.toarray(), nedge_vertex)
    edge_list = dictionary()
    foreach edge_combo, edge_combos do begin
        ; Prepare the key.
        edge_probes = edge_combo[sort(edge_combo)]
        key = strjoin(edge_probes,'_')

        ; Get info about the vertex.
        probe_index = []
        foreach probe, edge_probes do probe_index = [probe_index, where(probes eq probe)]
        vertex_obs_times = df_obs_times[probe_index]
        vertex_widths = df_widths[probe_index,*]

        ; Sort by obs_time.
        index = sort(vertex_obs_times)
        edge_probes = edge_probes[index]
        vertex_obs_times = vertex_obs_times[index]
        vertex_widths = vertex_widths[index]
        ; Must be positive or 0.
        time_lag_df = total(vertex_obs_times*[-1,1])

    ;---xcor.
        the_width = min(vertex_widths)  ; the mean width.
        section_time_range = xcor_section_ratio*the_width   ; expand the width.
        section_times = make_bins(section_time_range, common_data_rate)
        nsection_time = n_elements(section_times)

        ; Read data.
        nrec = n_elements(section_times)
        data = fltarr(nrec,nedge_vertex)
        for probe_id=0,nedge_vertex-1 do begin
            the_times = vertex_obs_times[probe_id]+section_times
            the_var = edge_probes[probe_id]+'_theta'
            data[*,probe_id]= get_var_data(the_var, at=the_times)
        endfor

        ; Do xcor.
        nlag = floor(the_width/common_data_rate)
        lags = findgen(nlag)-round(nlag/2)
        xcors = c_correlate(data[*,0], data[*,1], lags)
        xcor_max = max(xcors, index)
        xcor_time_lag = lags[index]*common_data_rate
        ; xcor error.
        xx = total(data,2)/nedge_vertex
        del_x = xx-mean(xx)
        dx_dt = deriv(xx)/common_data_rate
        xcor_err = sqrt(1./(nrec-1)*(1-xcor_max)/xcor_max*2*mean(del_x^2)/mean(dx_dt^2))
        if index eq 0 or index eq nlag-1 then begin
            fillval = !values.f_nan
            xcor_err = fillval
        endif

        time_lag_xcor = time_lag_df+xcor_time_lag
        time_lag_ratio = abs(xcor_time_lag)/max(abs([time_lag_df,time_lag_xcor]))
        if time_lag_df eq 0 and time_lag_xcor eq 0 then time_lag_ratio = 0

        edge_info = dictionary($
            'probes', edge_probes, $
            'time_lag_df', time_lag_df, $   ; = obs_time_p2-obs_time_p1.
            'time_lag_xcor', time_lag_xcor, $   ; = p2-p1.
            'xcor_max', xcor_max, $
            'xcor_err', xcor_err, $
            'time_lag_diff', xcor_time_lag, $
            'time_lag_ratio', time_lag_ratio)
        edge_list[key] = edge_info
    endforeach

    lprmsg, 'All edges: ', log_file
    foreach edge, edge_list do azim_df_edge_write, edge, filename=log_file, prefix=tab
    nedge = edge_list.length


;---Filter by timing of the edges.
    foreach key, edge_list.keys() do begin
        edge = edge_list[key]
        lprmsg, 'Checking edge: '+key+' ...', log_file
;        if ~finite(edge.xcor_err) then begin
;            lprmsg, 'Probes of the current edge do not correlate, skip ...', log_file
;            edge_list.remove, key
;            continue
;        endif
        if edge.xcor_max lt xcor_min_value then begin
            lprmsg, 'xcor_max of the current edge is too small, skip ...', log_file
            edge_list.remove, key
            continue
        endif
        if abs(edge.time_lag_diff) le max_xcor_error_abs or edge.time_lag_ratio le max_xcor_error then begin
            lprmsg, 'time lag is good ...', log_file
        endif else begin
            lprmsg, 'time lag too large, skip ...', log_file
            edge_list.remove, key
        endelse
    endforeach
    lprmsg, tab+'Found '+string(edge_list.length,format='(I0)')+' good edges ...', log_file
    if edge_list.length ne nedge then begin
        lprmsg, 'Not all edges are good, skip ...', log_file
        return, retval
    endif else begin
        lprmsg, 'All edges are good ...', log_file
    endelse


;---All triangles.
    min_triad_count = 2
    ntriad_vertex = 3
    triad_combos = choose_from(probe_list.toarray(), ntriad_vertex)
    triad_list = dictionary()
    foreach triad_probes, triad_combos do begin
        sorted_probes = triad_probes[sort(triad_probes)]
        key = strjoin(sorted_probes,'_')

        ; Get info about the vertex.
        probe_index = []
        foreach probe, triad_probes do probe_index = [probe_index, where(probes eq probe)]
        vertex_obs_times = df_obs_times[probe_index]
        vertex_r_sms = df_r_sms[probe_index,*]

        ; Sort by obs_time.
        index = sort(vertex_obs_times)
        triad_probes = triad_probes[index]
        vertex_obs_times = vertex_obs_times[index]
        vertex_r_sms = vertex_r_sms[index,*]

        triad_angles = triangle_angles(reform(transpose(vertex_r_sms), [1,ntriad_vertex,ndim]))
        triad_angles = reform(triad_angles)
        triad_info = dictionary($
            'probes', triad_probes, $
            'times', vertex_obs_times, $
            'r_sms', vertex_r_sms, $
            'angles', triad_angles)
        triad_list[key] = triad_info
    endforeach


;---Filter by geometry.
    foreach key, triad_list.keys() do begin
        triad_info = triad_list[key]
        index = lazy_where(triad_info.angles, '()', triad_angle_range, count=count)
        if count eq ntriad_vertex then continue
        triad_list.remove, key
    endforeach
    lprmsg, tab+'Found '+string(triad_list.count(),format='(I0)')+' triads with good geometry ...', log_file
    if triad_list.count() lt min_triad_count then begin
        lprmsg, tab+'No enough triads with good geometry, skip ...', log_file
        return, retval
    endif


;---Update probes and df_list.
    probes = list()
    foreach triad_info, triad_list do probes.add, triad_info.probes
    probes = sort_uniq(probes.toarray())
    flags = bytarr(df_list.length)
    foreach df, df_list, df_id do begin
        index = where(probes eq df.probe, count)
        if count ne 0 then continue
        flags[df_id] = 1
    endforeach
    index = where(flags eq 0)
    df_list = df_list[index]

    probes = strarr(df_list.length)
    foreach df, df_list, df_id do probes[df_id] = df.probe


;---Do 3-sc timing.
    re = constant('re')
    km_per_s_2_deg_per_min = 1d/re*60*constant('deg')
    foreach triad_info, triad_list do begin
        the_probes = triad_info.probes  ; already sorted by obs_time.
        obs_times = triad_info.times
        obs_r_sms = triad_info.r_sms[*,0:1]
        center_r_sm = reform(total(obs_r_sms,1)/ntriad_vertex)
        center_rxy = [center_r_sm[0:1],0]
        triad_info['center_r_sm'] = center_r_sm

        ; Solve for v_2d using obs_times.
        rr = transpose(obs_r_sms[1:*,*]-(obs_r_sms[0,*] ## [1,1]))
        tt = obs_times[1:*]-obs_times[0]
        vv = la_linear_equation(rr,tt)
        vhat = sunitvec(vv)
        rr_normal = dblarr(ndim-1)
        for ii=0, ndim-2 do rr_normal[ii] = sdot(rr[*,ii],vhat)
        fit_result = linfit(tt, rr_normal)
        vmag = fit_result[1]*re
        triad_info['vhat_obs_time'] = vhat
        triad_info['vmag_obs_time'] = vmag
        omega = vec_cross([vhat,0],sunitvec(center_rxy))*vmag/snorm(center_rxy)*km_per_s_2_deg_per_min
        triad_info['omega_obs_time'] = omega[2]     ; negative is eastward.

        tt = dblarr(ntriad_vertex-1)
        for ii=0,ntriad_vertex-2 do begin
            edge_probes = [the_probes[0],the_probes[ii+1]]
            key = strjoin(edge_probes[sort(edge_probes)],'_')
            the_edge = edge_list[key]
            tt[ii] = the_edge.time_lag_xcor
            if edge_probes[0] eq the_edge.probes[1] then tt[ii] = -tt[ii]
        endfor
        vv = la_linear_equation(rr,tt)
        vhat = sunitvec(vv)
        rr_normal = dblarr(ndim-1)
        for ii=0, ndim-2 do rr_normal[ii] = sdot(rr[*,ii],vhat)
        fit_result = linfit(tt, rr_normal)
        vmag = fit_result[1]*re
        triad_info['vhat_xcor'] = vhat
        triad_info['vmag_xcor'] = vmag
        omega = vec_cross([vhat,0],sunitvec(center_rxy))*vmag/snorm(center_rxy)*km_per_s_2_deg_per_min
        triad_info['omega_xcor'] = omega[2]
    endforeach


    if keyword_set(test_time) then stop


    lprmsg, 'Passed all criteria ...', log_file
    subgroup['probes'] = probes
    subgroup['df_list'] = df_list
    subgroup['edge_list'] = edge_list
    subgroup['triad_list'] = triad_list
    return, subgroup


end
