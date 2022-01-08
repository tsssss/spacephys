;+
; Study the timing of a subgroup.
;-

function azim_df_analyze_subgroup, subgroup, project=project

    tab = constant('4space')
    retval = !null
    max_probe_length = 5

    xcor_section_ratio = [-1.5,1.5]
    xcor_min_value = 0.5
    common_data_rate = project.time_step


;---Check if need to load data.
    search_name = subgroup.search_type
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
            azim_df_read_data, 'theta', time_range=full_time_range, probe=probe, project=project
            prefix = probe_infos[probe].prefix
            interp_time, prefix+'theta', common_times
            ;azim_df_read_data, 'r_sm', time_range=full_time_range, probe=probe, project=project
        endforeach

        store_data, the_var, 0, search_name
    endif


;---Gather info of the current subgroup.
    probes = subgroup.probes
    df_list = subgroup.df_list
    nprobe = n_elements(probes)

    ndim = 3
    df_rsms = fltarr(nprobe,ndim)
    foreach df, df_list, ii do df_rsms[ii,*] = df.arrival_r_sm
    df_observed_times = dblarr(nprobe)
    foreach df, df_list, ii do df_observed_times[ii] = df.arrival_time
    df_widths = fltarr(nprobe)
    foreach df, df_list, ii do df_widths[ii] = df.width


;---All triangles.
    ntriad_vertex = 3
    triad_combos = choose_from(probes, ntriad_vertex)
    triad_infos = dictionary()
    foreach triad_combo, triad_combos do begin
        sorted_combo = triad_combo[sort(triad_combo)]
        key = strjoin(sorted_combo,'_')
        combo_index = []
        foreach probe, triad_combo do combo_index = [combo_index, where(probes eq probe)]
        vertex_rsms = df_rsms[combo_index,*]
        vertex_center_rsms = reform(total(vertex_rsms,1)/ntriad_vertex)
        triad_angles = triangle_angles(reform(transpose(vertex_rsms), [1,ntriad_vertex,ndim]))
        triad_info = dictionary($
            'probes', triad_combo, $
            'times', vertex_times, $
            'rsms', vertex_rsms, $
            'center_rsm', vertex_center_rsms, $
            'angles', triad_angles)
        triad_infos[key] = triad_info
    endforeach


;---Filter by dt.
    min_dt = 1800.  ;sec.
    df_times = sort_uniq(df_observed_times)
    dtimes = df_times[1:*]-df_times[0:-2]
    index = where(dtimes ge min_dt, count)
    if count ne 0 then begin
        lprmsg, 'dT too large ...'
        return, retval
    endif
    

;---Filter by geometry.
    triad_angle_limit = 15.
    triad_angle_range = [0.,180]+[1,-1]*triad_angle_limit
    min_combo_count = 2
    geom_triads = dictionary()
    foreach key, triad_infos.keys() do begin
        triad_info = triad_infos[key]
        index = lazy_where(triad_info.angles, '()', triad_angle_range, count=count)
        if count lt ntriad_vertex then continue
        geom_triads[key] = triad_info
    endforeach
    if n_elements(geom_triads) lt min_combo_count then begin
        lprmsg, 'No enough triads with good geometry ...'
        return, retval
    endif


;---Collect all envolved probes.
    triad_probes = list()
    foreach triad_info, geom_triads do triad_probes.add, triad_info.probes
    triad_probes = sort_uniq(triad_probes.toarray())


;---De xcor for all edges.
    nedge_vertex = 2
    edge_combos = choose_from(triad_probes, nedge_vertex)
    edge_infos = dictionary()
    foreach edge_combo, edge_combos do begin
        sorted_combo = edge_combo[sort(edge_combo)]
        key = strjoin(sorted_combo,'_')
        combo_index = []
        foreach probe, edge_combo do combo_index = [combo_index, where(probes eq probe)]
        vertex_observed_times = df_observed_times[combo_index]
        vertex_widths = df_widths[combo_index,*]
        time_lag_df = abs(total(vertex_observed_times*[-1,1]))

    ;---xcor.
        the_width = max(vertex_widths)  ; the longer width.
        section_time_range = xcor_section_ratio*the_width   ; expand the width.
        section_times = make_bins(section_time_range, common_data_rate)

        ; Read data.
        nrec = n_elements(section_times)
        data = fltarr(nrec,nedge_vertex)
        for probe_id=0,nedge_vertex-1 do begin
            the_times = vertex_observed_times[probe_id]+section_times
            the_var = edge_combo[probe_id]+'_theta'
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
            xcor_max = fillval
            xcor_time_lag = fillval
            xcor_err = fillval
        endif

        edge_info = dictionary($
            'probes', edge_combo, $
            'time_lag_df', time_lag_df, $
            'time_lag_xcor', time_lag_df+xcor_time_lag, $
            'xcor_max', xcor_max, $
            'xcor_err', xcor_err)
        edge_infos[key] = edge_info
    endforeach

    lprmsg, 'All edges: '
    foreach edge_info, edge_infos do begin
        msg = strjoin(extend_string(edge_info.probes,length=5),' ')
        print, tab+msg, edge_info.time_lag_df, edge_info.time_lag_xcor, edge_info.xcor_err, edge_info.xcor_max
    endforeach

;---Filter by timing.
    min_xcor = 0.5
    xcor_edges = dictionary()
    foreach key, edge_infos.keys() do begin
        edge_info = edge_infos[key]
        if ~finite(edge_info.xcor_max) then continue
        if edge_info.xcor_max lt min_xcor then continue
        xcor_edges[key] = edge_info
    endforeach

    lprmsg, 'All correlated edges: '
    foreach edge_info, xcor_edges do begin
        msg = strjoin(extend_string(edge_info.probes,length=5),' ')
        print, tab+msg, edge_info.time_lag_df, edge_info.time_lag_xcor, edge_info.xcor_err, edge_info.xcor_max
    endforeach

    if edge_infos.count() ne xcor_edges.count() then begin
        lprmsg, 'Not all edges are correlated ...'
        return, retval
    endif


;---Check time lag diff.
    min_dt_ratio = 0.5
    have_bad_edge = 0
    foreach edge_info, xcor_edges do begin
        time_lag_df = edge_info.time_lag_df
        time_lag_xcor = edge_info.time_lag_xcor
        time_lags = [time_lag_df,time_lag_xcor]
        if abs(total(time_lags*[-1,1])/max(abs(time_lags))) ge min_dt_ratio then begin
            have_bad_edge = 1
            break
        endif
    endforeach
    if have_bad_edge then begin
        lprmsg, 'Not all edges have consistent timing ...'
        return, retval
    endif

    lprmsg, 'The given subgroup is good for 3-spacecraft timing ...'
    subgroup['edge_infos'] = edge_infos
    subgroup['edges'] = xcor_edges
    subgroup['triad_infos'] = triad_infos
    subgroup['triads'] = geom_triads
    return, subgroup

end


; Post-midn.
test_time = time_double('2014-08-28/10:30') ; good eg.
test_time = time_double('2016-10-13/12:30') ; good eg.
;test_time = time_double('2016-10-24/12:30') ; bad eg.
;test_time = time_double('2016-11-10/10:20') ; bad eg.
test_time = time_double('2016-11-13/14:30') ; bad eg. ok.
test_time = time_double('2016-12-11/10:00') ; good eg. ok.
;candidates = azim_df_search_post_midn_events()


; Pre-midn.
test_time = time_double('2014-12-26/01:10') ; good eg.
candidates = azim_df_search_pre_midn_events()



search_name = 'within_15Re'
test_candidate_id = list()
secofday = constant('secofday')
foreach candidate, candidates, candidate_id do begin
    if candidate.search_type ne search_name then continue
    the_time_range = candidate.time_range+[1,-1]*secofday
    if product(the_time_range-test_time) lt 0 then begin
        print, time_string(the_time_range)
        test_candidate_id.add, candidate_id
    endif
endforeach
test_candidate_id = test_candidate_id.toarray()
the_id = !null
if n_elements(test_candidate_id) eq 0 then message, 'No test candidate found ...'
if n_elements(test_candidate_id) gt 1 then message, 'More than 1 candidates ...' else the_id = test_candidate_id[0]
subgroup = candidates[the_id]

project = azim_df_load_project()
event = azim_df_analyze_subgroup(subgroup, project=project)

end
