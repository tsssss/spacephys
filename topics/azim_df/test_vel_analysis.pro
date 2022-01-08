;+
; Need a probes and ref_times, then figure out:
;   1. All combos of triads and v_3sc
;   2. v_fit
;-


;---Prepare inputs.
    inputs = [$
        'thd'  , '2014-08-28/10:11', $
        'the'  , '2014-08-28/10:13', $
        'g15'  , '2014-08-28/10:15', $
        'tha'  , '2014-08-28/10:20', $
        'rbspb', '2014-08-28/10:45', $
        'g13'  , '2014-08-28/10:47']

;    ; This is based on azim_prop_v2019_1112 Fig 2.
;    inputs = [$
;        'thd'  , '2014-08-28/10:10:00', $
;        'g15'  , '2014-08-28/10:13:15', $
;        'the'  , '2014-08-28/10:15:25', $
;        'tha'  , '2014-08-28/10:17:10', $
;        'rbspb', '2014-08-28/10:39:00', $
;        'g13'  , '2014-08-28/10:43:20']


    inputs = [$
        'rbspa', '2016-10-13/12:22', $
        'the'  , '2016-10-13/12:26', $
        'g15'  , '2016-10-13/12:29', $
        'thd'  , '2016-10-13/12:31', $
        'g14'  , '2016-10-13/12:37', $
        'g13'  , '2016-10-13/12:45']
        
    ; Event 1. 
    inputs = [$
        'the'  , '2008-02-29/08:26', $
        'thd'  , '2008-02-29/08:27', $
        'thc'  , '2008-02-29/08:29', $
        'tha'  , '2008-02-29/08:31']
        
    ; Event 2.
    inputs = [$
        'tha'  , '2008-03-23/02:02', $
        'thc'  , '2008-03-23/02:09', $
        'thd'  , '2008-03-23/02:13', $
        'the'  , '2008-03-23/02:14']
        
    ; Event 3.
    inputs = [$
        'thc'  , '2008-03-27/00:22', $
        'tha'  , '2008-03-27/00:24', $
        'the'  , '2008-03-27/00:35', $
        'thd'  , '2008-03-27/00:36']  
    
    ; Event 4.
    inputs = [$
        'rbspb', '2013-10-09/09:01', $
        'rbspa', '2013-10-09/09:04', $
        'tha'  , '2013-10-09/08:48', $
        'g15'  , '2013-10-09/08:39', $
        'thd'  , '2013-10-09/08:40']
    
    ; Event 5.
    inputs = [$
        'the'  , '2013-10-15/06:10', $
        'g15'  , '2013-10-15/06:13', $
        'tha'  , '2013-10-15/06:13', $
        'thd'  , '2013-10-15/06:22', $
        'rbspb', '2013-10-15/06:39']
    
    ; Event 7.
    inputs = [$
        'tha'  , '2016-06-15/04:39', $
        'g15'  , '2016-06-15/04:37', $
        'g14'  , '2016-06-15/04:33', $
        'g13'  , '2016-06-15/04:19', $
        'mms1' , '2016-06-15/04:20']
        
    ; Event 9.
    inputs = [$
        'tha'  , '2016-12-11/09:46', $
        'thd'  , '2016-12-11/09:50', $
        'the'  , '2016-12-11/09:52', $
        'g13'  , '2016-12-11/10:04']

    ; Event 11.
    inputs = [$
        'rbspa', '2017-04-01/15:08', $
        'rbspb', '2017-04-01/15:09', $
        'thd'  , '2017-04-01/15:04', $
        'tha'  , '2017-04-01/15:02', $
        'the'  , '2017-04-01/15:08']
    
    ; Event 12.
    inputs = [$
        'rbspa', '2017-04-01/15:50', $
        'rbspb', '2017-04-01/15:51', $
        'thd'  , '2017-04-01/15:40', $
        'tha'  , '2017-04-01/15:44', $
        'the'  , '2017-04-01/15:45']

    ; Event 13.
    inputs = [$
        'the'  , '2017-04-09/12:18', $
        'tha'  , '2017-04-09/12:24', $
        'thd'  , '2017-04-09/12:25', $
        'rbspb', '2017-04-09/12:37']    

    nprobe = n_elements(inputs)/2
    inputs = reform(inputs, 2,nprobe)
    probes = reform(inputs[0,*])
    ref_times = time_double(reform(inputs[1,*]))
    if n_elements(project) eq 0 then project = azim_df_load_project()
    time_range = minmax(ref_times)
    plot_time_range = minmax(make_bins(time_range,600))
    data_time_range = minmax(make_bins(time_range,constant('secofday')))

    event_info = dictionary()
    event_info['probes'] = probes
    event_info['ref_times'] = ref_times
    event_info['plot_time_range'] = plot_time_range
    event_info['data_time_range'] = data_time_range


;---Read position data.
    cdf_file = (ref_times[0] ge time_double('2010-01-01'))? $
        'azim_df_primitive_data.cdf': 'azim_df_primitive_data_themis.cdf'
    cdf_file = join_path([project.data_dir,cdf_file])
    foreach probe, probes do begin
        azim_df_read_data, 'r_sm', probe=probe, time_range=data_time_range, project=project
    endforeach

;---Get positions.
    ndim = 3
    ref_rsms = fltarr(nprobe,ndim)
    ref_mlts = fltarr(nprobe)
    foreach probe, probes, ii do begin
        ref_rsms[ii,*] = get_var_data(probe+'_r_sm', at=ref_times[ii])
        ref_mlts[ii] = azim_df_calc_pseudo_mlt(ref_rsms[ii,*])
    endforeach
    event_info['ref_mlts'] = ref_mlts
    event_info['ref_rsms'] = ref_rsms

;---Linear fit.
    hr_sec_2_deg_min = 15d*60   ; hr/sec to deg/min.
    xxs = ref_times
    yys = ref_mlts
    fit_result = linfit(xxs,yys, yfit=yfit)
    omega_fit = fit_result[1]*hr_sec_2_deg_min
    ss_res = total((yys-yfit)^2)
    ss_tot = total((yys-mean(yys))^2)
    fit_r2 = 1-ss_res/ss_tot
    event_info['omega_fit'] = omega_fit
    event_info['fit_r2'] = fit_r2

;---Combos.
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
        triad_angles = triangle_angles(reform(vertex_rsms, [1,nvertex,ndim]))
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

        combo_info = dictionary()
        combo_info['probes'] = combo
        combo_info['times'] = vertex_times
        combo_info['rsms'] = vertex_rsms
        combo_info['center_rsm'] = vertex_center_rsms
        combo_info['angles'] = triad_angles
        combo_info['time_diffs'] = triad_time_diffs
        combo_info['v_hat'] = v_hat
        combo_info['v_mag'] = v_mag
        combo_infos[key] = combo_info
    endforeach

    xrange = [-1,1]*5
    yrange = [-1,1]*5
    foreach combo_info, combo_infos do begin
        xrange = [xrange, combo_info.rsms[*,0]]
        yrange = [yrange, combo_info.rsms[*,1]]
    endforeach
    xstep = 5
    ystep = 5
    xrange = reverse(minmax(make_bins(xrange, xstep)))
    yrange = reverse(minmax(make_bins(yrange, ystep)))
    plot, xrange, yrange, /nodata, /iso, xrange=xrange, yrange=yrange

    foreach combo_info, combo_infos do begin
        rsms = combo_info.rsms
        center_rsm = combo_info.center_rsm
        v_mag = combo_info.v_mag
        v_hat = combo_info.v_hat
        plots, rsms[*,0], rsms[*,1], psym=1, color=sgcolor('red')
        x0 = center_rsm[0]
        y0 = center_rsm[1]
        scale = 0.075
        x1 = x0+v_hat[0]*scale*v_mag
        y1 = y0+v_hat[1]*scale*v_mag
        arrow, x0,y0, x1,y1, /data, hsize=5, /solid
    endforeach
    stop

end
