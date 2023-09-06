;+
; Detect ramps over the project.
;-

pro azim_df_detect_ramp, project=project, $
    settings=settings, return_setting=return_setting, reset=reset


    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(settings) eq 0 then begin
        overall_roi = project.overall_roi
        settings = dictionary($
            'boxcar_window', 120., $    ; sec.
            'boxcar_ratio', 0.8, $      ; #.
            'dtheta_nsigma', 1., $      ; to select significant peaks of positive slope.
            'df_min_duration', 60. , $
            'pdyn', overall_roi.pdyn, $
            'rxy_range', overall_roi.rxy_range, $
            'mlt_range', overall_roi.mlt_range, $
            'roi_min_duration', project.roi_min_duration )
    endif

    boxcar_window = settings.boxcar_window
    boxcar_ratio = settings.boxcar_ratio
    df_min_duration = settings.df_min_duration
    dtheta_nsigma = settings.dtheta_nsigma

    pdyn = settings.pdyn
    rxy_range = settings.rxy_range
    mlt_range = settings.mlt_range
    roi_min_duration = settings.roi_min_duration

    if keyword_set(return_setting) then return


;---Search according to type. Do not search all through because the data gap in between.
    search_types = project.search_types
    probe_infos = project.probe_infos
    foreach search_type, search_types do begin
        search_name = search_type.name
        probes = search_type.probes
        foreach probe, probes do begin
        ;---Check out file status.
            out_file = join_path([project.data_dir,project.name+'_detect_ramp_'+probe+'_'+search_name+'.txt'])
            log_file = -1

            ; Reset files.
            if keyword_set(reset) then begin
                lprmsg, 'Resetting DF search ...'
                file_delete, out_file, /allow_nonexistent
                lprmsg, 'Clear memory ...'
                del_data, '*'
            endif

            ; Check empty files.
            if file_test(out_file) eq 1 then begin
                if file_lines(out_file) eq 0 then file_delete, out_file
            endif

            ; Already done detection.
            if file_test(out_file) eq 1 then continue
;stop    ; To prevent re-detect ramps again by mistake.

        ;---New detection.
            if file_test(out_file) eq 0 then ftouch, out_file

            lprmsg, 'Settings for detecting ramps ...', log_file
            lprmsg, tab+'Probe: '+strupcase(probe)+' ...', log_file
            lprmsg, '', log_file
            lprmsg, tab+'Window size (sec): '+string(boxcar_window,format='(I0)'), log_file
            lprmsg, tab+'Points used within window: '+string(boxcar_ratio*100,format='(I0)')+'% around mean', log_file
            lprmsg, tab+'Min ramp duration (sec): '+string(df_min_duration,format='(I0)'), log_file
            lprmsg, tab+'Significant peaks in dtheta are above : '+string(dtheta_nsigma,format='(I0)')+' sigma', log_file
            lprmsg, '', log_file
            lprmsg, 'ROI: ', log_file
            lprmsg, tab+'MLT (hr): '+strjoin(string(mlt_range, format='(I0)'),' to ')+' ...', log_file
            lprmsg, tab+'Rxy (Re): '+strjoin(string(rxy_range, format='(I0)'),' to ')+' ...', log_file
            lprmsg, tab+'Pdyn (nPa): '+string(pdyn,format='(I0)')+' ...', log_file
            lprmsg, tab+'Min duration in ROI (hr): '+string(roi_min_duration/3600., format='(I0)')+' ...', log_file
            lprmsg, '', log_file
            lprmsg, 'Writing results to file: '+out_file+' ...', log_file


        ;---Time range to be searched.
            full_time_range = search_type.time_range
            time_step = project.time_step
            common_times = make_bins(full_time_range, time_step)
            ncommon_time = n_elements(common_times)
            lprmsg, '', log_file
            lprmsg, tab+'Time range: '+strjoin(time_string(full_time_range,tformat='YYYY-MM-DD'), ' to ')+' ...', log_file

        ;---Load xxx_theta.
            lprmsg, 'Loading orbit and magnetic field data ...'
            azim_df_read_data, 'theta', time_range=full_time_range, probe=probe, project=project
            prefix = probe_infos[probe].prefix
            interp_time, prefix+'theta', common_times
            azim_df_read_data, 'r_sm', time_range=full_time_range, probe=probe, project=project
            azim_df_read_data, 'r_gsm', time_range=full_time_range, probe=probe, project=project

        ;---Remove data outside the general ROI.
            roi_flags = bytarr(n_elements(common_times))+1

            ; Magnetopause.
            r_gsm_var = prefix+'r_gsm'
            r_gsm = get_var_data(r_gsm_var, at=common_times)
            index = where(check_if_in_magn(r_gsm, dynamic_pressure=pdyn) eq 0, count)
            if count ne 0 then roi_flags[index] = 0
            ; rxy.
            r_sm_var = prefix+'r_sm'
            r_sm = get_var_data(r_sm_var, at=common_times)
            rxy = snorm(r_sm[*,0:1])
            index = where_pro(rxy, '][', rxy_range, count=count)
            if count ne 0 then roi_flags[index] = 0
            ; mlt.
            mlt = azim_df_calc_pseudo_mlt(r_sm)
            index = where_pro(mlt, '][', mlt_range, count=count)
            if count ne 0 then roi_flags[index] = 0

            index = where(roi_flags eq 1, count)
            if count eq 0 then begin
                lprmsg, 'No data in ROI ...'
                continue
            endif

            roi_times = common_times[index]
            roi_time_ranges = time_to_range(roi_times, time_step=time_step)
            roi_durations = roi_time_ranges[*,1]-roi_time_ranges[*,0]
            index = where(roi_durations ge roi_min_duration, count)
            if count eq 0 then begin
                lprmsg, 'No data in ROI for long enough ...'
                continue
            endif
            roi_time_ranges = roi_time_ranges[index,*]
            nroi_time_range = n_elements(roi_time_ranges)*0.5

        ;---Detect dipolarizations.
            boxcar_width = boxcar_window/time_step
            for roi_id=0,nroi_time_range-1 do begin
                time_range = reform(roi_time_ranges[roi_id,*])
                lprmsg, 'Processing ROI section '+strjoin(time_string(time_range),' to ')+' ...'

                ; Load data.
                theta_var = prefix+'theta'
                azim_df_smooth_theta, theta_var, time_range, smooth_window=boxcar_window, stat_ratio=boxcar_ratio
                theta = get_var_data(theta_var, in=time_range, times=times)
                theta_smooth_combo_var = theta_var+'_smooth_combo'
                tmp = get_var_data(theta_smooth_combo_var)
                theta_smooth = tmp[*,0]
                theta_stddev = tmp[*,1]
                dtheta = tmp[*,2]
                dtheta_stddev = tmp[*,3]
                if keyword_set(test_time) then begin
                    theta_combo_var = prefix+'theta_combo'
                    store_data, theta_combo_var, times, [[theta],[theta_smooth],[theta_stddev],[-theta_stddev]], $
                        limits={ytitle:'(deg)',labels:tex2str('theta')+['orig',' smooth',' upper', ' lower'],colors:sgcolor(['silver','red','tan','tan'])}
                endif
                if keyword_set(test_time) then begin
                    dtheta_combo_var = prefix+'dtheta_combo'
                    store_data, dtheta_combo_var, times, [[dtheta],[dtheta_stddev*dtheta_nsigma]], $
                        limits={ytitle:'(deg/sec)',labels:'d'+tex2str('theta')+['',' stddev'+tex2str('times')+string(dtheta_nsigma,format='(I0)')],colors:sgcolor(['black','blue'])}
                endif

                ; Pick out the times when dtheta has significant peaks.
                index = where(dtheta gt dtheta_stddev*dtheta_nsigma, count)
                if count eq 0 then begin
                    lprmsg, tab+tab+'No ramp with significant positive slope ...'
                endif
                time_ranges = time_to_range(times[index], time_step=time_step)
                durations = time_ranges[*,1]-time_ranges[*,0]
                index = where(durations gt df_min_duration, ntime_range)
                if ntime_range eq 0 then begin
                    lprmsg, tab+tab+'No ramp lasts long enough ...'
                endif
                time_ranges = time_ranges[index,*]


                ; Pick out the ramps that have nodes.
                for ii=0, ntime_range-1 do begin
                    the_time_range = reform(time_ranges[ii,*])

                    ; Exclude the time ranges that are on the edges.
                    ;if the_time_range[0] eq time_range[0] then continue
                    ;if the_time_range[1] eq time_range[1] then continue
                    if (the_time_range[0]-time_range[0]) le boxcar_window then continue
                    if (time_range[1]-the_time_range[1]) le boxcar_window then continue

                    ; Find if there is a node.
                    index = where_pro(times, '[]', the_time_range, count=npoint)
                    min_value_index = index[0]
                    max_value_index = index[npoint-1]
                    min_value = theta_smooth[min_value_index]
                    max_value = theta_smooth[max_value_index]
                    min_value_stddev = theta_stddev[min_value_index]
                    max_value_stddev = theta_stddev[max_value_index]
                    if min_value ge  min_value_stddev then continue
                    if max_value le -max_value_stddev then continue
                    if (max_value-min_value) le (min_value_stddev+max_value_stddev) then continue   ; value change is too small.
                    the_value_range = [min_value, max_value]

                    obs_time = mean(the_time_range)
                    width = total(the_time_range*[-1,1])
                    height = total(the_value_range*[-1,1])
                    obs_r_sm = get_var_data(prefix+'r_sm', at=obs_time)

                    the_ramp = dictionary($
                        'probe', probe, $
                        'time_range', the_time_range, $
                        'value_range', the_value_range, $
                        'obs_time', obs_time, $
                        'obs_r_sm', obs_r_sm, $
                        'width', width, $
                        'height', height)

                    if keyword_set(test_time) then begin
                        vars = [theta_combo_var,dtheta_combo_var]
                        nvar = n_elements(vars)
                        poss = sgcalcpos(nvar)

                        theta_combo = get_var_data(theta_combo_var, times=times)
                        theta_smooth = theta_combo[*,1]
                        yrange = minmax(theta_combo)
                        options, theta_combo_var, 'yrange', yrange
                        options, theta_combo_var, 'ystyle', 1

                        tplot, vars, trange=time_range, position=poss

                        tpos = poss[*,0]
                        plot, time_range, yrange, position=tpos, $
                            xstyle=5, $
                            ystyle=5, $
                            /noerase, /nodata
                        the_time_range = the_ramp.time_range
                        index = where_pro(times, '[]', the_time_range)
                        xxs = times[index]
                        yys = theta_smooth[index]
                        plots, xxs,yys, color=sgcolor('green'), thick=2
                        stop
                    endif

                ;---Print output.
                    azim_df_ramp_write, the_ramp, filename=out_file
                endfor; Done loop time_ranges.
            endfor; Done loop roi_time_ranges.
        endforeach; Done loop probes.
    endforeach; Done loop search_types.

end


azim_df_detect_ramp
end
