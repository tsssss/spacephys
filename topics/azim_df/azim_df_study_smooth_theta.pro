;+
; Study the different methods to smooth theta: median/mean, for given probe and search_name.
;
; To focus is to select the method that provides the most smooth version.
;-

pro azim_df_study_smooth_theta, project=project, probe=probe, search_name=search_name, $
    test_time=test_time

    if n_elements(project) eq 0 then project = azim_df_load_project()

;---Load data for probe.
    if n_elements(probe) eq 0 then probe = 'rbspb'

    if n_elements(search_name) eq 0 then begin
        case probe of
            'rbspb': search_name = 'within_15Re'
            'tha': search_name = 'beyond_15Re'
            else: search_name = 'within_15Re'
        endcase
    endif

test = 0

    probe_infos = project.probe_infos
    probe_info = probe_infos[probe]
    prefix = probe_info.prefix
    time_step = project.time_step
    boxcar_window = 120.    ; sec.
    boxcar_width = boxcar_window/time_step
    boxcar_ratio = 0.8      ; use the value below this ratio.

    search_settings = project.search_settings
    foreach search_setting, search_settings do if search_setting.name eq search_name then break


    full_time_range = search_setting.time_range
    common_times = make_bins(full_time_range, time_step)
    theta_var = prefix+'theta'
    if check_if_update(theta_var, full_time_range) then begin
        azim_df_read_data, 'theta', time_range=full_time_range, probe=probe
    endif


;---Filter candidates.
    candidates = azim_df_search_candidate(project=project)

    ; Use test_time to check particular event.
;test_time = time_double('2016-10-13/12:30')
;test_time = time_double('2017-05-02/02:05')
    if keyword_set(test_time) then begin
        candidate_ids = list()
        foreach candidate, candidates do begin
            if product(candidate.time_range-test_time) gt 0 then continue else candidate_ids.add, candidate.id
        endforeach
        candidate_ids = candidate_ids.toarray()-1   ; id starts from 1.
        candidates = candidates[candidate_ids]
    endif else begin
        candidate_ids = list()
        foreach candidate, candidates do begin
            if candidate.search_type ne search_name then continue
            probes = candidate.all_probes
            index = where(probes eq probe, count)
            if count ne 0 then candidate_ids.add, candidate.id
        endforeach
        candidate_ids = candidate_ids.toarray()-1
        candidates = candidates[candidate_ids]
    endelse


;---Loop through.
    log_file = join_path([project.data_dir,project.name+'_study_mean_median_on_'+search_setting.name+'_'+probe+'.txt'])
    if keyword_set(test_time) or file_test(log_file) eq 0 then begin
        if keyword_set(test_time) then log_file = -1 else if file_test(log_file) eq 0 then ftouch, log_file
        tab = constant('4space')

        foreach candidate, candidates do begin
            time_range = candidate.time_range
            lprmsg, 'Processing candidate in '+strjoin(time_string(time_range,tformat='YYYY-MM-DD/hh:mm'),' to ')+' ...'

        ;---Fixed boxcar median.
            theta_median_var = prefix+'theta_median'
            if check_if_update(theta_median_var, time_range) then begin
                theta = get_var_data(prefix+'theta')
                boxcar_boundary_times = make_bins(time_range, boxcar_window)
                nboxcar = n_elements(boxcar_boundary_times)-1
                boxcar_center_times = boxcar_boundary_times[0:nboxcar-1]+boxcar_window*0.5
                theta_median = fltarr(nboxcar)+!values.f_nan
                for ii=0, nboxcar-1 do begin
                    index = lazy_where(common_times, '[]', boxcar_boundary_times[ii:ii+1], count=count)
                    if count eq 0 then continue
                    the_theta = theta[index]
                    index = where(finite(the_theta),count)
                    if count lt 2 then continue
                    theta_median[ii] = median(the_theta)
                endfor
                store_data, theta_median_var, boxcar_center_times, theta_median
                add_setting, theta_median_var, /smart, {$
                    display_type: 'scalar', $
                    short_name: tex2str('theta')+' median', $
                    boxcar_window: boxcar_window, $
                    unit: 'deg'}
            endif

        ;---Fixed boxcar mean.
            theta_mean_var = prefix+'theta_mean'
            if check_if_update(theta_mean_var, time_range) then begin
                theta = get_var_data(prefix+'theta')
                boxcar_boundary_times = make_bins(time_range, boxcar_window)
                nboxcar = n_elements(boxcar_boundary_times)-1
                boxcar_center_times = boxcar_boundary_times[0:nboxcar-1]+boxcar_window*0.5
                theta_mean = fltarr(nboxcar)+!values.f_nan
                for ii=0, nboxcar-1 do begin
                    index = lazy_where(common_times, '[]', boxcar_boundary_times[ii:ii+1], count=count)
                    if count eq 0 then continue
                    the_theta = theta[index]
                    index = where(finite(the_theta),count)
                    if count lt 2 then continue
                    theta_mean[ii] = mean(the_theta)
                endfor
                store_data, theta_mean_var, boxcar_center_times, theta_mean
                add_setting, theta_mean_var, /smart, {$
                    display_type: 'scalar', $
                    short_name: tex2str('theta')+' mean', $
                    boxcar_window: boxcar_window, $
                    unit: 'deg'}
            endif

        ;---Fixed boxcar mean with ratio.
            theta_mean_ratio_var = prefix+'theta_mean_ratio'
            if check_if_update(theta_mean_ratio_var, time_range) then begin
                theta = get_var_data(prefix+'theta')
                boxcar_boundary_times = make_bins(time_range, boxcar_window)
                nboxcar = n_elements(boxcar_boundary_times)-1
                boxcar_center_times = boxcar_boundary_times[0:nboxcar-1]+boxcar_window*0.5
                theta_mean_ratio = fltarr(nboxcar)+!values.f_nan
                for ii=0, nboxcar-1 do begin
                    index = lazy_where(common_times, '[]', boxcar_boundary_times[ii:ii+1], count=count)
                    if count eq 0 then continue
                    the_theta = theta[index]
                    index = where(finite(the_theta),count)
                    if count lt 2 then continue
                    the_theta = the_theta[index]
                    index = sort(abs(the_theta))
                    the_theta = the_theta[index]
                    theta_mean_ratio[ii] = mean(the_theta[0:count*boxcar_ratio])
                endfor
                store_data, theta_mean_ratio_var, boxcar_center_times, theta_mean_ratio
                add_setting, theta_mean_ratio_var, /smart, {$
                    display_type: 'scalar', $
                    short_name: tex2str('theta')+' mean ratio', $
                    boxcar_window: boxcar_window, $
                    unit: 'deg'}
            endif


        ;---Mean slide.
            theta_mean_slide_var = prefix+'theta_slide'
            if check_if_update(theta_mean_slide_var,time_range) then begin
                slide_time_range = time_range+[-1,1]*boxcar_window
                theta = get_var_data(prefix+'theta', in=slide_time_range, times=times)
                boxcar_width = boxcar_window/time_step
                theta_slide = sliding_boxcar(theta, boxcar_width, type='mean')
                index = lazy_where(times, '[]', time_range)
                theta_slide = theta_slide[index]
                times = times[index]
                store_data, theta_mean_slide_var, times, theta_slide
                add_setting, theta_mean_slide_var, /smart, {$
                    display_type: 'scalar', $
                    short_name: tex2str('theta')+' mean slide', $
                    boxcar_window: boxcar_window, $
                    unit: 'deg'}
            endif


        ;---Mean ratio slide.
            theta_mean_ratio_slide_var = prefix+'theta_mean_ratio_slide'
            if check_if_update(theta_mean_ratio_slide_var,time_range) then begin
                slide_time_range = time_range+[-1,1]*boxcar_window
                theta = get_var_data(prefix+'theta', in=slide_time_range, times=times)
                boxcar_width = boxcar_window/time_step
                theta_mean_ratio_slide = sliding_boxcar(theta, boxcar_width, type='mean', ratio=boxcar_ratio)
                index = lazy_where(times, '[]', time_range)
                theta_mean_ratio_slide = theta_mean_ratio_slide[index]
                times = times[index]
                store_data, theta_mean_ratio_slide_var, times, theta_mean_ratio_slide
                add_setting, theta_mean_ratio_slide_var, /smart, {$
                    display_type: 'scalar', $
                    short_name: tex2str('theta')+' mean ratio slide', $
                    boxcar_window: boxcar_window, $
                    unit: 'deg'}
            endif


        ;---Combine.
            theta_combo_var = prefix+'theta_combo'
            vars = [theta_var, theta_median_var, theta_mean_ratio_var, theta_mean_ratio_slide_var]
            nvar = n_elements(vars)
            labels = ['original','median','mean '+sgnum2str(boxcar_ratio*100)+'%','mean slide '+sgnum2str(boxcar_ratio*100)+'%']
            colors = sgcolor(['silver','purple','green','red'])
            times = make_bins(time_range, time_step)
            ntime = n_elements(times)
            combo = fltarr(ntime,nvar)
            for ii=0, nvar-1 do begin
                the_var = vars[ii]
                get_data, the_var, the_times, data
                if n_elements(the_times) eq ntime then begin
                    combo[*,ii] = data
                endif else begin
                    combo[*,ii] = interpol(data, the_times, times, /quadratic)
                    ;combo[*,ii] = interpol(data, the_times, times)
                endelse
            endfor
            store_data, theta_combo_var, times, combo
            add_setting, theta_combo_var, {$
                ytitle: '(deg)', $
                labels: labels, $
                colors: colors}

            if keyword_set(test_time) then begin
                tplot, theta_combo_var, trange=time_range
            endif

        ;---Calculate parameters to gauge the differences.
            mean_median_study = dictionary()
            for ii=1, nvar-1 do mean_median_study[vars[ii]] = dictionary('name', vars[ii])
            name_length = max(strlen(vars[1:*]))

            ; Deviation from original data: stddev(x-x0)/stddev(x0)
            stddevs = fltarr(nvar-1)
            stddev0 = stddev(combo[*,0])
            for ii=1, nvar-1 do begin
                tinfo = mean_median_study[vars[ii]]
                tinfo['deviation'] = stddev(combo[*,ii]-combo[*,0])/stddev0
            endfor

            ; The smoothness: standard deviation of derivative.
            stddev0 = stddev(deriv(combo[*,0]))
            for ii=1, nvar-1 do begin
                tinfo = mean_median_study[vars[ii]]
                tinfo['smoothness'] = stddev(deriv(combo[*,ii]))/stddev0
            endfor
            candidate['mean_median_study'] = mean_median_study


        ;---Print results to file.
            msg = string(candidate.id,format='(I4)')+tab+$
                strjoin(time_string(candidate.time_range,tformat='YYYY-MM-DD/hh:mm'),' ')+tab+$
                candidate.region+tab+$
                candidate.search_type+tab+$
                string(candidate.duration,format='(I4)')+tab+$
                strjoin(candidate.all_probes,',')
            lprmsg, msg, log_file

            foreach tinfo, mean_median_study do begin
                msg = tab+extend_string(tinfo.name,length=name_length)+tab+'deviation/smoothness'+tab+string(tinfo.deviation,format='(F6.3)')+'/'+string(tinfo.smoothness,format='(F6.3)')
                lprmsg, msg, log_file
            endforeach
        endforeach
    endif else begin
        lines = read_all_lines(log_file)
        ntest_key = 3
        ntest_info = n_elements(lines)/(ntest_key+1)
        test_infos = list(length=ntest_info)
        for ii=0, ntest_info-1 do begin
            i0 = ii*(ntest_key+1)
            tline = lines[i0+0]
            infos = strsplit(tline, ' ', /extract)
            ; The major line: candidate id | start_time | end_time | region | search | duration | nsection | probes
            candidate_id = fix(infos[0])
            time_range = time_double(infos[1:2])
            region_name = infos[3]
            search_name = infos[4]
            duration = float(infos[5])
            all_probes = strsplit(infos[6],',',/extract)
            the_info = dictionary($
                'id', candidate_id, $
                'time_range', time_range, $
                'duration', duration, $  ; in min.
                'search_type', search_name, $
                'region', region_name, $
                'all_probes', all_probes)
            tline = lines[i0+1]
            infos = strsplit(tline, ' ', /extract)
            test_probe = strmid(infos[0],0,strpos(infos[0],'_'))
            test_keys = strarr(ntest_key)
            test_deviation = fltarr(ntest_key)
            test_smoothness = fltarr(ntest_key)
            for jj=0,ntest_key-1 do begin
                tline = lines[i0+jj+1]
                infos = strsplit(tline, ' /', /extract)
                test_keys[jj] = infos[0]
                test_deviation[jj] = float(infos[3])
                test_smoothness[jj] = float(infos[4])
            endfor
            the_info['test_probe'] = test_probe
            the_info['test_keys'] = test_keys
            the_info['test_deviation'] = test_deviation
            the_info['test_smoothness'] = test_smoothness

            test_infos[ii] = the_info
        endfor


        all_deviation = fltarr(ntest_info,ntest_key)
        all_smoothness = fltarr(ntest_info,ntest_key)
        for ii=0, ntest_info-1 do begin
            all_deviation[ii,*] = test_infos[ii].test_deviation
            all_smoothness[ii,*] = test_infos[ii].test_smoothness
        endfor
        test_probe = (test_infos[0]).test_probe
        test_keys = (test_infos[0]).test_keys
        test_labels = strmid(test_keys,12)
        index_mean_ratio = where(test_keys eq test_probe+'_theta_mean_ratio')
        index_mean_ratio_slide = where(test_keys eq test_probe+'_theta_mean_ratio_slide')
        index_median = where(test_keys eq test_probe+'_theta_median')


    ;---Plot.
        nterms = 3
        xrange = 1+[-1,1]*0.4
        ytitle = 'Percent (%)'
        fwhm_coef = sqrt(-alog(0.5)*2)
        binsize = 0.02

        xxs_list = list()
        yys_list = list()
        xtitle_list = list()

    ;---deviation.
        mean_ratio_vs_mean_ratio_slide = all_deviation[*,index_mean_ratio]/all_deviation[*,index_mean_ratio_slide]
        yys_list.add, histogram(mean_ratio_vs_mean_ratio_slide, locations=xxs, binsize=binsize)*100./ntest_info
        xxs_list.add, xxs
        xtitle_list.add, 'Deviation. '+test_labels[index_mean_ratio]+'/'+test_labels[index_mean_ratio_slide]

        mean_ratio_vs_median = all_deviation[*,index_mean_ratio]/all_deviation[*,index_median]
        yys_list.add, histogram(mean_ratio_vs_median, locations=xxs, binsize=binsize)*100./ntest_info
        xxs_list.add, xxs
        xtitle_list.add, 'Deviation. '+test_labels[index_mean_ratio]+'/'+test_labels[index_median]

    ;---smoothness.
        mean_ratio_vs_mean_ratio_slide = all_smoothness[*,index_mean_ratio]/all_smoothness[*,index_mean_ratio_slide]
        yys_list.add, histogram(mean_ratio_vs_mean_ratio_slide, locations=xxs, binsize=binsize)*100./ntest_info
        xxs_list.add, xxs
        xtitle_list.add, 'Smoothness. '+test_labels[index_mean_ratio]+'/'+test_labels[index_mean_ratio_slide]

        mean_ratio_vs_median = all_smoothness[*,index_mean_ratio]/all_smoothness[*,index_median]
        yys_list.add, histogram(mean_ratio_vs_median, locations=xxs, binsize=binsize)*100./ntest_info
        xxs_list.add, xxs
        xtitle_list.add, 'Smoothness. '+test_labels[index_mean_ratio]+'/'+test_labels[index_median]


    ;---Plot.
        plot_file = join_path([project.plot_dir,'diagnostic_plot','study_smooth_theta',search_name+'_'+test_probe+'.pdf'])
        if keyword_set(test) then plot_file = test
        sgopen, plot_file, xsize=6, ysize=6, /inch
        margins = [8,5,2,3]
        poss = sgcalcpos(2,2, margins=margins, xchsz=xchsz, ychsz=ychsz, ypad=5, xpad=8)
        xticklen = -0.01
        yticklen = -0.01

        msg = 'Study smoothing '+strupcase(test_probe)+' '+tex2str('theta')+' using: '+strjoin(test_labels,', ')
        tx = poss[0,0,0]
        ty = 1-ychsz*constant('full_ychsz')*2
        xyouts, tx,ty,/normal, msg
        msg = '# of events studied: '+string(ntest_info,format='(I0)')+', search_name: '+search_name
        ty = 1-ychsz*(constant('full_ychsz')*3+constant('lineskip'))
        xyouts, tx,ty,/normal, msg

        fig_labels = letters(4)+'.'
        foreach xtitle, xtitle_list, ii do begin
            index = array_indices([2,2], ii, /dimensions)
            tpos = poss[*,index[1],index[0]]
            xxs = xxs_list[ii]
            yys = yys_list[ii]
            yfit = gaussfit(xxs, yys, coeff, nterms=nterms)
            yrange = sg_autolim([0,max(yys)])
            the_color = sgcolor('red')

            plot, xxs, yys, psym=1, symsize=0.5, position=tpos, /noerase, $
                xstyle=1, xtitle=xtitle, xrange=xrange, xticklen=xticklen, $
                ystyle=1, ytitle=ytitle, yrange=yrange, yticklen=yticklen
            plots, xrange, [0,0], linestyle=1
            oplot, xxs, yfit, color=the_color
            plots, coeff[1]+[0,0], [0,coeff[0]], linestyle=1
            plots, coeff[1]+[-1,1]*coeff[2]*fwhm_coef, coeff[0]*0.5+[0,0], linestyle=1

            msg = sgnum2str(coeff[1],ndec=2)+tex2str('pm')+sgnum2str(coeff[2]*fwhm_coef,ndec=2)
            tmp = convert_coord(coeff[1],coeff[0]*0.5, /data, /to_normal)
            tx = tmp[0]
            ty = tpos[3]-ychsz*constant('full_ychsz')
            xyouts, tx,ty,/normal,alignment=0.5, msg, charsize=constant('label_size'), color=the_color

            tx = tpos[0]-xchsz*5
            ty = tpos[3]-ychsz*constant('full_ychsz')
            xyouts, tx,ty,/normal, fig_labels[ii]
        endforeach
        if keyword_set(test) then stop
        sgclose
    endelse


end

project = azim_df_load_project()
tests = list()
tests.add, dictionary('probe', 'rbspb', 'search_name', 'within_15Re')
tests.add, dictionary('probe', 'tha', 'search_name', 'beyond_15Re')
foreach test, tests do begin
    probe = test.probe
    search_name = test.search_name
    azim_df_study_smooth_theta, project=project, probe=probe, search_name=search_name
endforeach
end
