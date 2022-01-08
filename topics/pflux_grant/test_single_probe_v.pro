;+
; Test to get Vsc from 4 single probe Vs, and use it to tell which ones work.
;-

test = 0
test_days = list()
if n_elements(project) eq 0 then project = pflux_grant_load_project()

;---RBSP-A.
;    probe = 'a'
;    ; Good days.
;    day_type = 'good_day'
;    test_days.add, dictionary($
;        'probe', probe, 'day_type', day_type, $
;        'date', time_double('2012-10-12'))
;    test_days.add, dictionary($
;        'probe', probe, 'day_type', day_type, $
;        'date', time_double('2014-05-12'))
;    ; Bad days.
;    day_type = 'bad_day'
;    test_days.add, dictionary($
;        'probe', probe, 'day_type', day_type, $
;        'date', time_double('2016-05-12'))  ; V34 bad.
;    test_days.add, dictionary($
;        'probe', probe, 'day_type', day_type, $
;        'date', time_double('2019-05-12'))  ; V4 bad.

    probe = 'a'
    ; Good days.
    day_type = 'good_day'
;    test_days.add, dictionary($
;        'probe', probe, 'day_type', day_type, $
;        'date', time_double('2012-10-12'))
;    test_days.add, dictionary($
;        'probe', probe, 'day_type', day_type, $
;        'date', time_double('2014-05-12')) ; there are offsets, but trends are the same.


    test_days = list()
    day_type = 'monthly_test'
    test_years = string(make_bins([2012,2019],1),format='(I4)')
    test_months = string(make_bins([1,12],1),format='(I02)')
    valid_range = time_double(['2012-10-01','2019-09-01'])
    foreach probe, ['a','b'] do begin
        foreach year, test_years do begin
            foreach month, test_months do begin
                test_day = time_double(year+'-'+month)
                if product(test_day-valid_range) gt 0 then continue
                test_days.add, dictionary($
                    'probe', probe, 'day_type', day_type, $
                    'date', test_day)
            endforeach
        endforeach
    endforeach

;    ; Bad days.
;    day_type = 'bad_day'
;    test_days.add, dictionary($
;        'probe', probe, 'day_type', day_type, $
;        'date', time_double('2016-05-12')) ; V1 bad.
;    test_days.add, dictionary($
;        'probe', probe, 'day_type', day_type, $
;        'date', time_double('2018-05-12')) ; all bad.


    secofday = constant('secofday')
    max_valid_v = 200.
    max_good_v = 150.  ; By inspecting monthly plots, this seems to be a good threshold for bad |V|.
    efield_time_step = 1d/16
    spinfit_time_step = 10d
    v_colors = sgcolor(['red','green','blue','black'])
    v_labels = 'V'+['1','2','3','4']
    tplot_options, 'labflag', -1

    foreach test_day, test_days do begin
        date = test_day.date
        date_time_range = date+[0,secofday]
        probe = test_day.probe
        prefix = 'rbsp'+probe+'_'
        day_type = test_day.day_type
        plot_file = join_path([project.plot_dir,'test_vsc',day_type,$
            prefix+'test_vsc_'+time_string(date,tformat='YYYY_MMDD')+'_v01.pdf'])
        if keyword_set(test) then plot_file = 0
        ;if file_test(plot_file) eq 1 then continue

    ;---Load Vsvy and orbit.
        datatype = (date ge time_double('2016-02-28'))? 'l2%vsvy-highres2': 'l2%vsvy-highres'
        rbsp_read_efw, date_time_range, id=datatype, probe=probe
        rbsp_read_orbit, date_time_range, probe=probe

        r_var = prefix+'r_gsm'
        dis = snorm(get_var_data(r_var, times=orbit_times))

        vsvy_var = prefix+'vsvy'
        rename_var, 'vsvy', to=vsvy_var
        get_data, vsvy_var, times, vsvy
        if n_elements(vsvy) le 6 then continue
        vsvy = vsvy[*,0:3]

        min_dis = 3.
        dis = interpol(dis, orbit_times, times)
        index = where(dis le min_dis, count)
        if count ne 0 then vsvy[index,*] = !values.f_nan

        index = where(abs(vsvy) ge max_valid_v, count)
        if count ne 0 then vsvy[index] = !values.f_nan
        store_data, vsvy_var, times, vsvy, limits={$
            ytitle: '(V)', $
            colors: v_colors, $
            labels: v_labels}
        uniform_time, vsvy_var, spinfit_time_step
        

    ;---Analysis Vsc.
        vsvy = get_var_data(prefix+'vsvy', times=common_times)
        index = where(abs(vsvy) ge max_good_v, count)
        if count ne 0 then vsvy[index] = !values.f_nan
        ncommon_time = n_elements(common_times)
        vsc_mean = fltarr(ncommon_time)
        vsc_median = fltarr(ncommon_time)
        for ii=0,ncommon_time-1 do begin
            the_vs = reform(vsvy[ii,*])
            the_vs = the_vs[sort(the_vs)]
            vsc_mean[ii] = mean(the_vs[1:2])
            vsc_median[ii] = median(the_vs[0:2])    ; remove the largest, median out of 3 is better than out of 4.
;            index = where(finite(the_vs), count)
;            if count lt 2 then begin
;                vsc_mean[ii] = !values.f_nan
;                vsc_median[ii] = !values.f_nan
;            endif
        endfor
        vsc_mean_var = prefix+'vsc_mean'
        store_data, vsc_mean_var, common_times, vsc_mean, limits={$
            ytitle: '(V)', $
            labels: 'Vsc mean'}

        for ii=0,ncommon_time-1 do begin
            the_vs = reform(vsvy[ii,*])
            the_vs = the_vs[sort(the_vs)]
            vsc_median[ii] = median(the_vs)
        endfor
        vsc_median_var = prefix+'vsc_median'
        store_data, vsc_median_var, common_times, vsc_median, limits={$
            ytitle: '(V)', $
            labels: 'Vsc median'}

        vsc_diff_var = prefix+'vsc_diff'
        store_data, vsc_diff_var, common_times, vsc_mean-vsc_median, limits={$
            ytitle: '(V)', $
            yrange: [-2,0.5], $
            labels: 'Vsc!C  mean-median'}

        sgopen, plot_file, xsize=5, ysize=5
        tplot, prefix+['vsvy','vsc_mean','vsc_median','vsc_diff'], trange=date_time_range
        if keyword_set(test) then stop
        sgclose
    endforeach

end
