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


;    test_days = list()
;    day_type = 'monthly_test'
;    test_years = string(make_bins([2012,2019],1),format='(I4)')
;    test_months = string(make_bins([1,12],1),format='(I02)')
;    valid_range = time_double(['2012-10-01','2019-09-01'])
;    foreach probe, ['a','b'] do begin
;        foreach year, test_years do begin
;            foreach month, test_months do begin
;                test_day = time_double(year+'-'+month)
;                if product(test_day-valid_range) gt 0 then continue
;                test_days.add, dictionary($
;                    'probe', probe, 'day_type', day_type, $
;                    'date', test_day)
;            endforeach
;        endforeach
;    endforeach


test_days = list()
valid_range = time_double(['2012-10-01','2019-09-01'])
dates = make_bins(valid_range, constant('secofday'))
foreach probe, ['a','b'] do begin
    foreach date, dates do begin
        year = time_string(date,tformat='YYYY')
        day_type = 'rbsp'+probe+'/'+year
        test_days.add, dictionary($
            'probe', probe, 'day_type', day_type, $
            'date', date)
    endforeach
endforeach


;    test_days = list()
;    test_days.add, dictionary($
;        'probe', 'a', 'day_type', 'known_days', $
;        'date', time_double('2013-03-01'))
;    test_days.add, dictionary($
;        'probe', 'b', 'day_type', 'known_days', $
;        'date', time_double('2013-06-07'))
;    test_days.add, dictionary($
;        'probe', 'a', 'day_type', 'known_days', $
;        'date', time_double('2012-11-01'))

    secofday = constant('secofday')
    max_valid_v = 200.
    max_good_v = 150.  ; By inspecting monthly plots, this seems to be a good threshold for bad |V|.
    efield_time_step = 1d/16
    spinfit_time_step = 10d
    spin_period = rbsp_info('spin_period')
    v_colors = sgcolor(['red','green','blue','black'])
    v_labels = 'V'+['1','2','3','4']
    tplot_options, 'labflag', -1
    half_boom_length = 50d

    foreach test_day, test_days do begin
        date = test_day.date
        date_time_range = date+[0,secofday]
        probe = test_day.probe
        prefix = 'rbsp'+probe+'_'
        day_type = test_day.day_type
        plot_file = join_path([project.plot_dir,'test_select_good_probes',day_type,$
            prefix+'test_vsc_diff_'+time_string(date,tformat='YYYY_MMDD')+'_v01.pdf'])
        if keyword_set(test) then plot_file = 0
        ;if file_test(plot_file) eq 1 then continue

    ;---Load Vsvy.
        datatype = (date ge time_double('2016-02-28'))? 'l2%vsvy-highres2': 'l2%vsvy-highres'
        rbsp_read_efw, date_time_range, id=datatype, probe=probe

        vsvy_var = prefix+'vsvy'
        rename_var, 'vsvy', to=vsvy_var
        get_data, vsvy_var, times, vsvy
        if n_elements(vsvy) le 6 then continue
        vsvy = vsvy[*,0:3]

        index = where(abs(vsvy) ge max_valid_v, count)
        if count ne 0 then vsvy[index] = !values.f_nan
        store_data, vsvy_var, times, vsvy, limits={$
            ytitle: '(V)', $
            colors: v_colors, $
            labels: v_labels}
        uniform_time, vsvy_var, efield_time_step
        vsvy = get_var_data(vsvy_var, times=common_times)

;    ;---Remove perigee data.
        min_dis = 3.
        rbsp_read_orbit, date_time_range, probe=probe
        r_var = prefix+'r_gsm'
        dis = snorm(get_var_data(r_var, times=orbit_times))
        dis = interpol(dis, orbit_times, common_times)
        index = where(dis le min_dis, count)
        if count ne 0 then vsvy[index,*] = !values.f_nan

    ;---Remove SDT, eclipse, maneuver.
        rbsp_read_eclipse_flag, date_time_range, probe=probe
        rbsp_read_sdt_flag, date_time_range, probe=probe

        flag_vars = prefix+['eclipse','sdt']+'_flag'
        flag_time_step = 60.    ; sec.
        pad_time = 120.
        foreach tvar, flag_vars do begin
            flags = get_var_data(tvar, times=times)
            index = where(flags eq 1, count)
            if count eq 0 then continue
            flag_time_ranges = time_to_range(times[index], time_step=flag_time_step)
            nflag_time_range = n_elements(flag_time_ranges)/2
            for ii=0, nflag_time_range-1 do begin
                index = lazy_where(common_times, '[]', flag_time_ranges[ii,*]+[-1,1]*pad_time, count=count)
                if count eq 0 then continue
                vsvy[index,*] = !values.f_nan
            endfor
        endforeach

        maneuver_times = rbsp_read_maneuver_time(date_time_range, probe=probe)
        foreach maneuver_time, maneuver_times do begin
            ; Do nothing.
        endforeach
        store_data, vsvy_var, common_times, vsvy


    ;---Low-res version.
        vsvy_lowres_var = prefix+'vsvy_lowres'
        copy_data, vsvy_var, vsvy_lowres_var
        uniform_time, vsvy_lowres_var, spinfit_time_step
        vsvy = get_var_data(vsvy_lowres_var, times=common_times)


    ;---Analysis Vsc.
        vsvy = get_var_data(vsvy_lowres_var, times=lowres_times)
        index = where(abs(vsvy) ge max_good_v, count)
        if count ne 0 then vsvy[index] = !values.f_nan
        nlowres_time = n_elements(lowres_times)
        vsc_mean = fltarr(nlowres_time)
        vsc_median = fltarr(nlowres_time)
        for ii=0,nlowres_time-1 do begin
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
        store_data, vsc_mean_var, lowres_times, vsc_mean, limits={$
            ytitle: '(V)', $
            labels: 'Vsc mean'}

        for ii=0,nlowres_time-1 do begin
            the_vs = reform(vsvy[ii,*])
            the_vs = the_vs[sort(the_vs)]
            vsc_median[ii] = median(the_vs)
        endfor
        vsc_median_var = prefix+'vsc_median'
        store_data, vsc_median_var, lowres_times, vsc_median, limits={$
            ytitle: '(V)', $
            labels: 'Vsc median'}


    ;---Get the E field. "Single-ended E field" does not seem correct.
        vsvy = get_var_data(vsvy_var, times=highres_times)
        vsc = get_var_data(vsc_median_var, at=highres_times)
        smooth_width = round(spin_period/efield_time_step)
        for ii=0,3 do begin
            id_str = string(ii+1,format='(I0)')
            dv = vsvy[*,ii]-vsc
            dv0 = smooth(dv, smooth_width, /edge_truncate, /nan)
            dv = dv-dv0
            dv0 = interpol(dv0, highres_times, lowres_times)
            store_data, prefix+'dv0_'+id_str, lowres_times, dv0, limits={$
                ytitle:'(V)', $
                labels:'(V'+id_str+'-Vsc)_BG', $
                ystyle: 1, $
                yrange: [-1,1]*5}
            de = dv/half_boom_length*1e3
            store_data, prefix+'e'+id_str, highres_times, de, limits={$
                ytitle:'(mV)', $
                labels:'E'+id_str, $
                ystyle: 1, $
                yrange: [-1,1]*200}
        endfor

    ;---Get the flags for probes.
        nprobe = 4
        max_valid_dv0 = 5.  ; V.
        max_valid_dv0_ratio = 0.2
        orbit_time_ranges = pflux_grant_calculate_orbit_time_range(date_time_range+[-1,1]*secofday, probe=probe)
        index = where(orbit_time_ranges[*,0] le date_time_range[0])
        orbit_time_ranges = orbit_time_ranges[index[-1]:*,*]
        index = where(orbit_time_ranges[*,1] ge date_time_range[1])
        orbit_time_ranges = orbit_time_ranges[0:index[0],*]
        norbit = n_elements(orbit_time_ranges)/2
        for ii=0, nprobe-1 do begin
            id_str = string(ii+1,format='(I0)')
            dv0_var = prefix+'dv0_'+id_str
            probe_flags = intarr(nlowres_time)
            dv0 = get_var_data(dv0_var, times=lowres_times)
            index = where(abs(dv0) lt max_valid_dv0, count)
            if count ne 0 then probe_flags[index] = 1

            flag_var = prefix+'v'+id_str+'_flag'
            store_data, flag_var, lowres_times, probe_flags, limits={$
                ytitle: '(#)', $
                labels: 'V'+id_str+' flag!C  1: good', $
                ystyle: 1, $
                yrange: [0,1]+[-1,1]*0.2, $
                ytickv: [0,1], $
                panel_size: 0.4, $
                yticks: 1, $
                yminor: 0}

;            for jj=0, norbit-1 do begin
;                current_time_range = reform(orbit_time_ranges[jj,*])
;                the_flags = get_var_data(flag_var, in=current_time_range, times=times)
;                ntime = n_elements(times)
;                index = where(the_flags eq 1, count)
;                if count le ntime*max_valid_dv0_ratio then continue
;                index = lazy_where(lowres_times, '[]', current_time_range)
;                probe_flags[index] = 1
;                store_data, flag_var, lowres_times, probe_flags
;            endfor
        endfor

        sgopen, plot_file, xsize=6, ysize=10
        tplot, prefix+['vsvy','vsc_median','dv0_?','v?_flag'], trange=date_time_range, title='RBSP-'+strupcase(probe)+' '+time_string(date, tformat='YYYY_MMDD')
        maneuver_times = rbsp_read_maneuver_time(date_time_range, probe=probe)
        if n_elements(maneuver_times) ne 0 then begin
            timebar, maneuver_times, color=sgcolor('red')
        endif
        if keyword_set(test) then stop
        sgclose
    endforeach

end
