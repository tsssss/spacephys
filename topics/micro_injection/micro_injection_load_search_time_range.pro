;+
; Load time range for event search.
;-
function micro_injection_load_search_time_range, test=test, root_dir=root_dir


    if n_elements(root_dir) eq 0 then root_dir = srootdir()
    plot_dir = root_dir
    data_dir = root_dir

;---Settings.
    search_tr = time_double(['2015','2018'])
    mlt_range = [6d,18]
    dis_range = [0d,12.5]

    base = 'micro_injection_load_search_time_range_v01.pdf'
    plot_file = join_path([plot_dir,base])
    base = 'micro_injection_load_search_time_range_v01.cdf'
    data_file = join_path([data_dir,base])

    search_time_range_var = 'search_time_ranges'
    probe = '1'
    prefix = 'mms'+probe+'_'
    apogee_mlt_var = prefix+'apogee_mlt'
    apogee_dis_var = prefix+'apogee_dis'
    if file_test(data_file) eq 0 then begin
    ;---Load orbit over many years and obtain apogee history.
        mission_probe = 'mms'+probe
        ;r_gsm_var = lets_read_this(func='mms_read_orbit', search_tr, probe=mission_probe)
        find_apogee, r_gsm_var, period=period, apogee=apogee, times=apogee_times
        mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)
        mlt_var = mlat_vars.mlt
        apogee_mlts = get_var_data(mlt_var, at=apogee_times)
        
        ; conver to [0,24].
        index = where(apogee_mlts lt 0, count)
        if count ne 0 then apogee_mlts[index] += 24
        store_data, apogee_mlt_var, apogee_times, apogee_mlts
        add_setting, apogee_mlt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'unit', 'h', $
            'short_name', strupcase(mission_probe), $
            'constant', [-6,0,6]+12 )
        set_ytick, apogee_mlt_var, yrange=[-1,1]*12+12, ystep=6, yminor=6
        apogee_dis_var = rename_var(prefix+'r_gsm_apogee', output=apogee_dis_var)
        add_setting, apogee_dis_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'unit', 'Re', $
            'short_name', strupcase(mission_probe), $
            'colors', sgcolor('black'), $
            'constant', [10,15,20,25,30] )
        set_ytick, apogee_dis_var, yrange=[10,30], ystep=10, yminor=5
        
    ;---Determine the search range.
        apogee_mlts = get_var_data(apogee_mlt_var, times=apogee_times)
        apogee_diss = get_var_data(apogee_dis_var)
        index = where((apogee_mlts ge mlt_range[0] and apogee_mlts le mlt_range[1]) and apogee_diss le dis_range[1], count)
        if count eq 0 then message, 'No time found ...'
        search_time_ranges = apogee_times[time_to_range(index,time_step=1)]
        nsearch_time_range = n_elements(search_time_ranges[*,0])
        cdf_save_var, search_time_range_var, value=search_time_ranges, save_as_one=1, filename=data_file
;    endif
    search_time_ranges = cdf_read_var(search_time_range_var, filename=data_file)
    nsearch_time_range = n_elements(search_time_ranges[*,0])
    ; make the time range from the start of a day.
    secofday = constant('secofday')
    search_time_ranges = search_time_ranges-(search_time_ranges mod secofday)
    search_time_ranges[*,1] += secofday

;---Make a plot.
    if file_test(plot_file) eq 0 then begin
        if keyword_set(test) then plot_file = 0
        sgopen, plot_file, size=[6,3], xchsz=xchsz, ychsz=ychsz
        label_size = 0.8
        plot_vars = [apogee_mlt_var,apogee_dis_var]
        nplot_var = n_elements(plot_vars)
        fig_labels = letters(nplot_var)+') Apogee!C     '+['MLT','|R|']
        margins = [12,3,5,2]
        poss = sgcalcpos(nplot_var,margins=margins)
        tpos = poss[*,0]
        tpos[1] = poss[1,-1]
        yrange = [0,1]
        xrange = search_tr
        set_axis, position=tpos, xrange=xrange, yrange=yrange
        for ii=0,nsearch_time_range-1 do begin
            txs = reform(search_time_ranges[ii,*])
            tys = yrange
            polyfill, txs[[0,1,1,0,0]], tys[[0,0,1,1,0]], data=1, color=sgcolor('wheat')
            tx = mean(txs)
            ty = yrange[1]
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            tx = tmp[0]
            ty = tmp[1]+ychsz*0.3
            msg = strjoin(time_string(txs,tformat='YYYY-MM-DD'),' to ')
            xyouts, tx,ty,msg, normal=1, alignment=0.5, charsize=label_size
        endfor
        tplot, plot_vars, position=poss, trange=xrange, noerase=1
        for pid=0,nplot_var-1 do begin
            msg = fig_labels[pid]
            tpos = poss[*,pid]
            tx = tpos[0]-xchsz*11
            ty = tpos[3]-ychsz*0.7
            xyouts, tx,ty,normal=1, msg
        endfor
        if keyword_set(test) then stop
        sgclose
    endif
    
    
    return, search_time_ranges
end


search_time_range = micro_injection_load_search_time_range()
print, time_string(search_time_range)
end