;+
; The Edot0 data around dawn are very large. Check why here.
;-

    project = pflux_survey_load_project()
    probes = ['a','b']
    plot_file = join_path([project.plot_dir,'figures','fig_de_vs_dedot0_vs_mlt.pdf'])
    if keyword_set(test) then plot_file = test
    mission_time_range = project.time_range

    de_vars = ['de','dedot0']
    fac_labels = ['b','w','o']
    fac_full_labels = ['parallel','westward','outward']
    ndim = 3
    de_step = 16*60.    ; To select data every minute.
    secofday = constant('secofday')

    ; Settings for edot0_flag.
    min_bw_ratio = 0.2
    edot0_flag_pad_time = 1800.


;---Check if data for plotting are saved.
    save_vars = list()
    data_file = join_path([project.data_dir,'pflux_survey_de_vs_dedot0.tplot'])
    ;tmp = pflux_grant_load_project()
    ;data_file2 = join_path([tmp.data_dir,'pflux_survey_de_vs_dedot0.tplot'])
    if file_test(data_file) eq 1 then tplot_restore, filename=data_file
    file_updated = file_test(data_file) eq 0

    foreach probe, probes do begin
        rbspx = 'rbsp'+probe
        prefix = rbspx+'_'

    ;---Load orbit data.
        r_gsm_var = prefix+'r_gsm'
        save_vars.add, r_gsm_var
        if check_if_update(r_gsm_var) then begin
            rbsp_read_orbit, mission_time_range, probe=probe
            file_updated = 1
        endif

    ;---Load E field.
        get_data, r_gsm_var, common_times
        ntime = n_elements(common_times)
        secofday = constant('secofday')
        nday = total(mission_time_range*[-1,1])/secofday
        days = smkarthm(mission_time_range[0],secofday,nday, 'x0')
        foreach de_var, de_vars do begin
            de_var = prefix+de_var+'_fac'
            de_save_var = de_var+'_interp'
            save_vars.add, de_save_var
            if ~check_if_update(de_save_var) then continue
            file_updated = 1
            de_data = fltarr(ntime,ndim)
            for day_id=0, nday-1 do begin
                time_range = days[day_id]+[0,secofday]
                pflux_grant_read_preprocessed_ebfield, time_range, probe=probe
                index = lazy_where(common_times,'[)', time_range)
                de_data[index,*] = (get_var_data(de_var))[0:*:de_step,*]
            endfor
            store_data, de_save_var, common_times, de_data
            short_name = (de_var eq 'de')? 'E': 'E!S!Udot0!R!N'
            add_setting, de_save_var, /smart, {$
                display_type: 'vector', $
                unit: 'mV/m', $
                short_name: short_name, $
                coord: '', $
                coord_labels: fac_labels}
        endforeach

    ;---Load Bw ratio.
        bw_var = prefix+'bw_ratio'
        bw_save_var = bw_var+'_interp'
        save_vars.add, bw_save_var
        if check_if_update(bw_save_var) then begin
            bw_data = fltarr(ntime)
            for day_id=0,nday-1 do begin
                time_range = days[day_id]+[0,secofday]
                pflux_grant_read_preprocessed_ebfield, time_range, probe=probe, id='bw_ratio'
                index = lazy_where(common_times,'[)', time_range)
                bw_data[index,*] = (get_var_data(bw_var))[0:*:de_step]
            endfor
            file_updated = 1

            store_data, bw_save_var, common_times, bw_data
            add_setting, bw_save_var, /smart, {$
                display_type: 'scalar', $
                unit: '#', $
                short_name: 'Bw/|B|'}
        endif

    ;---Index for good Edot0 data.
        edot0_flag_var = prefix+'edot0_flag'
        save_vars.add, edot0_flag_var
        if check_if_update(edot0_flag_var) then begin
            bw_ratio = get_var_data(bw_save_var)
            edot0_flag = abs(bw_ratio) le min_bw_ratio  ; 1 for bad data.
            index = where(edot0_flag eq 1, count)
            time_step = total(common_times[0:1]*[-1,1])
            drec = edot0_flag_pad_time/time_step
            ntime = n_elements(common_times)
            for ii=0,count-1 do begin
                i0 = (index[ii]-drec)>0
                i1 = (index[ii]+drec)<(ntime-1)
                edot0_flag[i0:i1] = 1
            endfor

            file_updated = 1
            store_data, edot0_flag_var, common_times, edot0_flag
            add_setting, edot0_flag_var, /smart, {$
                display_type: 'scalar', $
                unit: '#', $
                short_name: 'Edot0 flag', $
                min_bw_ratio: min_bw_ratio, $
                pad_time: edot0_flag_pad_time, $
                fieldnam: '1 for bad data'}
        endif
    endforeach

    save_vars = save_vars.toarray()
    if file_updated then tplot_save, save_vars, filename=data_file


;---Remove data when bw_ratio does not meet criteria.
    de_range = list()
    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        get_data, prefix+'edot0_flag', common_times, dedot0_flag
        dedot0_fac = get_var_data(prefix+'dedot0_fac_interp')
        de_fac = get_var_data(prefix+'de_fac_interp')
        r_gsm = get_var_data(prefix+'r_gsm')
        de_range.add, minmax(de_fac)
        de_range.add, minmax(dedot0_fac)

        if ~check_if_update(prefix+'de_fac') then continue
        index = where(dedot0_flag eq 0)
        times = common_times[index]
        dedot0_fac = dedot0_fac[index,*]
        de_fac = de_fac[index,*]
        r_sm = cotran(r_gsm[index,*], common_times, 'gsm2sm')
        mlt = pseudo_mlt(r_sm)
        store_data, prefix+'de_fac', times, de_fac
        store_data, prefix+'dedot0_fac', times, dedot0_fac
        store_data, prefix+'mlt', times, mlt
    endforeach
    de_range = minmax(de_range.toarray())
    de_tick_step = 50.
    de_range = [-1,1]*ceil(max(abs(de_range))/de_tick_step)*de_tick_step
    de_tickv = make_bins(de_range,de_tick_step*2,/inner)
    de_range = [-1,1]*180
    de_tickv = make_bins(de_range,de_tick_step*2,/inner)
    de_ticks = n_elements(de_tickv)-1
    de_minor = 5
    de_unit = 'mV/m'
    de_log = 0
    short_name = 'E'
    de_labels = fac_full_labels
    de_ticklen = -0.02
    fig_letters = letters(ndim)


;---MLT settings.
    mlt_range = [-6,6]
    mlt_step = 1
    mlt_tick_step = 6
    mlt_centers = make_bins(mlt_range,mlt_step)
    mlt_vertices = [[mlt_centers-mlt_step*0.5],[mlt_centers+mlt_step*0.5]]
    nmlt_bin = n_elements(mlt_centers)

    ct = 6  ; circular color table.
    ncolor = n_elements(mlt_centers)
    color_min = 0+60
    color_max = 255-60
    index_colors = smkarthm(color_min, color_max, ncolor, 'n')
    colors = index_colors
    foreach color, index_colors, ii do colors[ii] = sgcolor(color, ct=ct)

end
