;+
; Study the different between E and Edot0.
;-

test = 0

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project['pflux_calc_setting']
    probes = ['a','b']
    plot_dir = join_path([project.plot_dir,'pflux_survey','figures'])
    orbit_time_step = 60.

    de_vars = ['de','dedot0']
    fac_labels = ['earth','west','out']

    save_vars = list()
    data_file = join_path([project.data_dir,'pflux_survey_de_vs_dedot0.tplot'])
    if file_test(data_file) eq 1 then tplot_restore, filename=data_file

;---Read data for each probe.
    foreach probe, probes do begin
        rbspx = 'rbsp'+probe
        prefix = rbspx+'_'
        mission_time_range = pflux_calc_setting['rbspa'].time_range

    ;---Load orbit data.
        r_gsm_var = prefix+'r_gsm'
        save_vars.add, r_gsm_var
        if check_if_update(r_gsm_var) then begin
            rbsp_read_orbit, mission_time_range, probe=probe
        endif

    ;---Load E field.
        get_data, r_gsm_var, common_times
        ntime = n_elements(common_times)
        ndim = 3
        de_step = 16*60.
        secofday = constant('secofday')
        nday = total(mission_time_range*[-1,1])/secofday
        days = smkarthm(mission_time_range[0],secofday,nday, 'x0')
        foreach de_var, de_vars do begin
            de_var = prefix+de_var+'_fac'
            de_save_var = de_var+'_interp'
            save_vars.add, de_save_var
            if ~check_if_update(de_save_var) then continue
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
    endforeach

    save_vars = save_vars.toarray()
    if file_test(data_file) eq 0 then tplot_save, save_vars, filename=data_file


;---Pick out times when Edot0 are not NaNs.
    de_unit = 'mV/m'
    de_log = 0
    de_labels = ['earth','west','out']+'ward'
    de_range = [-1,1]*300
    xtickv = [-1,0,1]*200
    xticks = n_elements(xtickv)-1
    xminor = 4
    ytickv = xtickv
    yticks = xticks
    yminor = xminor
    ndim = 3
    fig_letters = letters(ndim)
    de_ticklen = -0.02

    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        get_data, prefix+'dedot0_fac_interp', times, dedot0
        get_data, prefix+'de_fac_interp', times, de
        index = where(finite(dedot0[*,0]))
        times = times[index]
        de = de[index,*]
        dedot0 = dedot0[index,*]
        mlt = azim_df_calc_pseudo_mlt(get_var_data(prefix+'r_gsm', at=times))

        fig_xsize = 6
        fig_ysize = 3
        plot_file = join_path([project.plot_dir,'pflux_survey','figures','fig_de_vs_dedot0_rbsp'+probe+'.pdf'])
        if keyword_set(test) then plot_file = test
        sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize

        poss = sgcalcpos(1,3, xchsz=xchsz, ychsz=ychsz, xpad=1)
        for ii=0,ndim-1 do begin
            tpos = poss[*,ii]
            txs = de[*,ii]
            tys = dedot0[*,ii]
            xtitle = 'E ('+de_unit+')'
            ytitle = 'E!Ddot0!N ('+de_unit+')'
            if ii ne 0 then ytitle = ' '
            if ii ne 0 then ytickformat='(A1)' else ytickformat=''
            plot, txs,tys, psym=1, symsize=0.5, /iso, $
                xstyle=1, xrange=de_range, xlog=de_log, xtitle=xtitle, xticklen=de_ticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
                ystyle=1, yrange=de_range, ylog=de_log, ytitle=ytitle, yticklen=de_ticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, $
                position=tpos, /noerase, ytickformat=ytickformat
            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            fig_label = fig_letters[ii]+'. FAC E '+de_labels[ii]
            plots, de_range, [0,0], linestyle=1
            plots, [0,0], de_range, linestyle=1
            plots, de_range, de_range, linestyle=1
            xyouts, tx,ty,/normal, fig_label
        endfor

        if keyword_set(test) then stop
        sgclose

    endforeach

end
