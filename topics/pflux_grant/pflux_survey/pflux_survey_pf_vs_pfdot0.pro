;+
; Study the different between pf and pfdot0.
;-

test = 1

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project['pflux_calc_setting']
    probes = ['a','b']
    plot_dir = join_path([project.plot_dir,'pflux_survey','figures'])
    orbit_time_step = 60.

    pf_vars = ['pf','pfdot0']
    fac_labels = ['earth','west','out']

    data_file = join_path([project.data_dir,'pflux_survey_pf_vs_pfdot0.tplot'])
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

    ;---Load pflux.
        get_data, r_gsm_var, common_times
        ntime = n_elements(common_times)
        ndim = 3
        pf_step = 16*60.
        secofday = constant('secofday')
        nday = total(mission_time_range*[-1,1])/secofday
        days = smkarthm(mission_time_range[0],secofday,nday, 'x0')
        foreach pf_var, pf_vars do begin
            pflux_var = prefix+pf_var+'_fac_norm'
            pflux_save_var = pflux_var+'_interp'
            save_vars.add, pflux_save_var
            if ~check_if_update(pflux_save_var) then continue
            pf_data = fltarr(ntime,ndim)
            for day_id=0, nday-1 do begin
                time_range = days[day_id]+[0,secofday]
                pflux_grant_read_preprocessed_pflux, time_range, probe=probe
                index = where_pro(common_times,'[)', time_range)
                pf_data[index,*] = (get_var_data(pflux_var))[0:*:pf_step,*]
            endfor
            store_data, pflux_save_var, common_times, pf_data
            short_name = (pf_var eq 'pf')? 'S': 'S!S!Udot0!R!N'
            add_setting, pflux_save_var, /smart, {$
                display_type: 'vector', $
                unit: 'mW/m!U2!N', $
                short_name: short_name, $
                coord: '', $
                coord_labels: fac_labels}
        endforeach
    endforeach

    save_vars = save_vars.toarray()
    if keyword_set(data_file) eq 0 then tplot_save, save_vars, filename=data_file


;---Pick out times when pfdot0 are not NaNs.
    pf_unit = 'mW/m!U2!N'
    pf_log = 0
    pf_labels = ['earth','west','out']+'ward'
    pf_range = [-1,1]*450
    xtickv = [-2,1,0,1,2]*200
    xticks = n_elements(xtickv)-1
    xminor = 4
    ndim = 3
    fig_letters = letters(ndim)
    pf_ticklen = -0.02

    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        get_data, prefix+'pfdot0_fac_norm_interp', times, pfdot0
        get_data, prefix+'pf_fac_norm_interp', times, pf
        index = where(finite(pfdot0[*,0]))
        times = times[index]
        pf = pf[index,*]
        pfdot0 = pfdot0[index,*]
        mlt = azim_df_calc_pseudo_mlt(get_var_data(prefix+'r_gsm', at=times))

        fig_xsize = 6
        fig_ysize = 3
        plot_file = join_path([project.plot_dir,'pflux_survey','figures','fig_pf_vs_pfdot0_rbsp'+probe+'.pdf'])
        if keyword_set(test) then plot_file = test
        sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize

        poss = sgcalcpos(1,3, xchsz=xchsz, ychsz=ychsz, xpad=1)
        for ii=0,ndim-1 do begin
            tpos = poss[*,ii]
            txs = pf[*,ii]
            tys = pfdot0[*,ii]
            xtitle = 'S ('+pf_unit+')'
            ytitle = 'S!Ddot0!N ('+pf_unit+')'
            if ii ne 0 then ytitle = ' '
            if ii ne 0 then ytickformat='(A1)' else ytickformat=''
            plot, txs,tys, psym=1, symsize=0.5, /iso, $
                xstyle=1, xrange=pf_range, xlog=pf_log, xtitle=xtitle, xticklen=pf_ticklen, xtickv=xtickv, xticks=xtickx, xminor=xminor, $
                ystyle=1, yrange=pf_range, ylog=pf_log, ytitle=ytitle, yticklen=pf_ticklen, $
                position=tpos, /noerase, ytickformat=ytickformat
            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            fig_label = fig_letters[ii]+'. FAC S '+pf_labels[ii]
            plots, pf_range, [0,0], linestyle=1
            plots, [0,0], pf_range, linestyle=1
            plots, pf_range, pf_range, linestyle=1
            xyouts, tx,ty,/normal, fig_label
        endfor

        if keyword_set(test) then stop
        sgclose

    endforeach

end
