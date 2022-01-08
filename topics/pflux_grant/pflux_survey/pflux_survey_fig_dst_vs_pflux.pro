;+
; Dst, pflux avg per orbit (-A and -B).
;-

dot0 = 0

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project['pflux_calc_setting']
    probes = ['a','b']
    plot_dir = join_path([project.plot_dir,'pflux_survey','figures'])

;---Read data for each probe.
    dst_time_range = list()
    foreach probe, probes do begin
        rbspx = 'rbsp'+probe
        prefix = rbspx+'_'

        mission_time_range = pflux_calc_setting[rbspx].time_range
        dst_time_range.add, mission_time_range
        
        orbit_var = prefix+'orbit_time_range'
        if tnames(orbit_var) eq '' then begin
            orbit_time_ranges = pflux_grant_calculate_orbit_time_range(mission_time_range, probe=probe)
            store_data, orbit_var, 0, orbit_time_ranges
        endif
        orbit_time_ranges = get_var_data(orbit_var)
        norbit = n_elements(orbit_time_ranges)/2
        orbit_times = (orbit_time_ranges[*,0]+orbit_time_ranges[*,1])*0.5

        pflux_avg_var = prefix+'pf_avg'
        if keyword_set(dot0) then pflux_avg_var = prefix+'pfdot0_avg'
        if ~check_if_update(pflux_avg_var, minmax(orbit_times)) then continue

        pflux_var = prefix+'pf_fac_norm'
        if keyword_set(dot0) then pflux_var = prefix+'pfdot0_fac_norm'
        pflux_avg = fltarr(norbit)
        for ii=0, norbit-1 do begin
            time_range = reform(orbit_time_ranges[ii,*])
            pflux_grant_read_preprocessed_pflux, time_range, probe=probe
            pflux = get_var_data(pflux_var)
            pflux_avg[ii] = mean(snorm(pflux),/nan)
        endfor
        store_data, pflux_avg_var, orbit_times, pflux_avg
    endforeach

;---Read Dst.
    dst_time_range = minmax(dst_time_range.toarray())
    dst_var = 'dst'
    if check_if_update(dst_var) then omni_read_index, dst_time_range

    vars = 'rbsp'+probes+'_pf_avg'
    if keyword_set(dot0) then vars = 'rbsp'+probes+'_pfdot0_avg'
    options, vars, 'ytitle', '(mW/m!U2!N)'
    foreach var, vars, ii do begin
        options, var, 'labels', 'RB'+strupcase(probes[ii])
    endforeach

    time_range_list = list()
    time_range_list.add, time_double(['2012-12-01','2013-10-01'])
    time_range_list.add, time_double(['2014-09-01','2015-07-01'])
    time_range_list.add, time_double(['2016-06-01','2017-04-01'])

    fig_xsize = 12
    fig_ysize = 3
    margins = [10,5,5,3]
    tplot_options, 'labflag', -1
    foreach time_range, time_range_list do begin
        vars = ['dst','rbsp'+probes+'_pf_avg']
        if keyword_set(dot0) then vars = ['dst','rbsp'+probes+'_pfdot0_avg']
        options, vars, 'xticklen', -0.02
        options, vars, 'yticklen', -0.005

        ofn = join_path([plot_dir,'fig_dst_vs_pflux_avg_'+strjoin(time_string(time_range,tformat='YYYY_MMDD'),'_')+'.pdf'])
        if keyword_set(dot0) then ofn = join_path([homedir(),'fig_dst_vs_pflux_dot0_avg_'+strjoin(time_string(time_range,tformat='YYYY_MMDD'),'_')+'.pdf'])
        if keyword_set(test) then ofn = test
        nvar = n_elements(vars)
        sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize
        tpos = sgcalcpos(nvar, margins=margins, xchsz=xchsz, ychsz=ychsz)
        tplot, vars, trange=time_range, position=tpos

        tx = tpos[0]
        ty = tpos[3]+ychsz*0.3
        xyouts, tx,ty,/normal, 'Dst and |S| (Poynting flux, normalized to 100 km, and averaged value per orbit)'
        labels = letters(nvar)+'.'
        for ii=0,nvar-1 do begin
            tx = xchsz*1
            ty = tpos[3,ii]-ychsz*0.8
            xyouts, tx,ty,/normal, labels[ii]
        endfor
        if keyword_set(test) then stop
        sgclose
    endforeach

end
