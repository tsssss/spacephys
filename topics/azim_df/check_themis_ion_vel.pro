
test = 0

test_time = time_double('2014-08-28/10:20')
test_time = !null
vel_type = 'u'

    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_find_dfgroup(project=project)

    foreach event, events do begin
        direction = (strsplit(event.region,'%',/extract))[1]
        lprmsg, direction
        lprmsg, time_string(event.time_range)

        if direction ne 'eastward' and direction ne 'westward' then continue
        if keyword_set(test_time) then begin
            if product(event.time_range-test_time[0]) gt 0 then continue
        endif
        
        
        probes = event.probes
        df_list = list()
        foreach df, event.df_list do begin
            probe = df.probe
            if strmid(probe,0,2) eq 'th' then df_list.add, df
        endforeach
        if df_list.length eq 0 then continue
        
        
        long_time_range = event.time_range+[-1,1]*600
        ndf = df_list.length
        probes = strarr(ndf)
        azim_vel = fltarr(ndf)+!values.f_nan
        time_step = 3
        foreach df, df_list, df_id do begin
            probe = df.probe
            probes[df_id] = probe
            time_range = df.time_range
            the_probe = strmid(probe,2,1)
            check_themis_vel, long_time_range, probe=the_probe
            prefix = probe+'_'
            v_fac_var = prefix+vel_type+'_fac'
            stplot_split, v_fac_var, newnames=prefix+vel_type+['r','w','z']

            val = get_var_data(prefix+vel_type+'w', in=time_range)
            nval = total(time_range*[-1,1])/time_step
            if n_elements(val) lt nval*0.9 then continue
            azim_vel[df_id] = mean(val,/nan)
            
            azim_dp_read_theta, long_time_range, probe=probe
        endforeach
        
        index = where(finite(azim_vel), count)
        if count eq 0 then continue
        df_list = df_list[index]
        probes = probes[index]
        azim_vel = azim_vel[index]
        
        lprmsg, time_string(event.time_range)
        print, direction
        print, probes
        print, azim_vel
        
        vars = []
        foreach probe, probes do begin
            vars = [vars, probe+'_'+['theta',vel_type+'w']]
            options, probe+'_theta', 'labels', strupcase(probe)+' '+tex2str('theta')
        endforeach
        
        options, probes+'_'+vel_type+'w', 'labels', 'V_west'
        tplot_options, 'labflag', -1
        tplot_options, 'constant', 0
        
        plot_file = join_path([homedir(),'check_themis_ion_vel','fig_event_'+$
            strjoin(time_string(event.time_range,tformat='YYYY_MMDD_hh'),'_')+'.pdf'])
        if keyword_set(test) then plot_file = 0
        
        nvar = n_elements(vars)
        fig = fig_init(plot_file)
        margins = [8d,4,6,1]
        fig_replace_region, fig, region_init(fig, $
            nypan=nvar, xsize=4, pansize=[4,1], margins=margins)
        size = fig_pos_calc_size(fig)
        
        region = fig.regions[0]
        sgopen, fig.id, xsize=region.xsize, ysize=region.ysize, $
            xchsz=xchsz, ychsz=ychsz
        poss = region.abs_pos
        poss[[0,2],*,*] /= region.xsize
        poss[[1,3],*,*] /= region.ysize
        poss = reform(poss)
        tplot, vars, trange=long_time_range, position=poss, /novtitle
        
        fig_labels = letters(nvar)+'.'
        for panel_id=0, nvar-1 do begin
            tpos = poss[*,panel_id]
            tx = xchsz*2
            ty = tpos[3]-ychsz*0.8
            xyouts, tx,ty,/normal, fig_labels[panel_id]
        endfor
        
        foreach df, df_list, df_id do begin
            tpos = poss[*,df_id*2+1]
            yrange = [0,1]
            plot, long_time_range, yrange, xstyle=5, ystyle=5, $
                position=tpos, nodata=1, noerase=1
            foreach tx, df.time_range do plots, tx+[0,0], yrange
            tx = tpos[0]+0.5*xchsz
            ty = tpos[1]+0.3*ychsz
            xyouts, tx,ty,/normal, string(azim_vel[df_id],format='(F7.1)')+' km/s'
        endforeach
        
        tpos = poss[*,0]
        tx = tpos[0]+0.5*xchsz
        ty = tpos[3]-ychsz*0.9
        xyouts, tx,ty,/normal, direction
        
        if keyword_set(test) then stop
        sgclose
    endforeach

end
