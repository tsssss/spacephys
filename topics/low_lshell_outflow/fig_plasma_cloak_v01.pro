
function fig_plasma_cloak_v01, plot_dir, event_info=event_info, test=test

    version = 'v01'
    id = '2015_0317'
    if n_elements(event_info) eq 0 then event_info = low_lshell_outflow_load_data(id)
    label_size = 0.8
    
    time_range = time_double(['2015-03-17','2015-03-18'])
    rbsp_probes = ['a','b']
    plasma_cloak_times = dictionary($
        'a', time_double([$
            '2015-03-17/13:46','2015-03-17/14:06',$
            '2015-03-17/15:48','2015-03-17/16:34',$
            '2015-03-17/22:44','2015-03-17/23:18']), $
        'b', time_double([$
            '2015-03-17/09:44','2015-03-17/10:24',$
            '2015-03-17/12:34','2015-03-17/12:58',$
            '2015-03-17/19:28','2015-03-17/19:45',$
            '2015-03-17/21:16','2015-03-17/21:46']) )

    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'

        p_en_var = prefix+'p_en_spec'
        o_en_var = rbsp_read_en_spec(time_range, probe=probe, species='o')
        o_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='o', errmsg=errmsg, energy_range=[50,5e4], update=1)

        the_vars = [p_en_var,o_en_var]
        options, the_vars, 'ytitle', 'Energy!C(eV)'
        the_vars = [p_en_var,o_en_var]
        options, the_vars, ytickname=['10!U2','10!U3','10!U4'], $
            ytickv=[100,1000,10000], yticks=2, yminor=9, $
            yrange=[20,3e4]

        the_vars = o_pa_var
        options, the_vars, 'ytitle', 'PA!C(deg)'

        zticklen = -0.5
        the_vars = [p_en_var,o_en_var,o_pa_var]
        options, the_vars, zticklen=zticklen, zminor=9

        ; dis.
        lshell_var = prefix+'lshell'
        var = lshell_var
        options, var, 'yrange', [1,7]
        options, var, 'ytickv', [2,4,6]
        options, var, 'yticks', 2
        options, var, 'yminor', 2
        options, var, 'constant', [2,4,6]
        options, var, 'ytitle', '(#)'
        options, var, 'labels', 'L-shell'
    endforeach

    plot_vars = []
    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'
        plot_vars = [plot_vars, prefix+['o_en_spec','o_pa_spec','lshell']]
    endforeach
    nvar = n_elements(plot_vars)
    fig_letters = letters(nvar)
    tmp = ['O!U+!N '+['EN','PA'],'L-shell']
    tmp = [['EN','PA'],'L-shell']
    fig_labels = []
    foreach probe, rbsp_probes do begin
        fig_labels = [fig_labels,probe+'-'+string(findgen(n_elements(tmp))+1,format='(I0)')+') '+tmp]
    endforeach
    tmp = [1,1,0.6]
    ypans = [tmp,tmp]


    pansize = [6,0.7]
    margins = [10,3,7,1]
    
    if n_elements(plot_dir) eq 0 then plot_dir = event_info.plot_dir
    plot_file = join_path([plot_dir,'fig_plasma_cloak_'+id+'_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    poss = panel_pos(plot_file, ypans=ypans, fig_size=fig_size, $
        nypan=nvar, pansize=pansize, panid=[0,1], margins=margins)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz, inch=1

    xticklen_chsz = -0.25   ; in ychsz.
    yticklen_chsz = -0.40   ; in xchsz.
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        options, plot_vars[ii], 'xticklen', xticklen
        options, plot_vars[ii], 'yticklen', yticklen
    endfor

    tplot_options, 'tickinterval', 3600*4
    tplot_options, 'version', 3
    var_labels = ''
    tplot, plot_vars, trange=time_range, position=poss
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*(margins[0]-1)
        ty = tpos[3]-ychsz*0.7
        msg = fig_labels[pid]
        xyouts, tx,ty,normal=1, msg
    endfor
    
    ; Add labels
    foreach probe, rbsp_probes, probe_id do begin
        msg = strupcase('rbsp-'+probe)+' O+'
        pid = probe_id*3
        tpos = poss[*,pid]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,normal=1, msg, color=sgcolor('white')
    endforeach

    ; Add pp_times
    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'

        pp_times = ((event_info['rbsp'])['rbsp'+probe])['pp_times']
        spec_var = prefix+'o_en_spec'
        pid = where(plot_vars eq spec_var)
        tpos = poss[*,pid]
        tpos[1] = poss[1,pid+2]

        xr = time_range
        yr = [0,1]
        plot, xr,yr, $
            xstyle=5, xrange=xr, ystyle=5, yrange=yr, $
            position=tpos, nodata=1, noerase=1
        foreach pp_time, pp_times do begin
            oplot, pp_time+[0,0], yr, linestyle=2, color=sgcolor('red')
        endforeach
    endforeach


    plasma_cloak_color = sgcolor('red')
    void_color = sgcolor('light_blue')
    thick = (keyword_set(test))? 2: 8
    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'
        foreach pa_var, prefix+['o_pa_spec','o_en_spec'] do begin
            pid = where(plot_vars eq pa_var, count)
            if count eq 0 then continue
            tpos = poss[*,pid]
            xrange = time_range
            yrange = [0,1]
            plot, xrange, yrange, $
                xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                position=tpos, nodata=1, noerase=1
            the_times = plasma_cloak_times[probe]
            nsector = n_elements(the_times)*0.5

            tys = mean(yrange)+[0,0]
            for ii=0,nsector-1 do begin
                txs = the_times[ii*2:ii*2+1]
                plots, txs, tys, data=1, color=plasma_cloak_color
                
                if pa_var eq prefix+'o_pa_spec' then begin
                    tx = mean(txs)
                    ty = tys[0]
                    tmp = convert_coord(tx,ty, data=1, to_normal=1)
                    tx = tmp[0]
                    ty = tmp[1]+ychsz*0.2
                    msg = 'PC'+string(ii+1,format='(I0)')
                    xyouts,tx,ty,normal=1, msg, alignment=0.5, color=plasma_cloak_color, charsize=label_size
                endif
                
                foreach tx, txs do begin
                    tmp = convert_coord(tx,tys[0],data=1,to_normal=1)
                    tx = tmp[0]
                    ty = tmp[1]
                    plots, tx+[0,0], ty+[-1,1]*ychsz*0.15, normal=1, color=plasma_cloak_color
                endforeach
            endfor
        endforeach
        
        lshell_var = prefix+'lshell'
        pid = where(plot_vars eq lshell_var, count)
        if count eq 0 then continue
        tpos = poss[*,pid]
        xrange = time_range
        yrange = get_var_setting(lshell_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        the_times = plasma_cloak_times[probe]
        nsector = n_elements(the_times)*0.5
        for ii=0,nsector-1 do begin
            tr = the_times[ii*2:ii*2+1]
            tys = get_var_data(lshell_var, in=tr, times=txs)
            plots, txs,tys, color=plasma_cloak_color, thick=thick
        endfor
        
        
        ; Get the times for L in [3-4].
        lshell_range = [4,5]
        lshells = get_var_data(lshell_var, times=times)
        index = where_pro(lshells, '[]', lshell_range)
        trs = times[time_to_range(index,time_step=1)]
        for ii=1,1 do begin
            tr = reform(trs[ii,*])
            
            pid = where(plot_vars eq lshell_var, count)
            if count eq 0 then continue
            tpos = poss[*,pid]
            xrange = time_range
            yrange = get_var_setting(lshell_var, 'yrange')
            plot, xrange, yrange, $
                xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                position=tpos, nodata=1, noerase=1
            tys = get_var_data(lshell_var, in=tr, times=txs)
            plots, txs,tys, color=void_color, thick=thick
            
            foreach pa_var, prefix+['o_pa_spec','o_en_spec'] do begin
                pid = where(plot_vars eq pa_var, count)
                if count eq 0 then continue
                tpos = poss[*,pid]
                xrange = time_range
                yrange = [0,1]
                plot, xrange, yrange, $
                    xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
                    position=tpos, nodata=1, noerase=1
                

                tys = mean(yrange)+[0,0]
                txs = tr
                plots, txs, tys, data=1, color=void_color
                
                if pa_var eq prefix+'o_pa_spec' then begin
                    tx = mean(txs)
                    ty = tys[0]
                    tmp = convert_coord(tx,ty, data=1, to_normal=1)
                    tx = tmp[0]
                    ty = tmp[1]+ychsz*0.2
                    msg = 'NoPC';+string(ii,format='(I0)')
                    xyouts,tx,ty,normal=1, msg, alignment=0.5, color=void_color, charsize=label_size
                endif
                
                
                foreach tx, txs do begin
                    tmp = convert_coord(tx,tys[0],data=1,to_normal=1)
                    tx = tmp[0]
                    ty = tmp[1]
                    plots, tx+[0,0], ty+[-1,1]*ychsz*0.15, normal=1, color=void_color
                endforeach
            endforeach
        endfor
    endforeach

    if keyword_set(test) then stop
    sgclose

    return, plot_file

end

print, fig_plasma_cloak_v01(event_info=event_info, test=0)
end