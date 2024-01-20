;+
; Show SAPS E field and the association of auroral boundary.
;-


function fig_saps_v01, plot_dir, event_info=event_info, test=test

    version = 'v01'
    id = '2015_0317'
    event_info = low_lshell_outflow_load_data(id)

    time_range = event_info.time_range

    if keyword_set(test) then plot_file = 0
    fig_size = [6,6]
    sgopen, plot_file, size=fig_size


    saps_times = dictionary()
    saps_times['rbspa'] = ['2015-03-17/15:30','2015-03-17/16:00']
    saps_times['rbspb'] = ['2015-03-17/21:05','2015-03-17/21:35']
    colors = sgcolor(['purple','red'])

    margins = [5,4,2,1]
    tpos = sgcalcpos(1, margins=margins)
    xrange = [2,3.5]
    xtitle = 'L-shell (#)'
    xrange = [45d,60]
    xtitle = 'ILat'
    yrange = [-5,15]
    ytitle = 'E!D'+tex2str('perp')+',out!N (mV/m)'
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, $
        nodata=1, noerase=1, position=tpos
    foreach mission_probe, saps_times.keys(), probe_id do begin
        probe_info = resolve_probe(mission_probe)
        prefix = probe_info['prefix']
        probe = probe_info['probe']
        lshell_var = prefix+'lshell'
        e_var = prefix+'edot0_fac'
        tr = time_double(saps_times[mission_probe])
        lshell = get_var_data(lshell_var, in=tr, times=times)
        ilat = lets_calc_ilat(lshell)
        mlt_var = prefix+'mlt'
        mlt = get_var_data(mlt_var, at=times)
        efac = get_var_data(e_var, at=times)

        plots, ilat, efac[*,2], color=colors[probe_id]
        print, minmax(mlt), mean(mlt)
    endforeach
    
    plot, xrange, yrange, $
        xstyle=1, ystyle=1, $
        xtitle=xtitle, ytitle=ytitle, $
        nodata=1, noerase=1, position=tpos

    

    if keyword_set(test) then stop
    sgclose
    return, plot_file

end

print, fig_saps_v01(event_info=event_info, test=1)
end