;+
; Inverted-V on Jupiter and Earth.
;-

function fig_2015_0416_0800_jupiter_v01, event_info=event_info

    test = 0
    version = 'v01'
    plot_dir = join_path([googledir(),'works','pflux_grant','alfven_arc','plot'])

    label_color = sgcolor('green')
    ct = 55
    label_color = sgcolor('red')
    ct = 64   

;---Load Jupiter inverted-V data.
    data_file = join_path([srootdir(),'data','JEDI_Jupiter_potential_drop.d2s'])
    if file_test(data_file) eq 0 then message, 'Data file does not exist ...'
    lines = read_all_lines(data_file)

    infos = strsplit(lines[12],'=", ',extract=1)
    energy_bins = float(infos[1:*])
    nenergy_bin = n_elements(energy_bins)
    energy_unit = 'keV'
    
    flux_unit = '#/s-sr-cm!U2!N-keV'
    ntime = 600
    times = dblarr(ntime)
    fluxs = fltarr(ntime,nenergy_bin)
    for ii=0,ntime-1 do begin
        infos = strsplit(lines[ii+16],' ',extract=1)
        
        times[ii] = time_double(infos[0],tformat=':01:YYYY-MM-DDThh:mm:ss.fff')
        fluxs[ii,*] = float(infos[1:*])
    endfor
    
    var_juno = 'juno_eflux'
    store_data, var_juno, times, fluxs, energy_bins
    add_setting, var_juno, smart=1, dictionary($
        'display_type', 'spec', $
        'unit', flux_unit, $
        'ylog', 1, $
        'zlog', 1, $
        'zrange', [1e5,3e7], $
        'ytitle', 'Energy ('+energy_unit+')', $
        'short_name', '' )

    log_ytickv = [2,3]
    ytickn = '10!U'+string(log_ytickv,format='(I0)')+'!N'
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    yminor = 9
    yrange = [2e1,1e3]
    options, var_juno, 'ytickv', ytickv
    options, var_juno, 'yticks', yticks
    options, var_juno, 'yminor', yminor
    options, var_juno, 'ytickname', ytickn
    options, var_juno, 'yrange', yrange
    
    juno_plot_time_range = time_double(['2017-02-02/13:38:30','2017-02-02/13:39:20'])

    

    


;---Load DMSP inverted-V data.
    dmsp_plot_time_range = time_double(['2015-04-16/08:01','2015-04-16/08:06'])
    probe = 'f19'
    var_dmsp = dmsp_read_en_spec(dmsp_plot_time_range, probe=probe, species='e', errmsg=errmsg)
    get_data, var_dmsp, times, fluxs, energys
    index = where(fluxs eq 0 or finite(fluxs,nan=1), count)
    if count ne 0 then begin
        fluxs[index] = 1e-10
        store_data, var_dmsp, times, fluxs, energys
    endif    
    
    yrange = [3e2,3e4]
    log_ytickv = [3,4]
    ytickn = '10!U'+string(log_ytickv,format='(I0)')+'!N'
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    yminor = 9
    options, var_dmsp, 'yrange', yrange
    options, var_dmsp, 'ytickv', ytickv
    options, var_dmsp, 'yticks', yticks
    options, var_dmsp, 'yminor', yminor
    options, var_dmsp, 'ytickname', ytickn
    

    xpans = [1,1]
    xpad = 15
    margins = [8,2.5,8,1]
    plot_file = join_path([plot_dir,'fig_2015_0416_0800_jupiter_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    poss = panel_pos(plot_file, $
        xpans=xpans, xpad=xpad, pansize=[2.5,1], fig_size=fig_size, margins=margins)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    plot_vars = [var_dmsp,var_juno]
    nvar = n_elements(plot_vars)

    uniform_ticklen = -ychsz*0.15*fig_size[1]
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, plot_vars[pid], 'xticklen', xticklen
        options, plot_vars[pid], 'yticklen', yticklen
        var = plot_vars[pid]
        display_type = get_setting(var, 'display_type')
        if display_type eq 'spec' then begin
            zticklen = uniform_ticklen/xchsz*1/fig_size[0]
            options, var, 'zticklen', zticklen
        endif
        options, var, 'color_table', ct
        options, var, 'xstyle', 5
        
        tr = (pid eq 1)? juno_plot_time_range: dmsp_plot_time_range
        tplot, var, position=tpos, noerase=1, trange=tr, single_line_uttick=1        
    endfor


    fig_labels = letters(nvar)+') '+['Earth','Jupiter']
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*6
        if pid eq 1 then tx = tpos[0]-xchsz*9
        ty = tpos[3]-ychsz*0.6
        msg = fig_labels[pid]
        xyouts, tx,ty, msg, normal=1
    endfor
    
    ; Add labels.
    label_size = 0.8
    hsize = (keyword_set(test))? 6: 120
    thick = (keyword_set(test))? 1: 2

    dmsp_invertedv_times = time_double('2015-04-16/'+['08:02:45','08:03:22','08:03:37'])
    juno_invertedv_times = time_double('2017-02-02/'+['13:38:51','13:38:56'])
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        var = plot_vars[pid]
        is_juno = pid
        
        
        xrange = (is_juno eq 1)? juno_plot_time_range: dmsp_plot_time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            nodata=1, noerase=1, position=tpos
            
        txs = (is_juno eq 1)? juno_invertedv_times: dmsp_invertedv_times
        tx1 = mean(minmax(txs))
        tmp = convert_coord(tx1,yrange[0], data=1, to_normal=1)
        tx0 = tmp[0]
        ty0 = (is_juno eq 1)? tpos[3]-ychsz*1.5: tpos[1]+ychsz*1.5
        ty1 = (is_juno eq 1)? ty0-ychsz*2: ty0+ychsz*2
        foreach tx,txs do begin
            tmp = convert_coord(tx,yrange[0], data=1, to_normal=1)
            tx1 = tmp[0]
            arrow, tx0,ty0, tx1,ty1, normal=1, solid=1, color=label_color, hsize=hsize, thick=thick
        endforeach
        ty2 = (is_juno eq 1)? ty0+ychsz*0.2: ty0-ychsz*1
        msg = 'Inverted-V electron'
        xyouts, tx0,ty2, msg, normal=1, alignment=0.5, color=label_color;, charsize=label_size
    endfor
    
    
;---Manual x-axis.
    pid = where(plot_vars eq var_juno, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        var = plot_vars[pid]
        is_juno = pid

        xrange = (is_juno eq 1)? juno_plot_time_range: dmsp_plot_time_range
        yrange = get_setting(var, 'yrange')
        xticklen = get_setting(var, 'xticklen')
        
        xstep = 20
        xtickv = smkarthm(xrange[0],xrange[1],xstep,'dx')
        xticks = n_elements(xtickv)-1
        xminor = xstep/10
        xtickn = time_string(xtickv,tformat='hh:mm:ss')
        xtickn[0] = time_string(xtickv[0],tformat='MTH DD')
        
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            nodata=1, noerase=1, position=tpos
        axis, xaxis=0, save=1, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xticklen=xticklen
        axis, xaxis=1, save=1, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, $
            xtickname=xtickn, xticklen=xticklen, xtickformat='(A1)'
    endif
    
    pid = where(plot_vars eq var_dmsp, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        var = plot_vars[pid]
        is_juno = pid

        xrange = (is_juno eq 1)? juno_plot_time_range: dmsp_plot_time_range
        yrange = get_setting(var, 'yrange')
        xticklen = get_setting(var, 'xticklen')

        xstep = 120
        xtickv = smkarthm(xrange[0],xrange[1],xstep,'dx')
        xticks = n_elements(xtickv)-1
        xminor = xstep/60
        xtickn = time_string(xtickv,tformat='hh:mm')
        xtickn[0] = time_string(xtickv[0],tformat='MTH DD')

        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            nodata=1, noerase=1, position=tpos
        axis, xaxis=0, save=1, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xticklen=xticklen
        axis, xaxis=1, save=1, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, $
            xtickname=xtickn, xticklen=xticklen, xtickformat='(A1)'
    endif

    if keyword_set(test) then stop
    sgclose

    return, plot_file



end

print, fig_2015_0416_0800_jupiter_v01(event_info=event_info)
end
