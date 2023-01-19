;+
; An overview figure, showing the conjunction of wave and arc, and ion beams.
;-

function fig_2013_0501_field_v01, plot_file, event_info=event_info

test = 1



;---Load data.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
    plasma_param = event_info['plasma_param']

;---Settings.
    prefix = event_info['prefix']
    probe = event_info['probe']
    time_range = event_info['time_range']
    snapshot_time = event_info['snapshot_time']
    fac_labels = event_info['fac_labels']
    bar_times = make_bins(time_range,600, inner=1)
    
    
    xticklen_chsz = -0.2
    yticklen_chsz = -0.35
    
    ; n var.
    n_var = prefix+'plot_density'
    get_data, prefix+'density_hope', times, data, limits=lim
    index = where(finite(data))
    data = interpol(data[index], times[index], times)
    store_data, n_var, times, data, limits=lim
    
    var = n_var
    yrange = [0.1,3]
    ytickv = [0.1,1]
    yticks = n_elements(ytickv)-1
    yminor = 10
    constant = [1]
    options, var, 'ystyle', 1
    options, var, 'ylog', 1
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant


    ; vexb var.
    vexb_var = prefix+'vexbdot0_fac'
    vexb_vars = prefix+'vexbdot0_'+['para','west','out']
    vexb_labels = 'ExB V!D'+fac_labels+'!N'
    stplot_split, vexb_var, newnames=vexb_vars, labels=vexb_labels
    
    var = vexb_vars[1]
    yrange = [-1,1]*600
    ystep = 400
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 4
    constant = [0,-1,1]*ystep
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant
    
    var = vexb_vars[2]
    yrange = [-1,1]*300
    ystep = 200
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 2
    constant = [0,-1,1]*ystep
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant
    
    foreach var, vexb_vars do begin
        uniform_time, var, 5
    endforeach


    ; e vars.
    e_var = prefix+'edot0_fac'
    e_vars = prefix+'edot0_'+['para','west','out']
    e_labels = 'E!D'+fac_labels+'!N'
    stplot_split, e_var, newnames=e_vars, labels=e_labels
    
    var = e_vars
    yrange = [-1,1]*190
    ystep = 100
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 5
    constant = [0,-1,1]*ystep
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant
    
    
    ; dB vars.
    b_var = prefix+'b1_fac'
    var = b_var
    yrange = [-1,1]*60
    ystep = 50
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 5
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', 0



    ; O+ outflow.
    ;vars = rbsp_read_en_spec_combo(time_range, probe=probe, species='o')
    spec_vars = prefix+[$
        'o_en_spec_'+['anti','para'], $
        'p_en_spec_'+['anti','para'] ]
    foreach var, spec_vars do begin
        get_data, var, times, data, vals
        index = where(data eq 0 or finite(data,nan=1), count)
        if count ne 0 then begin
            data[index] = 1e-2
            store_data, var, times, data, vals
        endif
    endforeach
    yrange = [5,5e4]
    log_yrange = alog10(yrange)
    log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
    log_ytickv = make_bins(log_yrange,1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    yminor = 10
    foreach tx, ytickv, ii do begin
        if tx eq 1 then begin
            ytickn[ii] = '1'
        endif else if tx eq 10 then begin
            ytickn[ii] = '10'
        endif
    endforeach

    var = [spec_vars,prefix+['p_en_spec','e_en_spec']]
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'ytickname', ytickn
    options, var, 'constant', ytickv
    options, var, 'color_table', 64
    options, var, 'ytitle', 'Energy!C(eV)'
    
    var = [spec_vars,prefix+'p_en_spec']
    options, var, 'zrange', [1e4,1e6]
    
    
    var = prefix+'e_en_spec'
    yrange = [15,5e4]
    log_yrange = alog10(yrange)
    log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
    log_ytickv = make_bins(log_yrange,1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    yminor = 10
    foreach tx, ytickv, ii do begin
        if tx eq 1 then begin
            ytickn[ii] = '1'
        endif else if tx eq 10 then begin
            ytickn[ii] = '10'
        endif
    endforeach
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'ytickname', ytickn
    options, var, 'constant', ytickv
    options, var, 'color_table', 65
    options, var, 'ytitle', 'Energy!C(eV)'
    zrange = [1e5,1e10]
    log_zrange = alog10(zrange)
    log_ztickv = make_bins(log_zrange,2,inner=1)
    ztickv = 10.^log_ztickv
    zticks = n_elements(ztickv)-1
    options, var, 'zrange', zrange
    options, var, 'ztickv', ztickv
    options, var, 'zticks', zticks
    
    
    var = prefix+['p_en_spec','p_en_spec_'+['anti','perp','para']]
    options, var, 'color_table', 63
    
    
;---Plot.
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,'fig_2013_0501_field_v01.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz

    plot_vars = [n_var,e_vars[1:2],b_var,vexb_vars[1:2]]
    nplot_var = n_elements(plot_vars)
    nvar = nplot_var
    fig_letters = letters(nvar*2)
    fig_labels = fig_letters[0:nplot_var-1]+') '+['N','E!D'+fac_labels[1:2],$
        'dB','V!D'+fac_labels[1:2]]
        

    margins = [10,4,10,1]
    label_size = 0.8

    tmp = panel_pos(plot_file, pansize=[3,1], $
        xpans=[1,1], ypans=1+fltarr(nplot_var), margins=margins, $
        fig_size=fig_size, xpad=17)

    poss_left = reform(tmp[*,0,*])
    poss_right = reform(tmp[*,1,*])

    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    

;---tplot panels.
    plot_poss = poss_left

    for ii=0,nplot_var-1 do begin
        tpos = plot_poss[*,ii]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        options, plot_vars[ii], 'xticklen', xticklen
        options, plot_vars[ii], 'yticklen', yticklen
    endfor

    tplot, plot_vars, trange=time_range, position=plot_poss, noerase=1
    for ii=0,nplot_var-1 do begin
        tpos = plot_poss[*,ii]
        fig_label = fig_labels[ii]
        tx = 2*xchsz
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,fig_label, normal=1
    endfor
    timebar, bar_times, linestyle=1, color=sgcolor('black')
    timebar, snapshot_time, color=sgcolor('red')

    
    spec_poss = poss_right
    spec_vars = prefix+['e_en_spec','p_en_spec', $
        'p_en_spec_'+['para','anti'], $
        'o_en_spec_'+['para','anti'] ]
    spec_labels = fig_letters[nplot_var:-1]+') '+['e-',$
        'H+','H+ '+['para','anti'], $
        'O+ '+['para','anti'] ]
    
    for ii=0,nplot_var-1 do begin
        tpos = spec_poss[*,ii]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        options, spec_vars[ii], 'xticklen', xticklen
        options, spec_vars[ii], 'yticklen', yticklen
    endfor
    
    tplot, spec_vars, trange=time_range, position=spec_poss, noerase=1
    for ii=0,nplot_var-1 do begin
        tpos = spec_poss[*,ii]
        fig_label = spec_labels[ii]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,fig_label, normal=1
    endfor
    timebar, bar_times, linestyle=1, color=sgcolor('black')
    timebar, snapshot_time, color=sgcolor('red')
    stop
    
    
    if keyword_set(test) then stop
    sgclose
    
    

end

print, fig_2013_0501_field_v01(event_info=event_info)
end