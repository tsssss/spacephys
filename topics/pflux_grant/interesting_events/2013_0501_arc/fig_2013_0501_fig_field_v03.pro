;+
; Show wave and particle data.
;-

function fig_2013_0501_fig_field_v03, plot_file, event_info=event_info

test = 0



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

    label_size = 0.8
    
    
    ; n var.
    n_var = prefix+'density_plot'
    get_data, prefix+'density_hope', times, data, limits=lim
    index = where(finite(data))
    data = interpol(data[index], times[index], times)
    store_data, n_var, times, data, limits=lim
    options, n_var, 'labels', 'Density!C  >200 eV'
    options, n_var, 'ytitle', '(cm!U-3 !N)'
    
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
    yminor = 4
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
    db_var = duplicate_var(prefix+'b1_fac', output=prefix+'b1_fac_plot')
    var = db_var
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
    options, var, 'labels', 'FAC dB!D'+fac_labels+'!N'

    ; B var.
    b_var = prefix+'b_sm'
    var = prefix+'b_gsm'
    get_data, var, times, b_gsm, limits=lim
    b_sm = cotran(b_gsm, times, 'gsm2sm')
    store_data, b_var, times, b_sm
    add_setting, b_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'B', $
        'unit', 'nT', $
        'coord', 'SM', $
        'coord_labels', constant('xyz') )
    yrange = [-90,190]
    ystep = 100
    yminor = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    constant = ytickv[where_pro(ytickv,'()', yrange)]
    var = b_var
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', constant


    
    ; electron energy spectrogram.
    e_spec_var = duplicate_var(prefix+'e_en_spec', output=prefix+'e_en_spec_plot')
    var = e_spec_var
    ct_electron = 65
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
    options, var, 'ytitle', 'Energy!C(eV)'
    zrange = [1e5,1e10]
    log_zrange = alog10(zrange)
    log_ztickv = make_bins(log_zrange,2,inner=1)
    ztickv = 10.^log_ztickv
    zticks = n_elements(ztickv)-1
    options, var, 'zrange', zrange
    options, var, 'ztickv', ztickv
    options, var, 'zticks', zticks
    options, var, 'zcharsize', label_size
    options, var, 'color_table', ct_electron




;---Plot.
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,'fig_2013_0501_fig_field_v03.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz

    plot_vars = [n_var,e_spec_var,e_vars[1:2],db_var,b_var,vexb_vars[1:2]]
    nplot_var = n_elements(plot_vars)
    nvar = nplot_var
    fig_letters = letters(nvar)
    fig_labels = fig_letters+') '+['N!De!N','e!E-!N','E!I'+fac_labels[1:2],$
        'dB','SM B','V!I'+fac_labels[1:2]]
    fig_labels = fig_letters+') '
        

    margins = [10,4,10,1]


    poss_left = panel_pos(plot_file, pansize=[4,0.8], $
        ypans=1+fltarr(nplot_var), margins=margins, $
        fig_size=fig_size)


    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    uniform_ticklen = -ychsz*0.2*fig_size[0]
    

;---tplot panels.
    plot_poss = poss_left

    for ii=0,nplot_var-1 do begin
        tpos = plot_poss[*,ii]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, plot_vars[ii], 'xticklen', xticklen
        options, plot_vars[ii], 'yticklen', yticklen
    endfor

    tplot, plot_vars, trange=time_range, position=plot_poss, noerase=1

    ; Add fig label.
    for ii=0,nplot_var-1 do begin
        tpos = plot_poss[*,ii]
        fig_label = fig_labels[ii]
        tx = 2*xchsz
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,fig_label, normal=1
    endfor
    timebar, bar_times, linestyle=1, color=sgcolor('black')
    timebar, snapshot_time, color=sgcolor('red')


    ; Add electron for density and en_spec.
    foreach var, [n_var,e_spec_var] do begin
        pid = where(plot_vars eq var, count)
        if count eq 0 then continue
        tpos = plot_poss[*,pid]
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = 'Electron'
        xyouts, tx,ty,normal=1, msg, alignment=1, charsize=label_size
    endforeach

    ; Add electron temperature.
    var = e_spec_var
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = plot_poss[*,pid]
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*1

        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        
        get_data, prefix+'e_temp', times, data
        oplot, times, data
    endif


    if keyword_set(test) then stop
    sgclose

    return, plot_file
    

end



print, fig_2013_0501_fig_field_v03(event_info=event_info)
end