;+
; An overview figure, showing the conjunction of wave and arc, and ion beams.
;-

function fig_2013_0501_arc_overview, plot_file, event_info=event_info

test = 1



;---Load data.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()


;---Settings.
    prefix = event_info['prefix']
    time_range = event_info['time_range']

    e_var = prefix+'edot0_fac'
    var = e_var
    yrange = [-1,1]*190
    ystep = 100
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 5
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', [0,-100,100]
    get_data, var, limits=lim
    labels = lim.labels
    ndim = n_elements(labels)
    for ii=0,ndim-1 do begin
        label = labels[ii]
        index = strpos(label,'dot0')
        if index[0] eq -1 then continue
        labels[ii] = strmid(label,0,index)+strmid(label,index+4)
    endfor
    options, var, 'labels', labels

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

    pf_var = prefix+'pfdot0_fac_map'
    var = pf_var
    yrange = [-100,500]
    ystep = 200
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 2
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', [0,200,400]

    pf_spec_var = prefix+'pfdot0_fac_mor_spec_1'
    var = pf_spec_var
    options, var, 'constant', [10,100,1000]
    
    
    
    
;---Plot.
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, xsize=6, ysize=6
    
    label_size = 0.8
    
    margins = [10,4,12,1]
    plot_vars = [e_var,b_var,pf_var,pf_spec_var]
    nvar = n_elements(plot_vars)
    poss = sgcalcpos(nvar, margins=margins, xchsz=xchsz, ychsz=ychsz)
    tplot, plot_vars, trange=time_range, position=poss
    
    pid = where(plot_vars eq pf_var)
    tpos = poss[*,pid]
    xrange = time_range
    yrange = get_setting(pf_var, 'yrange')
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = 'Norm to 110 km altitude'
    xyouts, tx,ty,normal=1, msg, charsize=label_size
    
    tmp = convert_coord(xrange[0],0, data=1, to_normal=1)
    tx = tmp[0]+xchsz*0.5
    ty = tmp[1]+ychsz*0.4
    msg = 'To N-hem'
    xyouts, tx,ty,normal=1, msg, charsize=label_size
    tx = tmp[0]+xchsz*0.5
    ty = tmp[1]-ychsz*0.8
    msg = 'To S-hem'
    xyouts, tx,ty,normal=1, msg, charsize=label_size
    
    bar_times = make_bins(time_range,600, inner=1)
    timebar, bar_times, linestyle=1, color=sgcolor('black')
    stop

end

print, fig_2013_0501_arc_overview(event_info=event_info)
end