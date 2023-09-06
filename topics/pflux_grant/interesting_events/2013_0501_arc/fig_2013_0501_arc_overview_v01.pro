;+
; An overview figure, showing the conjunction of wave and arc, and ion beams.
;-

function fig_2013_0501_arc_overview_v01, plot_file, event_info=event_info

test = 1



;---Load data.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()


;---Settings.
    prefix = event_info['prefix']
    probe = event_info['probe']
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
    foreach label, labels, ii do begin
        index = strpos(label,'dot0')
        if index[0] eq -1 then continue
        labels[ii] = strmid(label,0,index)+strmid(label,index+4)
    endforeach
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
    
    
    vars = [e_var,b_var,pf_var]
    foreach var, vars do begin
        labels = get_setting(var,'labels')
        foreach label, labels, ii do begin
            index = strpos(label,'FAC')
            if index[0] eq -1 then continue
            labels[ii] = strmid(label,index[0]+4)
        endforeach
        options, var, 'labels', labels
    endforeach


;---Fpt ASI count.
    model_setting = event_info['model_setting']
    external_model = model_setting['external_model']
    internal_model = model_setting['internal_model']
    fmlt_var = prefix+'fmlt_'+internal_model
    fmlat_var = prefix+'fmlat_'+internal_model
    models = model_setting['models']
    model_index = where(models eq external_model)

    mlt_image_var = 'thg_asf_mlt_image_rect'
    get_data, mlt_image_var, times, mlt_images, limits=lim
    mlt_bins = lim.mlt_bins
    mlat_bins = lim.mlat_bins
    ntime = n_elements(times)
    asi_count = fltarr(ntime)

    fmlt = (get_var_data(fmlt_var, at=times))[*,model_index]
    fmlat = (get_var_data(fmlat_var, at=times))[*,model_index]
    dmlat = 1d      ; deg.
    dmlt = dmlat/15 ; h.

    foreach time, times, time_id do begin
        mlt_index = where_pro(mlt_bins, '[]', fmlt[time_id]+[-1,1]*dmlt*0.5, count=count)
        if count eq 0 then continue
        mlat_index = where_pro(mlat_bins, '[]', fmlat[time_id]+[-1,1]*dmlat*0.5, count=count)
        if count eq 0 then continue

        asi_count[time_id] = mean(mlt_images[time_id,mlt_index,mlat_index])
    endforeach
    asi_count_var = 'asi_count'
    store_data, asi_count_var, times, asi_count
    add_setting, asi_count_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'ASI count', $
        'unit', '#')
    yrange = [0,2.5e4]
    ytickv = [1,2]*1e4
    yticks = n_elements(ytickv)-1
    yminor = 5
    var = asi_count_var
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', ytickv
    
;---Plot.
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, xsize=4.5, ysize=6
    
    label_size = 0.8
    
    margins = [10,4,12,1]
    plot_vars = [e_var,b_var,pf_var,asi_count_var,pf_spec_var]
    nvar = n_elements(plot_vars)
    poss = sgcalcpos(nvar, margins=margins, xchsz=xchsz, ychsz=ychsz)
    tplot, plot_vars, trange=time_range, position=poss
    
    ; Add labels for pflux.
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
    
    
    ; Add label for ASI count.
    pid = where(plot_vars eq asi_count_var)
    tpos = poss[*,pid]
    xrange = time_range
    yrange = get_setting(asi_count_var, 'yrange')
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    tmp = string(dmlat,format='(I0)')
;    msg = tmp+tex2str('times')+tmp+' deg around!CRBSP-'+strupcase(probe)+' F/point ('+string(external_model)+')'
    msg = tmp+tex2str('times')+tmp+' deg around!Cfootpoint ('+string(external_model)+')'
    xyouts, tx,ty,normal=1, msg, charsize=label_size
    
    ; Add vertical grids and times.
    bar_times = make_bins(time_range,600, inner=1)
    timebar, bar_times, linestyle=1, color=sgcolor('black')
    
    snapshot_time = time_double('2013-05-01/07:38:09')
    timebar, snapshot_time, color=sgcolor('red'), linestyle=3

end

print, fig_2013_0501_arc_overview_v01(event_info=event_info)
end