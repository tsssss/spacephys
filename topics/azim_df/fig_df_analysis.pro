;+
; DF processing steps.
;-

    test = 0
    time_range = time_double(['2014-08-28/09:30','2014-08-28/11:30'])

;---Check inputs.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(time_range) ne 2 then begin
        errmsg = handle_error('No input time_range ...')
        return
    endif


;---Find the event.
    ;events = azim_df_search_all_events(project=project)
    events = azim_df_find_dfgroup(project=project)
    foreach event, events do if product(event.time_range-time_range) lt 0 then break
    azim_df_load_basic_data, project=project


    probe = 'thd'
    prefix = probe+'_'
    df_list = event.df_list
    foreach df, df_list do if df.probe eq probe then break


;---Get the size of the figure, and position of the panels.
    plot_file = join_path([project.plot_dir,'fig_df_analysis.pdf'])
    if keyword_set(test) then plot_file = test
    magnify = keyword_set(test)? 2: 1
    sgopen, plot_file, xsize=4, ysize=3, /inch, xchsz=xchsz, ychsz=ychsz, magnify=magnify

    xpad = 2
    ypad = 0.4
    margins = [8,3,2,2]
    nypanel = 3
    poss = sgcalcpos(nypanel, ypad=ypad, margins=margins)
    labels = letters(nypanel)+'.'
    for ii=0,nypanel-1 do begin
        tpos = poss[*,ii]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,labels[ii], /normal
    endfor
    label_size = constant('label_size')
    half_ychsz = constant('half_ychsz')
    full_ychsz = constant('full_ychsz')
    tmp = findgen(21)/20*2*!dpi
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs,tys,/fill

;---Common x-axis.
    xrange = time_range
    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    xminor = 10 ; min.
    xstep = 30*60
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    secofday = 86400d
    for ii=0, xticks do begin
        the_time = xtickv[ii]
        xtickn[ii] = time_string(the_time,tformat='hh:mm')
        date = time_string(the_time,tformat='YYYY-MM-DD')+'        '
        if ii eq 0 then begin
            xtickn[ii] = date
            continue
        endif
    endfor


;---Settings.
    boxcar_window = 120.    ; sec.
    boxcar_ratio = 0.8
    df_min_duration = 60.   ; sec.
    dtheta_nsigma = 1.      ; to select significant peaks of positive slope.


    theta_var = prefix+'theta'
    azim_df_smooth_theta, theta_var, time_range, smooth_window=boxcar_window, stat_ratio=boxcar_ratio
    theta = get_var_data(theta_var, in=time_range, times=times)
    theta_smooth_combo_var = theta_var+'_smooth_combo'
    tmp = get_var_data(theta_smooth_combo_var)
    theta_smooth = tmp[*,0]
    theta_stddev = tmp[*,1]
    dtheta = tmp[*,2]
    dtheta_stddev = tmp[*,3]

    azim_df_read_data, 'b_tilt', probe=probe, time_range=time_range, project=project
    theta = get_var_data(prefix+'theta', times=times, in=time_range)
    alpha = get_var_data(prefix+'b_tilt', at=times)
    theta_bg = alpha-theta


;---Panel 1. The original theta.
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    second_color = sgcolor('blue')

    xtickformat = '(A1)'
    yrange = minmax(theta)
    yrange = [-1,1]*max(abs(make_bins(yrange,2)))
    yticks = 2
    ytitle = '(deg)'
    yminor = 5
    ytickv = [-1,0,1]*floor(max(yrange)/5)*5
    yoffset = round(mean(minmax(alpha))-mean(yrange))
    yrange += yoffset
    ytickv += yoffset

    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    oplot, times, alpha
    oplot, times, theta_bg, color=second_color

    ; Add box.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        position=tpos, /nodata, /noerase

    ; Add legend.
    tx0 = tpos[2]-xchsz*0.5
    ty0 = tpos[1]+ychsz*0.3
    alignment = 1
    tx = tx0
    ty = ty0+ychsz*full_ychsz*1
    msg = tex2str('alpha')+' = arcsin(B!Dz!N/|B|), tilt angle'
    xyouts, tx,ty,/normal, msg, charsize=label_size, alignment=alignment
    ty = ty0+ychsz*full_ychsz*0
    msg = 'Background of '+tex2str('alpha')
    xyouts, tx,ty,/normal, msg, charsize=label_size, color=second_color, alignment=alignment


;---Panel 2. The detrended theta.
    tpos = poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xtickformat = '(A1)'
    yrange = minmax(theta)
    yrange = [-1,1]*max(abs(make_bins(yrange,2)))
    yticks = 2
    ytitle = '(deg)'
    yminor = 5
    ytickv = [-1,0,1]*floor(max(yrange)/5)*5
    second_color = sgcolor('red')
    third_color = sgcolor('blue')

    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    xxs = times
    yys = theta
    oplot, xxs, yys;, color=sgcolor('gray')
    plots, xrange, [0,0], linestyle=1
    oplot, times, theta_smooth, color=second_color

    ; Plot height.
    txs = df.time_range
    foreach tx,txs do plots, tx+[0,0], yrange, linestyle=1
    tys = df.theta_range
    foreach ty,tys do plots, txs,ty+[0,0], linestyle=1
    ;foreach tx,tys, kk do plots, txs[kk],tys[kk],/data, psym=1, color=third_color, symsize=0.5

    tx0 = (convert_coord(df.time_range[0],yrange[0], /data, /to_normal))[0]-xchsz*0.5
    txs = tx0+[0,0]
    tys = df.theta_range
    foreach ty,tys, kk do begin
        tmp = convert_coord(xrange[0],ty, /data, /to_normal)
        tys[kk] = tmp[1]
    endforeach
    plots, txs,tys,/normal
    for kk=0,1 do plots, txs[kk]+[-1,1]*xchsz*0.2, tys[kk]+[0,0], /normal
    ty = mean(tys)-ychsz*half_ychsz*label_size
    tx = tx0-xchsz*0.5*label_size
    xyouts, tx,ty,/normal, tex2str('Delta')+tex2str('theta'), charsize=label_size, alignment=1

    ; Add box.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        position=tpos, /nodata, /noerase

    ; Add legend.
    tx0 = tpos[2]-xchsz*0.5
    ty0 = tpos[1]+ychsz*0.3
    alignment = 1
    tx = tx0
    ty = ty0+ychsz*full_ychsz*1
    msg = tex2str('theta')+': detrended tilt angle'
    xyouts, tx,ty,/normal, msg, charsize=label_size, alignment=alignment
    ty = ty0+ychsz*full_ychsz*0
    msg = tex2str('theta')+'!Dsmth!N: '+tex2str('theta')+' smoothed over 2 min'
    xyouts, tx,ty,/normal, msg, charsize=label_size, color=second_color, alignment=alignment


;---Panel 3. dtheta/dt.
    tpos = poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    xtickformat = ''
    yrange = [-0.05,0.15]
    yticks = 2
    ytitle = '(deg/sec)'
    yminor = 5
    second_color = sgcolor('green')

    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Plot dtheta.
    xxs = times
    yys = dtheta
    oplot, xxs, yys
    plots, xrange, [0,0], linestyle=1

    ; Plot stddev.
    foreach sign, [1] do oplot, xxs, sign*dtheta_stddev, color=second_color


    ; Plot width.
    txs = df.time_range
    foreach tx,txs do plots, tx+[0,0], yrange, linestyle=1
    foreach tx,txs, kk do begin
        tmp = convert_coord(tx,yrange[1], /data, /to_normal)
        txs[kk] = tmp[0]
    endforeach
    tys = tpos[3]-ychsz*0.3+[0,0]
    plots, txs, tys, /normal
    for ii=0,1 do plots, txs[ii]+[0,0], tys[ii]+[-1,1]*ychsz*0.1, /normal
    xyouts, txs[1]+xchsz*0.5, tys[0]-ychsz*half_ychsz*label_size, tex2str('Delta')+'T', /normal, charsize=label_size


    ; Add ramp time.
    ty = tpos[1]+ychsz*3
    tx = (convert_coord(df.obs_time, yrange[1], /data, /to_normal))[0]
    plots, tx,ty,/normal, psym=1, symsize=0.5
    plots, tx+[0,0],tpos[[1,3]],/normal, linestyle=1
    xyouts, tx+xchsz*1.5*label_size, ty-ychsz*0.25*label_size, /normal, $
        'ramp time', alignment=0, charsize=label_size


    ; Add box.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        position=tpos, /nodata, /noerase

    ; Add legend.
    tx0 = tpos[2]-xchsz*0.5
    ty0 = tpos[3]-ychsz*0.9
    alignment = 1
    tx = tx0
    ty = ty0-ychsz*full_ychsz*0
    msg = 'd'+tex2str('theta')+'!Dsmth!N/dt'
    xyouts, tx,ty,/normal, msg, charsize=label_size, alignment=alignment
    ty = ty0-ychsz*full_ychsz*1
    msg = 'Local stddev over 2 min'
    xyouts, tx,ty,/normal, msg, charsize=label_size, color=second_color, alignment=alignment


    if keyword_set(test) then stop
    sgclose
end
