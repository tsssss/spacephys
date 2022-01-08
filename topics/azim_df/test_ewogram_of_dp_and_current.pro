
long_time_range = time_double(['2013-06-07/02:30','2013-06-07/07:00'])
sorted_probes = ['tha','thd','the','g13','rbspa','rbspb','g15']

;long_time_range = time_double(['2008-03-14/05:00','2008-03-14/08:00'])
;sorted_probes = ['tha','thd','the']


    if n_elements(project) eq 0 then project = azim_df_load_project()
    azim_df_load_basic_data, project=project
    mlt_range = [-9,9]
    probe_xys = [-1,1]
    asi_mlt_range = [-1.5,2]

    line_linestyle = 2
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 1
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1
    rxy_range = [4.,30]
    root_dir = join_path([homedir(),'Downloads'])
    ofn = join_path([root_dir,'fig_ewogram_of_dp_and_current.pdf'])

test = 1


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ; Get the sizes of the left panels.
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1,1]
    nypanel = n_elements(ypans)
    ypad = 0.4
    panel_x_ratio = 1 ; inch per hour.
    panel_xsize = abs(total(long_time_range*[-1,1]))/3600*panel_x_ratio
    margins = [10.,4,8,1]
    fig_xsize = panel_xsize+(total(margins[[0,2]]))*abs_xchsz
    fig_ysize = panel_ysize*total(ypans)+(total(margins[[1,3]])+ypad*(nypanel-1))*abs_ychsz


;---Load data.
    if check_if_update('ae',long_time_range) then omni_read_index, long_time_range


;---Common x-axis.
    xrange = long_time_range
    xticklen_chsz = -0.2
    yticklen_chsz = -0.4
    xminor = 6 ; hr.
    xstep = 3600.
    xtickv = make_bins(xrange, xstep, /inner) ; make time line up.
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    secofday = 86400d
    for ii=0, xticks do begin
        the_time = xtickv[ii]
        xtickn[ii] = time_string(the_time,tformat='hh:mm')
        date = time_string(the_time,tformat='YYYY-MM-DD')
        if ii eq 0 then begin
            xtickn[ii] += '!C'+date
            continue
        endif
        if the_time mod secofday ne 0 then continue
        xtickn[ii] += '!C'+date
    endfor


    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz
    poss = sgcalcpos(nypanel, ypans=ypans, ypad=ypad, margins=margins)


    ; Label.
    fig_labels = letters(nypanel)+'. '+['AE','Dipolari-!C     zation','Upward!C    current']
    for ii=0,nypanel-1 do begin
        tpos = poss[*,ii]
        tx = label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, fig_labels[ii]
    endfor



;---AE.
    the_var = 'ae'
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xtickformat = '(A1)'
    ystep = 1000.
    constants = [500]
    yys = get_var_data(the_var, in=long_time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    yminor = 5
    ytitle = '(nT)'

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Add data.
    oplot, xxs, yys

    ; Add grid.
    grid_color = sgcolor('gray')
    grid_linestyle = 1
    y_constants = [1000]
    x_constants = make_bins(long_time_range, 3600)
    x_constants = x_constants[1:-2]
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


;---Panel b. UT-MLT long.
    tpos = poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    ytitle = 'MLT (hr)'
    yrange = mlt_range
    yminor = 3
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0

    data_var_suffix = 'scaled_theta'
    ztitle = 'Scaled '+tex2str('theta')+' (deg)'
    tilt_range = [-1,1]*64
    zrange = tilt_range
    ztickv = [-64,-16,-4,0,4,16,64]
    ztickn = string(ztickv,format='(I0)')
    ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
    zrange = alog10(abs(zrange))*[-1,1]
    zticks = n_elements(ztickv)-1

    ; Symbol.
    psym = 8
    symsize = 0.5
    usersym, [1,1,-1,-1,1]*0.1, [-1,1,1,-1,-1]*1
    ct = 70
    reverse_ct = 1

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase


    time_step = 10.
    times = make_bins(long_time_range+[0,-1]*time_step, time_step)
    ntime = n_elements(times)
    thetas = fltarr(ntime,nprobe)
    mlts = fltarr(ntime,nprobe)
    rxys = fltarr(ntime,nprobe)
    rsms = fltarr(ntime,3,nprobe)
    foreach probe, sorted_probes, probe_id do begin
        prefix = probe+'_'
        the_var = prefix+'scaled_theta'
        if tnames(the_var) eq '' then continue
        thetas[*,probe_id] = get_var_data(the_var, at=times)
        mlts[*,probe_id] = get_var_data(prefix+'pseudo_mlt', at=times)
        rsms[*,*,probe_id] = get_var_data(prefix+'r_sm', at=times)
        rxys[*,probe_id] = snorm(reform(rsms[*,0:1,probe_id]))
    endforeach


    foreach time, times, time_id do begin
        zzs = thetas[time_id,*]
        yys = mlts[time_id,*]
        rxy = rxys[time_id,*]
        rsm = rsms[time_id,*,*]

        index = where(finite(zzs),count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = azim_df_normalize_theta(zzs[index], zrange=tilt_range, ct=ct, /reverse_ct)
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Filter spatially.
        index = lazy_where(yys, '[]', yrange, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = lazy_where(rxy, '[]', rxy_range, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = where(check_if_in_magn(rsm, dynamic_pressure=pdyn) eq 1, count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Sort by tilt.
        index = reverse(sort(zzs))
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        for ii=0, count-1 do plots, time, yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach
        

    ; Add grid.
    y_constants = [0,3,-3]
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle


    ; Add box.
    xtickformat = '(A1)'
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Draw color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;---Panel c, EWOgram of current.
    tpos = poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz

    j_ewo_var = 'thg_j_up_ewo'
    if check_if_update(j_ewo_var) then begin
        themis_read_upward_current_ewo, long_time_range, mlt_range=mlt_range
    endif
    j_ewo_var1 = j_ewo_var+'1'
    copy_data, j_ewo_var, j_ewo_var1
    get_data, j_ewo_var1, xx,yy,zz
    yy *= 1e-3
    store_data, j_ewo_var1, xx,yy,zz
    options, j_ewo_var1, 'ztitle', 'Upward current (kA)'
    options, j_ewo_var1, 'zrange', [0,50]
    options, j_ewo_var1, 'zcharsize', 0.8
    options, j_ewo_var1, 'zposition', cbpos
    j_ewo_color_table = 62
    options, j_ewo_var1, 'xticklen', xticklen
    options, j_ewo_var1, 'yticklen', yticklen
    options, j_ewo_var1, 'color_table', j_ewo_color_table
    tplot, j_ewo_var1, trange=long_time_range, position=tpos, /noerase

    yrange = mlt_range
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, noerase=1, nodata=1

    ; Add grid.
    grid_color = sgcolor('gray')
    grid_linestyle = 1
    y_constants = [0,3,-3]
    x_constants = make_bins(long_time_range, 3600)
    x_constants = x_constants[1:-2]
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle



    if keyword_set(test) then stop
    sgclose
end
