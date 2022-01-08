
long_time_range = time_double(['2013-06-07/04:00','2013-06-07/07:00'])

    if n_elements(project) eq 0 then project = azim_df_load_project()
    azim_df_load_basic_data, project=project
    mlt_range = [-1,1]*6
    probe_xys = [-1,1]
    asi_mlt_range = [-1.5,2]

    line_linestyle = 2
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 1
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1
    sorted_probes = ['thd','the','g13','rbspa','rbspb','g15']
    rxy_range = [4.,30]
    root_dir = join_path([homedir(),'Downloads'])
    ofn = join_path([root_dir,'fig_ewogram_of_dp.pdf'])

test = 0


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ; Get the sizes of the left panels.
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1]
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
    fig_labels = letters(nypanel)+'. '+['AE','Dipolari-!C     zation']
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

    pos_var_suffix = 'pseudo_mlt'
    flags = list()
    foreach probe, reverse(sorted_probes) do begin
        prefix = probe+'_'
        pos_var = prefix+pos_var_suffix
        data_var = prefix+data_var_suffix
        zzs = get_var_data(data_var, in=long_time_range, times=xxs)
        yys = get_var_data(pos_var, at=xxs)
        mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
        zzs = azim_df_normalize_theta(zzs, zrange=tilt_range, ct=ct, /reverse_ct)

        ; Remove data out of range.
        r_sm = get_var_data(prefix+'r_sm', at=xxs)
        rxy = snorm(r_sm[*,0:1])
        ; Exclude data outside the MLT range.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = lazy_where(rxy, '][', rxy_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(r_sm)
        index = where(magn_flags eq 0, count)
        if count ne 0 then xxs[index] = !values.f_nan

        index = where(finite(xxs), count)
        the_flag = (count gt 0)
        flags.add, the_flag
        if count eq 0 then continue

        nxx = count
        xxs = xxs[index]
        yys = yys[index]
        zzs = zzs[index]
        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach

    ; Add grid.
    y_constants = [0,3,-3]
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle


    ; Add box.
    xtickformat = ''
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Draw color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks




    if keyword_set(test) then stop
    sgclose
end
