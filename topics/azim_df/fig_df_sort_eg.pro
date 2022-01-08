;+
; A example to sort DF.
;-



;---Test event times, prepare DF and the tilt angle data.
    plot_time_range = time_double(['2014-08-28/09:10:00','2014-08-28/13:10:00'])
    region_name = 'post_midn'

    project = azim_df_load_project()
    file = join_path([project.data_dir,'azim_df_search_'+region_name+'_search_vertex.txt'])
    if ~file_test(file) then stop
    events = azim_df_subgroup_read_file(file)

    df_list = list()
    foreach event, events do begin
        foreach df, event.df_list do begin
            if product(df.obs_time-plot_time_range) gt 0 then continue
            df_list.add, df
        endforeach
    endforeach

    df_probes = list()
    foreach df, df_list do begin
        the_probe = df.probe
        if df_probes.where(the_probe) ne !null then continue
        df_probes.add, the_probe
    endforeach

    azim_df_load_basic_data, project=project



test = 0
    plot_file = join_path([project.plot_dir,'fig_df_sort_eg.pdf'])
    if keyword_set(test) then plot_file = test

;---Plot settings.
    xtickformat = ''
    xticklen_chsz = -0.15   ; in ychsz.
    yticklen_chsz = -0.30   ; in xchsz.
    full_ychsz = constant('full_ychsz')
    half_ychsz = full_ychsz*0.5
    lineskip = constant('lineskip')
    label_size = 0.7
    label_xshift = 10
    label_yshift = full_ychsz
    secofday = constant('secofday')
    bar_thick = keyword_set(test)? 0.5: 4

    rad = constant('rad')
    deg = constant('deg')

    mlt_range = [-1,1]*9
    rxy_range = [4.,30]


    panel_x_ratio = 1.     ; inch per hour.
    panel_xsize = abs(total(plot_time_range*[-1,1]))/3600*panel_x_ratio
    panel_ysize = 2
    margins = [8,4,15,4]

    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=xchsz, ychsz=ychsz
    sgclose, /wdelete
    abs_xchsz = xchsz
    abs_ychsz = ychsz
    fig_xsize = panel_xsize+total(margins[[0,2]])*abs_xchsz
    fig_ysize = panel_ysize+total(margins[[1,3]])*abs_ychsz
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch

    pos_list = list()
    pos_list.add, sgcalcpos(1, xchsz=xchsz, ychsz=ychsz, margins=margins)


    tpos = pos_list[0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

;---The common x-axis setting.
    xrange = plot_time_range
    xstep = constant('secofhour')
    if total(xrange*[-1,1]) le xstep then xstep = 600.
    xticks = floor(total(xrange*[-1,1])/xstep)
    xtickv = smkarthm(xrange[0], xstep, xticks+1, 'x0')
    xminor = 6
    xtickn = strarr(xticks+1)
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

;---Common z-range, colorbar, symbol.
    theta_var = 'theta'
    ztitle = 'Scaled detr.tilt '+tex2str('theta')+' (deg)'
    theta_range = [-1,1]*64
    zrange = theta_range
    ztickv = [-64,-16,-4,0,4,16,64]
    ztickn = string(ztickv,format='(I0)')
    ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
    zrange = alog10(abs(zrange))*[-1,1]
    zticks = n_elements(ztickv)-1
    spec_psym = 8
    spec_symsize = 0.3
    usersym, [1,1,-1,-1,1]*0.5, [-1,1,1,-1,-1]*1, /fill
    spec_ct = 70


;---Common x data.
    time_step = project.time_step
    xxs = make_bins(xrange, time_step)
    nxx = n_elements(xxs)

;---MLT-UT.
    pos_var = 'pseudo_mlt'
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    ytitle = 'MLT (hr)'
    yminor = 3
    yrange = []
    foreach df, df_list do yrange = [yrange,df.obs_mlt]
    yrange = [floor(min(yrange)),ceil(max(yrange))]+[-1,1]*1.5
    ytickv = make_bins(yrange, yminor)
    if n_elements(ytickv) le 1 then ytickv = smkarthm(yrange[0],yrange[1],yminor, 'dx')
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    constants = 0


    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Plot data.
    foreach probe, df_probes do begin
        prefix = probe+'_'
        yys = get_var_data(prefix+pos_var, at=xxs)
        zzs = get_var_data(prefix+theta_var, at=xxs)
        mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
        zzs = azim_df_scale_theta(zzs, mlt) ; apply MLT scaling.
        zzs = azim_df_normalize_theta(zzs, zrange=theta_range, ct=spec_ct, /reverse_ct)

        ; Remove data outside ROI.
        index = lazy_where(mlt, '][', yrange, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        rsm = get_var_data(prefix+'r_sm', at=xxs)
        rxy = snorm(rsm[*,0:1])
        index = lazy_where(rxy, '][', rxy_range, count=count)
        if count ne 0 then yys[index] = !values.f_nan
        index = where(finite(zzs,/nan), count)
        if count ne 0 then yys[index] = !values.f_nan

        ; Plot data.
        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=spec_psym, symsize=spec_symsize
    endforeach


    ; Add large DF.
    foreach df, df_list do begin
        ty = df.obs_mlt
        tx = df.obs_time
        plots, tx,ty,/data, psym=6, symsize=0.5, color=sgcolor('red')
    endforeach

    ; Add notations.
    foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Add label.
    tx = tpos[0]-label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    ;xyouts, tx,ty,/normal, 'a. UT-MLT'


    line_color = sgcolor('silver')
    ; Add DF per probe.
    nprobe = df_probes.length
    probe_infos = project.probe_infos
    dy = (tpos[3]-tpos[1])/(nprobe+1)
    ty2s = tpos[1]+(findgen(nprobe)+1)*dy   ; The msg.
    msgs = list()
    foreach probe, df_probes do begin
        msgs.add, strupcase(probe_infos[probe].short_name)
    endforeach
    ndfs = list()
    foreach probe, df_probes do begin
        ndf = 0
        foreach df, df_list do if df.probe eq probe then ndf += 1
        ndfs.add, ndf
    endforeach
    
    ty1s = list()               ; The data.
    foreach probe, df_probes do begin
        prefix = probe+'_'
        yys = get_var_data(prefix+pos_var, times=xxs, in=plot_time_range)
        tmp = convert_coord(xxs[-1], yys[-1], /data, /to_normal)
        ty1s.add, tmp[1]
    endforeach

    ; Sort according to end the data.
    ty1s = ty1s.sort(indices=index)
    msgs = msgs[index]
    ndfs = ndfs[index]
    tx1 = tpos[2]+xchsz*1
    tx2 = tpos[2]+xchsz*4
    foreach msg, msgs, ii do begin
        xyouts, tx2+xchsz*2, ty2s[ii]-ychsz*half_ychsz, /normal, msg, alignment=0.5
        xyouts, tx2+xchsz*7, ty2s[ii]-ychsz*half_ychsz, /normal, string(ndfs[ii],format='(I0)'), alignment=0.5
        plots, [tx1+[-1,0]*xchsz,tx2], [ty1s[ii]+[0,0],ty2s[ii]], /normal, color=line_color
    endforeach
    
    ty = tpos[3]-ychsz
    xyouts, tx2+xchsz*2, ty, /normal, 'SC', alignment=0.5
    xyouts, tx2+xchsz*7, ty, /normal, '# of DP', alignment=0.5
    plots, tx2+[0,10]*xchsz, ty+[0,0]-ychsz*0.2, /normal
    
    
    ; SC sequence.
    ndf = df_list.length
    tmp = convert_coord(the_times, findgen(ndf), /data, /to_normal)
    tx1s = reform(tmp[0,*])
    tx0 = tpos[0]+xchsz*6
    dx = (tpos[2]-tx0)/(ndf+1)
    tx2s = tx0+(findgen(ndf)+1)*dx
    
    ty1 = tpos[1]+ychsz*3
    ty2 = tpos[1]+ychsz*1
    xyouts, tpos[0]+xchsz*1, ty2, /normal, 'SC seq:'
    foreach df, df_list, ii do begin
        probe = df.probe
        tmp = convert_coord(df.obs_time, df.obs_mlt, /data, /to_normal)
        tx1 = tmp[0]
        ty0 = tmp[1]
        msg = strupcase(probe_infos[probe].short_name)
        if ii ne df_list.length-1 then msg = msg+','
        xyouts, tx2s[ii],ty2,/normal, alignment=0.5, msg
        
        plots, [tx1+[0,0],tx2s[ii]], [ty0,ty1,ty2+ychsz*1], /normal, color=line_color
    endforeach

    

    ; Colorbar.
    cbpos = tpos[[0,3,2,3]]+[0,0.5,0,1]*ychsz
    sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks, horizontal=1


    if keyword_set(test) then stop
    sgclose

end
