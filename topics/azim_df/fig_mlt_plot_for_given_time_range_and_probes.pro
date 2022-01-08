;+
; Load data and generate plot showing DF prop for a given time range and probes.
;-



    test = 1
    ; Angelopoulos event.
    short_time_range = time_double(['2007-03-23/11:00','2007-03-23/12:00'])
    probes = 'th'+letters('e')
    ; Ogasawara Figure 3 event.
;    short_time_range = time_double(['2008-02-26/03:50','2008-02-26/04:50'])
;    probes = 'th'+['a','d','e']
    ; Ogasawara Figure 6 event.
;    short_time_range = time_double(['2008-03-03/07:50','2008-03-03/08:40'])
;    probes = 'th'+['d','e']
    
    nprobe = n_elements(probes)


;---Check inputs.
    if n_elements(short_time_range) ne 2 then begin
        errmsg = handle_error('No input short_time_range ...')
        return
    endif

    if nprobe eq 0 then begin
        errmsg = handle_error('No input probes ...')
        return
    endif

;---Load data.
    mean_time = mean(short_time_range)
    event_id = time_string(mean_time,tformat='YYYY_MMDD_hhmm')
    common_data_rate = 10.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    local_root = project.data_dir
    data_file = join_path([local_root,'data','manual_event','fig_df_prop_'+event_id+'.tplot'])

    if file_test(data_file) eq 0 then begin
        vars = ['r_sm','theta']
        save_vars = list()
        foreach probe, probes do begin
            prefix = probe+'_'
            data_time_range = mean_time+[-1,1]*3600*5
            foreach var, vars do begin
                azim_df_read_data, var, time_range=data_time_range, probe=probe, project=project
                the_var = prefix+var
                uniform_time, the_var, common_data_rate
                save_vars.add, the_var
            endforeach

            mlt_var = prefix+'mlt'
            r_sms = get_var_data(prefix+'r_sm', times=times)
            mlt = azim_df_calc_pseudo_mlt(r_sms)
            store_data, mlt_var, times, mlt
            add_setting, mlt_var, /smart, {$
                unit: 'hr', $
                display_type: 'scalar', $
                short_name: 'MLT'}
            save_vars.add, mlt_var
        endforeach

        ; AE and Dst.
        omni_read_index, data_time_range
        save_vars.add, ['ae','dst'], /extract

        save_vars = save_vars.toarray()
        tplot_save, save_vars, filename=data_file
        lprmsg, 'Data saved to '+data_file+' ...'
    endif
    if file_test(data_file) then tplot_restore, filename=data_file


;---Plot settings.
    mlts = dblarr(nprobe)
    foreach probe, probes, ii do mlts[ii] = get_var_data(probe+'_mlt', at=mean_time)
    index = sort(abs(mlts))
    sorted_probes = probes[index]
    probe_colors = smkarthm(50,250,nprobe, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)

    mlt_range = []
    foreach probe, probes, ii do mlt_range = [mlt_range, minmax(get_var_data(probe+'_mlt', in=short_time_range))]
    mlt_range = minmax(make_bins(mlt_range, 1))
    rx_range = []
    ry_range = []
    foreach probe, probes, ii do begin
        rsm = get_var_data(probe+'_r_sm', in=short_time_range)
        rx_range = [rx_range,minmax(rsm[*,0])]
        ry_range = [ry_range,minmax(rsm[*,1])]
    endforeach
    rx_range = minmax(make_bins(rx_range,1))
    ry_range = minmax(make_bins(ry_range,1))

    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 2
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=xchsz, ychsz=ychsz
    sgclose, /wdelete
    abs_xchsz = xchsz
    abs_ychsz = ychsz
    ; Get the sizes of the left panels.
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1]
    nypanel = n_elements(ypans)
    panel_ypads = [0.4]
    panel_x_ratio = 3.0     ; inch per hour.
    panel_xsize = abs(total(short_time_range*[-1,1]))/3600*panel_x_ratio
    ; Get the size of the right panel.
    right_pan_ysize = total(panel_ypads)*abs_ychsz+total(ypans)*panel_ysize
    right_pan_xsize = right_pan_ysize


    margins = [12.,5,4,2]
    panel_xpad = 14

    fig_ysize = right_pan_ysize+total(margins[[1,3]])*abs_ychsz
    fig_xsize = panel_xsize+panel_xpad*abs_xchsz+right_pan_xsize+total(margins[[0,2]])*abs_xchsz
    file = join_path([project.plot_dir,'manual_event','fig_df_prop_'+event_id+'.pdf'])
    if keyword_set(test) then file = test
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch

    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    lower_pos = pos
    poss = sgcalcpos(1,2, position=lower_pos, xpans=[panel_xsize,right_pan_xsize], xpad=panel_xpad)
    left_pos = poss[*,0]
    left_poss = sgcalcpos(nypanel, position=left_pos, ypans=ypans, ypad=panel_ypads, xchsz=xchsz, ychsz=ychsz)
    right_pos = poss[*,1]

    ;for ii=0, nypanel-1 do plot, [0,1], [0,1], position=left_poss[*,ii], /noerase, /nodata
    ;plot, [0,1], [0,1], position=right_pos, /noerase, /nodata
    ;stop

;---Common x-axis.
    xrange = short_time_range
    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    xminor = 10 ; min.
    xtickv = make_bins(xrange, 60*xminor, /inner) ; make time line up.
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



;---Panel a. Dst/AE.
    tpos = left_poss[*,0]
    xtickformat = '(A1)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    the_var = 'dst'
    ystep = 50.
    constants = [-50,0]
    yys = get_var_data(the_var, in=short_time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    ytitle = '(nT)'

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Add data.
    oplot, xxs, yys
    foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=9, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    the_var = 'ae'
    ystep = 500.
    constants = [0,500]
    yys = get_var_data(the_var, in=short_time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    ytitle = '(nT)'

    ; Set up coord.
    ae_color = sgcolor('red')
    axis, yaxis=1, /save, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        color=ae_color

    ; Add data.
    oplot, xxs, yys, color=ae_color
    foreach constant, constants do oplot, xrange, constant+[0,0], color=ae_color, linestyle=1

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'a. Dst/'
    xyouts, tx,ty,/normal, '           AE', color=ae_color



;---Panel b. UT-MLT short.
    tpos = left_poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = ''
    ytitle = 'MLT (hr)'
    yrange = mlt_range
    yminor = 4.
    ystep = 2.
    ytickv = make_bins(yrange,ystep)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0

    data_var_suffix = 'theta'
    ztitle = 'Detrended tilt '+tex2str('theta')+' (deg)'
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

    pos_var_suffix = 'mlt'
    flags = list()
    foreach probe, sorted_probes do begin
        prefix = probe+'_'
        pos_var = prefix+pos_var_suffix
        data_var = prefix+data_var_suffix
        get_data, data_var, xxs, zzs
        get_data, pos_var, xxs, yys
        if n_elements(zzs) ne n_elements(xxs) then begin
            errmsg = handle_error(data_var+': inconsistent data and time ...')
            stop
        endif
        mlt = get_var_data(prefix+'mlt')
        zzs = azim_df_scale_theta(zzs, mlt) ; apply MLT scaling.
        zzs = azim_df_normalize_theta(zzs, zrange=tilt_range, ct=ct, /reverse_ct)

        ; Remove data out of range.
        get_data, prefix+'r_sm', times, r_sm
        rx = r_sm[*,0]
        ; Exclude data outside the MLT range.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = lazy_where(rx, '][', rx_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(cotran(r_sm, times, 'sm2gsm'))
        index = where(magn_flags eq 0, count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside time range.
        index = lazy_where(xxs, '][', short_time_range, count=count)
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
    oplot, xrange, constant+[0,0], linestyle=constant_linestyle

    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'b. UT-'+strmid(ytitle,0,strpos(ytitle,' '))


;---Right panel.
    tpos = right_pos
    if n_elements(xy_range) ne 2 then xy_range = [5,-15]  ; sun to the left.
    xstep = 5
    ystep = 5
    panel_xrange = minmax(make_bins([rx_range,1,-1],xstep))
    panel_yrange = minmax(make_bins([ry_range,1,-1],ystep))
    xrange = reverse(panel_xrange)
    xtickv = make_bins(panel_xrange,xstep)
    xticks = n_elements(xtickv)
    xminor = xstep
    yrange = reverse(panel_yrange)
    ytickv = make_bins(panel_yrange,ystep)
    yticks = n_elements(ytickv)
    yminor = ystep
    
    xtitle = 'SM X (Re)'
    ytitle = 'SM Y (Re)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    plot, xrange, yrange, /nodata, /noerase, /isotropic, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, ytitle=ytitle, yticklen=yticklen, $
        position=tpos


    ; Add figure label.
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1.2
    label = 'd. SM x-y plane'
    xyouts, tx,ty,/normal,alignment=0, label

    ; Add earth.
    tmp = 50
    tmp = findgen(tmp)/(tmp-1)*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    polyfill, xs<0, ys, /line_fill, orientation=45, spacing=ychsz*2
    plots, xs, ys
    foreach r, make_bins(minmax(abs(xy_range)),5,/inner) do oplot, xs*r, ys*r, linestyle=1
    
    foreach probe, sorted_probes, ii do begin
        r_sm = get_var_data(probe+'_r_sm', in=short_time_range)
        xxs = r_sm[*,0]
        yys = r_sm[*,1]
        plots, xxs, yys, color=probe_colors[ii]
        tx = median(xxs)
        ty = median(yys)
        xyouts, tx,ty,/data, strupcase(probe), color=probe_colors[ii], charsize=0.7
    endforeach

    if keyword_set(test) then stop
    sgclose

end
