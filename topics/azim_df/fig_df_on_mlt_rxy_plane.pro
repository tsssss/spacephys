;+
; Plot properties like event extent, DF width on the MLT-Rxy plane.
;-



test = 1
    magnify = (keyword_set(test))? 2: 1

    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_find_dfgroup(project=project)


;test_time = time_double('2014-08-28/10:30')
;if keyword_set(test_time) then begin
;    foreach event, events, event_id do begin
;        if product(event.time_range-test_time) lt 0 then break
;    endforeach
;    events = list(events[event_id])
;endif


    fig_xsize = 6
    fig_ysize = 5
    plot_file = join_path([project.plot_dir,'fig_df_on_mlt_rxy_plane.pdf'])
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize
    margins = [8,5,2,2]

    xticklen_chsz = -0.30
    yticklen_chsz = -0.40

    tpos = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz, margins=margins)
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xrange = [-1,1]*9   ; MLT.
    yrange = [4,30]     ; Re.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, nodata=1, noerase=1

    ncirc = 21
    angles = findgen(ncirc)/(ncirc-1)*2*!dpi
    circ_x = cos(angles)
    circ_y = sin(angles)
    usersym, circ_x, circ_y, /fill


    arrow_hsize = keyword_set(test)? 4:80
    arrow_solid = 1
    ndim = 3
    shape_color = sgcolor('light_cyan')
    rad = constant('rad')
    deg = constant('deg')
    foreach event, events do begin
        current_info = dictionary()

        ; Only check the azim DFs.
        direction = (strsplit(event.region,'%',/extract))[1]
        if direction eq 'earthward' then continue
        if direction eq 'outward' then continue

    ;---Sort based on obs_time.
        df_list = event.df_list
        ndf = df_list.length
        obs_times = dblarr(ndf)
        foreach df, df_list, ii do obs_times[ii] = df.obs_time
        index = sort(obs_times)
        obs_times = obs_times[index]
        df_list = df_list[index]
        probes = strarr(ndf)
        foreach df, df_list, ii do probes[ii] = df.probe
        sorted_probes = probes
        obs_mlts = fltarr(ndf)
        foreach df, df_list, ii do obs_mlts[ii] = df.obs_mlt
        obs_rxys = fltarr(ndf)
        foreach df, df_list, ii do obs_rxys[ii] = df.obs_rxy
        df_widths = fltarr(ndf)
        foreach df, df_list, ii do df_widths[ii] = df.width
        df_mean_width = mean(df_widths)

    ;---v_2d and omega.
        triad_list = event.triad_list
        vmags = fltarr(triad_list.length)
        foreach triad, triad_list, ii do vmags[ii] = triad.vmag_obs_time
        mean_vmag = mean(vmags)

        omegas = fltarr(triad_list.length)
        foreach triad, triad_list, ii do omegas[ii] = triad.omega_obs_time
        omega_2d = mean(omegas)
        omega_2d_error = stddev(omegas)

    ;---Properties.
        rxy_extent = minmax(obs_rxys)
        mlt_extent = minmax(obs_mlts)
        angle_extent = mlt_extent*15
        if obs_mlts[0] gt obs_mlts[ndf-1] then angle_extent = reverse(angle_extent)
        mean_rxy = mean(rxy_extent)
        mean_angle = mean(angle_extent)

    ;---Width.
        ww_angle = abs(omega_2d)*min(df_widths)/60
        ww_deg = ww_angle*deg


    ;---Plot.
        ; Shape.
        shape_color = sgcolor('light_cyan')
        rr = 0
        the_angle = angle_extent[0]*(1-rr)+angle_extent[1]*rr
        txs = (the_angle+[-1,1]*ww_angle)/15
        tys = rxy_extent
        txs = txs[[0,1,1,0,0]]
        tys = tys[[0,0,1,1,0]]
        ;polyfill, txs, tys, /data, color=shape_color
        plots, txs, tys

        ; MLT extent.
        rr = 1
        the_rxy = rxy_extent[0]*(1-rr)+rxy_extent[1]*rr
        txs = angle_extent/15
        tys = the_rxy+[0,0]
        plots, txs, tys, /data
        arrow, txs[-2],tys[-2], txs[-1],tys[-1], /data, hsize=arrow_hsize*1.5, solid=arrow_solid
        tmp = convert_coord(txs[0],tys[0], /data, /to_normal)
;        txs = tmp[0]+[0,0]
;        tys = tmp[1]+[-1,1]*ychsz*0.25
;        plots, txs,tys, /normal
        plots, txs[0],tys[0], psym=8
    endforeach

    xtitle = 'MLT (hr)'
    ytitle = 'Rxy (Re)'
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticklen=xticklen, xtitle=xtitle, $
        ystyle=1, yrange=yrange, yticklen=yticklen, ytitle=ytitle, $
        position=tpos, nodata=1, noerase=1

    stop


    if keyword_set(test) then stop
    sgclose


end
