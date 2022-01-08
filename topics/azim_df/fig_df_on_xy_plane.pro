;+
; Plot properties like event extent, DF width on the MLT-Rxy plane.
;-



test = 0
    magnify = (keyword_set(test))? 2: 1

    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_find_dfgroup(project=project)


    rsm_xrange = [10,-30]
    rsm_yrange = [1,-1]*20
    xsize = abs(total(rsm_xrange*[-1,1]))
    ysize = abs(total(rsm_yrange*[-1,1]))
    fig_xsize = 5
    margins = [8.,5,2,2]
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    xpanel_size = fig_xsize-total(margins[[0,2]])*abs_xchsz
    ypanel_size = xpanel_size/xsize*ysize
    fig_ysize = ypanel_size+total(margins[[1,3]])*abs_ychsz

    tpos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    tpos = [tpos[0:1],1-tpos[2:3]]

    plot_file = join_path([project.plot_dir,'fig_df_on_xy_plane.pdf'])
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize

    xticklen_chsz = -0.30
    yticklen_chsz = -0.40

;---Plot settings.
    xticklen_chsz = -0.15
    yticklen_chsz = -0.40
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    label_size = constant('label_size')

    xrange = rsm_xrange
    xstep = 10
    xminor = 5
    xtickv = make_bins(xrange,xstep,/inner)
    xticks = n_elements(xtickv)-1
    xlog = 0
    xtitle = 'SM X (Re)'

    yrange = rsm_yrange
    ystep = 10
    yminor = 5
    ytickv = make_bins(yrange,ystep,/inner)
    yticks = n_elements(ytickv)-1
    ylog = 0
    ytitle = 'SM Y (Re)'

    ; Set up coord.
    plot, xrange,yrange, /nodata, /noerase, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, ytitle=ytitle, yticklen=yticklen, $
        position=tpos

    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    plots, circle_x, circle_y

    ; Add ROI and dash lines.
    pdyn = 10.
    rxy_range = [4,30]
    mlt_range = [-1,1]*9
    rad = constant('rad')
    deg = constant('deg')
    angle_range = (mlt_range*15+180)*rad
    roi_color = sgcolor('gray')
    roi_thick = keyword_set(test)? 2: 8

    magn_test = fltarr(nangle,3)
    magn_test[*,0] = circle_x*30
    magn_test[*,1] = circle_y*30
    tmp = check_if_in_magn(magn_test, magn_pos=magn_pos, dynamic_pressure=pdyn)
    magn_pos = magn_pos[*,[0,2]]
    magn_rxy = snorm(magn_pos)
    magn_angle = atan(magn_pos[*,1],magn_pos[*,0])
    p_mlt = sinterpol(magn_pos, magn_angle, min(angle_range))
    p_rxy = sinterpol(magn_pos, magn_rxy, max(rxy_range))

    ; Magnetopause dash line.
    foreach sector_flag, [-1,1] do begin
        ; -1: pre-midn; 1: post-midn.
        txs = magn_pos[*,0]
        tys = magn_pos[*,1]*sector_flag
        oplot, txs,tys, linestyle=1
    endforeach

    ; MLT dash line.
    rs = [1,max(rxy_range)]
    foreach t, [angle_range,!dpi] do begin
        oplot, cos(t)*rs, sin(t)*rs, linestyle=1
    endforeach

    ; Rxy dash line.
    foreach r, make_bins(rxy_range,10,/inner) do begin
        oplot, circle_x*r, circle_y*r, linestyle=1
    endforeach


    ; Add label.
    label_x = tpos[0]+xchsz*1
    label_y = tpos[3]-ychsz*1
    msg = 'Region of interest (ROI)!C(1) Magnetopause @10 nPa;!C(2) Rxy [4,30] Re;!C(3) MLT [-9,9] hr'
    xyouts, label_x,label_y,/normal, msg

    ; Magnetopause.
    index = where(magn_rxy le max(rxy_range) and magn_angle ge min(angle_range))
    magn_pos = magn_pos[index,*]
    magn_xs = [magn_pos[*,0],p_mlt[0],p_rxy[0]]
    magn_ys = [magn_pos[*,1],p_mlt[1],p_rxy[1]]
    index = sort(magn_ys)
    magn_xs = magn_xs[index]
    magn_ys = magn_ys[index]

    foreach sector_flag, [-1,1] do begin
        ; -1: pre-midn; 1: post-midn.
        txs = magn_xs
        tys = magn_ys*sector_flag
        oplot, txs,tys, color=roi_color, thick=roi_thick
    endforeach

    ; MLT.
    mlt_rxy = snorm(p_mlt)
    rs = [min(rxy_range),mlt_rxy]
    foreach t, angle_range do begin
        oplot, cos(t)*rs, sin(t)*rs, color=roi_color, thick=roi_thick
    endforeach
    ; Rxy.
    rxy_angle = atan(p_rxy[1],p_rxy[0])
    foreach r, rxy_range do begin
        if r eq min(rxy_range) then begin
            ts = smkarthm(angle_range[0],angle_range[1],nangle, 'n')
        endif else begin
            ts = smkarthm(rxy_angle,2*!dpi-rxy_angle,nangle, 'n')
        endelse
        oplot, cos(ts)*r, sin(ts)*r, color=roi_color, thick=roi_thick
    endforeach


    ct = 64
    nevent = events.length
    index_colors = smkarthm(100,250,nevent,'n')
    colors = fltarr(nevent)
    for ii=0,nevent-1 do colors[ii] = sgcolor(index_colors[ii], ct=ct)


    arrow_hsize = keyword_set(test)? 4:80
    arrow_solid = 1
    ndim = 3
    shape_color = sgcolor('light_cyan')
    rad = constant('rad')
    deg = constant('deg')
    foreach event, events, event_id do begin
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
        angle_extent = (mlt_extent*15+180)*rad
        if obs_mlts[0] gt obs_mlts[ndf-1] then angle_extent = reverse(angle_extent)
        mean_rxy = mean(rxy_extent)
        mean_angle = mean(angle_extent)

    ;---Width.
        ww_angle = abs(omega_2d)*min(df_widths)/60*rad


    ;---Plot.
        rr = 0.5
        the_angle = angle_extent[0]*(1-rr)+angle_extent[1]*rr
        angles = smkarthm(the_angle-ww_angle*0.5, the_angle+ww_angle*0.5,nangle, 'n')
        trs = smkarthm(rxy_extent[1],rxy_extent[0],nangle,'n')
        txs = [rxy_extent[1]*cos(angles),trs*cos(angles[-1]),rxy_extent[0]*reverse(cos(angles)),reverse(trs)*cos(angles[0])]
        tys = [rxy_extent[1]*sin(angles),trs*sin(angles[-1]),rxy_extent[0]*reverse(sin(angles)),reverse(trs)*sin(angles[0])]
        ;polyfill, txs, tys, /data, color=shape_color
        plots, txs,tys, color=colors[event_id]


;        ; MLT extent.
;        rr = 0.8
;        the_rxy = rxy_extent[0]*(1-rr)+rxy_extent[1]*rr
;        angles = smkarthm(angle_extent[0], angle_extent[1], nangle, 'n')
;        txs = the_rxy*cos(angles)
;        tys = the_rxy*sin(angles)
;        plots, txs, tys, /data, color=colors[event_id]
;        arrow, txs[-2],tys[-2], txs[-1],tys[-1], /data, hsize=arrow_hsize*1.5, solid=arrow_solid, color=colors[event_id]
;        tmp = convert_coord(0,ychsz, /normal,/to_data)-convert_coord(0,0, /normal,/to_data)
;        drxy = tmp[1]*0.3
;        tmp = the_rxy+[-1,1]*drxy
;        plots, tmp*cos(angles[0]), tmp*sin(angles[0]), /data, color=colors[event_id]

    endforeach


    if keyword_set(test) then stop
    sgclose


end
