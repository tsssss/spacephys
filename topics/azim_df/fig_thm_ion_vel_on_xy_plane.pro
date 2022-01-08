;+
; Plot themis velocity in SM in km/s around the time of DF obs_time
;-


test = 0
reset = 1
    if n_elements(project) eq 0 then project = azim_df_load_project()

    themis_probes = 'th'+letters('e')
    ndim = 3
    events = azim_df_find_dfgroup(project=project)
    if reset then begin
        df_list = list()
        foreach event, events do begin
            direction = (strsplit(event.region,'%',/extract))[1]
            ;if direction eq 'earthward' then continue
            ;if direction eq 'outward' then continue
            foreach df, event.df_list do begin
                probe = df.probe
                index = where(themis_probes eq probe, count)
                if count eq 0 then continue
                ;time_range = df.obs_time+[-1,1]*df.width*2.5
                time_range = df.obs_time+[-1,1]*7.5*60
                the_probe = strmid(probe,2,1)
                themis_read_ion_vel, time_range, probe=the_probe
                the_var = probe+'_u_gsm'
                u_gsm = get_var_data(the_var, in=time_range, times=times)
                if n_elements(u_gsm) eq 0 then continue
                index = where(finite(u_gsm[*,0]), ntime)
                if ntime eq 0 then continue
                u_sm = cotran(u_gsm[index,*], times[index], 'gsm2sm')
                df.u_sm = total(u_sm,1)/ntime
                df_list.add, df
            endforeach
        endforeach
    endif

    ndf = df_list.length
    mlts = fltarr(ndf)
    rxys = fltarr(ndf)
    usms = fltarr(ndf,ndim)
    foreach df, df_list, ii do begin
        mlts[ii] = df.obs_mlt
        rxys[ii] = df.obs_rxy
        usms[ii,*] = df.u_sm
    endforeach




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


    file = join_path([project.plot_dir, 'fig_ion_vel_on_xy_plane.pdf'])
    if keyword_set(test) then file = 0
    magnify = keyword_set(test)? 2:1
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch, magnify=magnify



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


;---Add arrow.
    nvertex = 3
    dis_scale = 1.5
    vel_scale = 20
    hsize = keyword_set(test)? 3: 80
    arrow_solid = 0
    fill_alpha = 0.5
    fill_color = sgcolor('beige')
    fill_color2 = sgcolor('bisque')
    fill_index = [0,1,2,0]
    fill_psym = 8
    fill_symsize = 0.5
    tmp = 20
    tmp = findgen(tmp)/(tmp-1)*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    usersym, xs, ys, /fill
    
    ; Add scale.
    scale_x = tpos[2]-xchsz*2
    scale_y = tpos[3]-ychsz*1
    tx = scale_x
    ty = scale_y-ychsz*0.5
    tmp = convert_coord(tx,ty, /normal, /to_data)
    tmp = convert_coord(tmp[0]-dis_scale,tmp[1], /data, /to_normal)
    dx = abs(tmp[0]-tx)
    xs = tx+[-1,0]*dx
    plots, xs, ty+[0,0], /normal
    foreach tx, xs do plots, tx+[0,0], ty+[-1,1]*ychsz*0.1, /normal
    tx = mean(xs)
    ty = scale_y
    xyouts, tx,ty,/normal, alignment=0.5, sgnum2str(vel_scale)+' km/s'



    for df_id=0,ndf-1 do begin
        mlt = mlts[df_id]
        rxy = rxys[df_id]
        angle = (mlt*15+180)*constant('rad')
        xsm = rxy*cos(angle)
        ysm = rxy*sin(angle)

        vel = usms[df_id,0:1]
        v_hat = sunitvec(vel)
        v_mag = snorm(vel)
        if v_mag le 10 then continue
        ;v_mag = vel_scale

        x0 = xsm
        y0 = ysm
        scale = dis_scale/vel_scale*v_mag
        x1 = x0+v_hat[0]*scale
        y1 = y0+v_hat[1]*scale
        arrow, x0,y0,x1,y1,/data, solid=arrow_solid, hsize=hsize, color=color
    endfor
    
    sgclose

end
