
;---Input.
    probes = ['thd','g15','1991-080','g13','1994-084','LANL-01A','LANL-02A','LANL-04A']
    pos_times = time_double('2014-08-28/'+['10:10','10:20'])
    coord = 'sm'

test = 0
    ct = 40
    top_color = 250
    bottom_color = 100
    nprobe = n_elements(probes)
    index_colors = bytscl(findgen(nprobe), top=top_color-bottom_color)+bottom_color
    probe_colors = fltarr(nprobe)
    foreach color, index_colors, color_id do $
        probe_colors[color_id] = sgcolor(color, ct=ct)

    time_colors = sgcolor(['red','orange'])
    time_psyms = [6,7,1]

    short_names = strupcase(probes)

    _2014_0828_10_load_data
    pos_coord = list()
    foreach probe, probes do begin
        r_gsm = get_var_data(probe+'_r_gsm', at=pos_times)
        r_coord = cotran(r_gsm, pos_times, 'gsm2'+coord)
        pos_coord.add, r_coord
    endforeach
    pos_coord = pos_coord.toarray()

    x_range = minmax(pos_coord[*,*,0])
    y_range = minmax(pos_coord[*,*,1])
    z_range = minmax(pos_coord[*,*,2])
    rxy_step = 2
    z_step = 1
    xrange = (minmax(make_bins([-1,1,x_range+[-1,1]],rxy_step)))
    yrange = (minmax(make_bins([-1,1,y_range+[-1,1]],rxy_step)))
    zrange = [0,3]

    xstep = 5
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    ystep = 5
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = ystep
    xrange = reverse(xrange)
    yrange = reverse(yrange)

    zstep = 2
    ztickv = make_bins(zrange, zstep, /inner)
    zticks = n_elements(ztickv)-1
    zminor = zstep

    xtitle = strupcase(coord)+' X (Re)'
    ytitle = strupcase(coord)+' Y (Re)'
    ztitle = strupcase(coord)+' Z (Re)'

    xticklen_chsz = -0.20
    yticklen_chsz = -0.40

    sgopen, 0, xsize=1, ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    xpan_size = 3
    ypad = 0.4
    ypan_size = abs(xpan_size/total(xrange*[-1,1])*total(yrange*[-1,1]))
    margins = [8,5,2,2]
    fig_xsize = xpan_size+total(margins[[0,2]])*abs_xchsz
    fig_ysize = ypan_size+total(margins[[1,3]])*abs_ychsz


    root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
    ofn = join_path([root_dir,'case2_fig_injection_pos.pdf'])
    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz

    poss = sgcalcpos(1, margins=margins, ypad=ypad)


;---Setup coord.
    tpos = poss
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    plot, xrange, yrange, /iso, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos

    plots, xrange, [0,0], linestyle=1
    plots, [0,0], yrange, linestyle=1

    ; Add earth and circles.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    polyfill, circle_x>0, circle_y, color=sgcolor('white')
    plots, circle_x, circle_y
    foreach r, [5,10] do oplot, circle_x*r, circle_y*r, linestyle=1
    
    ; Add Geosync orbit.
    r = 6.6
    oplot, circle_x*r, circle_y*r, linestyle=0, color=sgcolor('silver')
    rad = constant('rad')
    t = -0
    tx = r*cos(t*rad)
    ty = r*sin(t*rad)
    xyouts, tx,ty,/data, ' '+string(r,format='(F3.1)')+' Re', color=sgcolor('silver')
    
    
    ; Add injections.
    thick = keyword_set(test)? 4: 10

    msg = 'Electron!Cinjection 1'
    mlt_range = [0.4,4.2]
    the_color = sgcolor('salmon')
    r = 6.6
    angle_range = mlt_range*15+180
    angles = make_bins(angle_range,5)*rad
    tx = r*cos(angles)
    ty = r*sin(angles)
    oplot, tx,ty, color=the_color, thick=thick
    mean_angle = mean(angles)
    tx = r*cos(mean_angle)
    ty = r*sin(mean_angle)
    tmp = convert_coord(tx,ty,/data,/to_normal)
    tx = tmp[0]
    ty = tmp[1]
    xyouts, tx+xchsz*4,ty+ychsz*1,/normal, alignment=0.5, msg, color=the_color
    
    msg = 'Electron!Cinjection 2'
    mlt_range = [2.6,3.9]
    the_color = sgcolor('wheat')
    r = 6.0
    angle_range = mlt_range*15+180
    angles = make_bins(angle_range,1)*rad
    tx = r*cos(angles)
    ty = r*sin(angles)
    oplot, tx,ty, color=the_color, thick=thick
    mean_angle = mean(angles)
    tx = r*cos(mean_angle)
    ty = r*sin(mean_angle)
    tmp = convert_coord(tx,ty,/data,/to_normal)
    tx = tmp[0]
    ty = tmp[1]
    xyouts, tx-xchsz*4,ty,/normal, alignment=0.5, msg, color=the_color
    
    msg = 'Ion!Cinjection 1'
    mlt_range = [-5.4,0.4]
    the_color = sgcolor('light_blue')
    r = 6.6
    angle_range = mlt_range*15+180
    angles = make_bins(angle_range,1)*rad
    tx = r*cos(angles)
    ty = r*sin(angles)
    oplot, tx,ty, color=the_color, thick=thick
    mean_angle = mean(angles)
    tx = r*cos(mean_angle)
    ty = r*sin(mean_angle)
    tmp = convert_coord(tx,ty,/data,/to_normal)
    tx = tmp[0]
    ty = tmp[1]
    xyouts, tx+xchsz*3,ty,/normal, alignment=0.5, msg, color=the_color


    ; Add probes.
    line_color = sgcolor('silver')
    psym = 8
    symsize = 0.5
    nangle = 20
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    usersym, circle_x, circle_y, /fill
    foreach probe, probes, probe_id do begin
        xx = reform(pos_coord[probe_id,*,0])
        yy = reform(pos_coord[probe_id,*,1])
        zz = reform(pos_coord[probe_id,*,2])
        color = probe_colors[probe_id]
        plots, xx, yy, color=line_color
        foreach time, pos_times, time_id do begin
            ;plots, xx[time_id], yy[time_id], color=color, $
            ;    psym=psym, symsize=symsize
            plots, xx[time_id], yy[time_id], color=color, $
                psym=time_psyms[time_id], symsize=symsize
        endforeach

        tmp = max(sqrt(xx^2+yy^2), index)
        tmp = min(yy, index)
        tmp = convert_coord(xx[index], yy[index], /data, /to_normal)
        tx = tmp[0]+xchsz*0.0
        ty = tmp[1]+ychsz*0.3
        xyouts, tx,ty,/normal,alignment=0.5, strupcase(short_names[probe_id]), color=color
    endforeach

    ; Draw box.
    plot, xrange, yrange, /iso, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, $
        ystyle=1, yrange=yrange, ytitle=ytitle, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
        nodata=1, noerase=1, position=tpos

    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    xyouts, tx,ty, /normal, 'a. XY plane'
    
    label_size = 0.8
    foreach time, pos_times, time_id do begin
        msg = time_string(time)
        tx = tpos[0]+xchsz*1
        ty = tpos[1]+ychsz*(0.5+n_elements(pos_times)-time_id-1)
        plots, tx,ty+ychsz*0.2*label_size,/normal, psym=time_psyms[time_id], symsize=label_size*0.5
        xyouts, tx+xchsz*1,ty,/normal, msg, charsize=label_size
    endforeach
    


    if keyword_set(test) then stop
    sgclose

end
