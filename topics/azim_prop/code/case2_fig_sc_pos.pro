;+
; Plot SC position in SM in 3D.
;-


test = 0

;---Input.
    probes = ['thd','the','g15','tha','g13','rbspb']
    pos_times = time_double('2014-08-28/'+['10:10','10:20','10:50'])
    ;pos_times = time_double('2014-08-28/'+['10:20'])
    coord = 'sm'


    ct = 40
    top_color = 250
    bottom_color = 100
    nprobe = n_elements(probes)
    index_colors = bytscl(findgen(nprobe), top=top_color-bottom_color)+bottom_color
    probe_colors = fltarr(nprobe)
    foreach color, index_colors, color_id do $
        probe_colors[color_id] = sgcolor(color, ct=ct)

    probe_colors = smkarthm(90,240,nprobe, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)



    time_colors = sgcolor(['red','orange','yellow'])
    time_psyms = [6,7,1]

    short_names = strarr(nprobe)
    project = azim_df_load_project()
    probe_infos = project.probe_infos
    foreach probe, probes, probe_id do begin
        short_name = probe_infos[probe].short_name
        short_names[probe_id] = short_name
    endforeach

    _2014_0828_10_load_data
    pos_coord = list()
    foreach probe, probes do begin
        r_gsm = get_var_data(probe+'_r_gsm', at=pos_times)
        r_coord = cotran(r_gsm, pos_times, 'gsm2'+coord)
        if n_elements(r_coord) eq 3 then r_coord = reform(r_coord,1,3)
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
    ypan_size = abs(xpan_size/total(xrange*[-1,1])*total(yrange*[-1,1]))+$
        abs(xpan_size/total(xrange*[-1,1])*total(zrange*[-1,1]))+abs_ychsz*ypad
    margins = [8,5,2,2]
    fig_xsize = xpan_size+total(margins[[0,2]])*abs_xchsz
    fig_ysize = ypan_size+total(margins[[1,3]])*abs_ychsz


    root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
    ofn = join_path([root_dir,'case2_fig_sc_pos_v1.pdf'])
    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz

    poss = sgcalcpos(2, ypans=[abs(total(zrange*[-1,1])),abs(total(yrange*[-1,1]))], margins=margins, ypad=ypad)


;---Setup coord.
    tpos = poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    plot, xrange, yrange, /iso, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos

    plots, xrange, [0,0], linestyle=1
    plots, [0,0], yrange, linestyle=1

    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    polyfill, circle_x>0, circle_y, color=sgcolor('white')
    plots, circle_x, circle_y
    foreach r, [5,10] do oplot, circle_x*r, circle_y*r, linestyle=1

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
        tmp = min(xx, index)
        tmp = convert_coord(xx[index], yy[index], /data, /to_normal)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.0
        xyouts, tx,ty,/normal, strupcase(short_names[probe_id]), color=color
    endforeach
    
    ; Add time.
    foreach time, pos_times, time_id do begin
        msg = time_string(time,tformat='hh:mm')
        tx = tpos[0]+xchsz*1
        ty = tpos[1]+ychsz*(5+time_id)
        plots, tx,ty+ychsz*0.3,/normal, psym=time_psyms[time_id], symsize=symsize
        xyouts, tx+xchsz*1,ty,/normal, msg
    endforeach

    ; Draw box.
    plot, xrange, yrange, /iso, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, $
        ystyle=1, yrange=yrange, ytitle=ytitle, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
        nodata=1, noerase=1, position=tpos

    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    xyouts, tx,ty, /normal, 'b. XY plane'



;---XZ plane.
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    plot, xrange, zrange, /iso, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=zrange, $
        nodata=1, noerase=1, position=tpos

    ;oplot, xrange, [0,0], linestyle=1
    oplot, [0,0], zrange, linestyle=1

    foreach probe, probes, probe_id do begin
        xx = reform(pos_coord[probe_id,*,0])
        yy = reform(pos_coord[probe_id,*,1])
        zz = reform(pos_coord[probe_id,*,2])
        color = probe_colors[probe_id]
        plots, xx, zz, color=line_color
        foreach time, pos_times, time_id do begin
            ;plots, xx[time_id], yy[time_id], color=color, $
            ;    psym=psym, symsize=symsize
            plots, xx[time_id], zz[time_id], color=color, $
                psym=time_psyms[time_id], symsize=symsize
        endforeach
    endforeach

    ; Draw box.
    plot, xrange, zrange, /iso, $
        xstyle=1, xrange=xrange, xtitle='', xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickformat='(A1)', $
        ystyle=1, yrange=zrange, ytitle=ztitle, ytickv=ztickv, yticks=zticks, yminor=zminor, yticklen=zticklen, $
        nodata=1, noerase=1, position=tpos

    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    xyouts, tx,ty, /normal, 'a. XZ plane'

    if keyword_set(test) then stop
    sgclose

end
