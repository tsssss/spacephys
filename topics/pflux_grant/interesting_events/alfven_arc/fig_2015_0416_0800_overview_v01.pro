;+
; Plot schematics, position, panels.
;-

pro set_circ, fill=fill

    tmp = smkarthm(0,2*!dpi,20,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs, circ_ys, fill=fill
    
end


function fig_2015_0416_0800_overview_v01_cartoon_panel, my_pos, event_info=event_info

    ypans = [2,5,3]
    nypan = n_elements(ypans)
    poss = sgcalcpos(nypan, ypans=ypans, position=my_pos, ypad=1)
    
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    symsize = 0.5
    label_size = 0.8
    hsize = (!d.name eq 'ps')? 40: 6
    thick = (!d.name eq 'ps')? 40: 4
    
    axis_pos = 0.
    axis_color = sgcolor('silver')
    
    left_label_alignment = 0.5
    left_label_xpos = my_pos[0]+xchsz*6
    
    
    
    
;---ground, aurora, and dmsp.
    tpos = poss[*,0]
    xrange = [-1,1]
    yrange = [1000d,0]
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos
    
    ; axis.
    ;plots, axis_pos+[0,0], yrange, color=axis_color
    tmp = convert_coord(0,0, data=0, to_normal=1)
    foreach tx, tmp[0]+[-1,0,1]*xchsz*0.8 do begin
        txs = tx+[0,0]
        tys = my_pos[[3,1]]
        ;plots, txs, tys, normal=1, color=axis_color
        arrow, txs[0],tys[0],txs[1],tys[1], normal=1, color=axis_color, solid=1, hsize=hsize
    endforeach
    tx = tmp[0]
    ty = my_pos[1]+ychsz*3
    msg = 'B'
    xyouts, tx,ty,msg, alignment=0.5, normal=1, color=axis_color
    
        
    ; Ground.
    ground_xs = [-1,1]*0.618*0.5
    ground_ys = [0,0]
    plots, ground_xs, ground_ys, data=1
    sep = 0.06
    len = ychsz*0.3
    txs = make_bins(ground_xs, sep, inner=1)
    foreach tx,txs do begin
        tmp = convert_coord(tx,ground_ys[0], data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        ttx = tx+[0,1]*len
        tty = ty+[0,1]*len
        plots, ttx, tty, normal=1
    endforeach
    tmp = convert_coord(min(ground_xs),ground_ys[0], data=1, to_normal=1)
    tx = left_label_xpos
    ty = tmp[1]-ychsz*0.2
    msg = 'Earth'
    xyouts, tx,ty,msg, normal=1, color=color, alignment=left_label_alignment
    msg = 'H = 0 km'
    xyouts, tx,ty-ychsz*1,msg, normal=1, color=color, alignment=left_label_alignment, charsize=label_size
    
    ; aurora.
    color = sgcolor('medium_aquamarine')
    arc_ys = [100,150]
    arc_xs = [-1,1]*0.1
    polyfill, arc_xs[[0,0,1,1,0]], arc_ys[[0,1,1,0,0]], color=color
    tmp = convert_coord(0,arc_ys[0], data=1, to_normal=1)
    tx = tmp[0]+xchsz*2.5
    ty = tmp[1]-ychsz*0.8
    msg = 'a) Auroral arc'
    xyouts, tx,ty,msg, normal=1, color=color, alignment=0
    msg = 'H ~ '+string(arc_ys[0],format='(I0)')+' km'
    xyouts, tx,ty-ychsz*1, msg, normal=1, alignment=0, color=color, charsize=label_size
    

    
    ; dmsp.
    color = sgcolor('teal')
    dmsp_ys = [0,0]+800
    dmsp_xs = [1,-1.5]*0.15
    ; arrow, dmsp_xs[0], dmsp_ys[0], dmsp_xs[1], dmsp_ys[1], data=1, color=color, solid=1, hsize=hsize
    plots, dmsp_xs, dmsp_ys, data=1, color=color
    tx = 0
    ty = dmsp_ys[0]
    set_circ, fill=1
    plots, tx,ty, psym=8, symsize=symsize, color=color
    tmp = convert_coord(dmsp_xs[1],dmsp_ys[0], data=1, to_normal=1)
    tx = left_label_xpos
    ty = tmp[1]-ychsz*0.3
    msg = 'DMSP-F19'
    xyouts, tx,ty,msg, normal=1, color=color, alignment=left_label_alignment
    msg = 'H = 800 km'
    xyouts, tx,ty-ychsz*1, msg, normal=1, alignment=left_label_alignment, color=color, charsize=label_size
    
    
;---potential drop, e- and O+.
    tpos = poss[*,1]
    xrange = [-1,1]
    yrange = [2.2,1+2000/constant('re')]  ; Re.
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos
    
    ; axis.
    ;plots, axis_pos+[0,0], yrange, color=axis_color
    symsize2 = 1.2
    
    ; electron.
    color = sgcolor('salmon')
    set_circ, fill=0
    xxs = [-0.032,0.031]
    yys = [1.2,1.4]
    rr = xchsz*0.5
    foreach xx, xxs, pt_id do begin
        yy = yys[pt_id]
        plots, xx,yy, data=1, psym=8, color=color, symsize=symsize2
        tmp = convert_coord(xx,yy, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        tr = rr*0.8
        plots, tx+[-1,1]*tr, ty+[0,0], color=color, normal=1
        ;plots, tx+[0,0], ty+[-1,1]*tr, color=color, normal=1
        ;plots, tx+[0,0], ty-[rr,ychsz*2], color=color, normal=1
        arrow, tx,ty+rr, tx,ty+ychsz*2, color=color, normal=1, solid=1, hsize=hsize
    endforeach
    tmp = convert_coord(0,mean(yys), data=1, to_normal=1)
    tx = tmp[0]+xchsz*2.5
    ty = tmp[1]-ychsz*0.8
    msg = 'b) Inverted-V e-'
    xyouts, tx,ty,msg, normal=1, color=color, alignment=0
    msg = '1-10 keV'
    xyouts, tx,ty-ychsz*1, msg, normal=1, color=color, alignment=0, charsize=label_size
    
    ; O+.
    color = sgcolor('blue')
    xxs = [-0.03,0.04]
    yys = [1.85,2.0]
    rr = xchsz*0.5
    foreach xx, xxs, pt_id do begin
        yy = yys[pt_id]
        plots, xx,yy, data=1, psym=8, color=color, symsize=symsize2
        tmp = convert_coord(xx,yy, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        tr = rr*0.8
        plots, tx+[-1,1]*tr, ty+[0,0], color=color, normal=1
        plots, tx+[0,0], ty+[-1,1]*tr, color=color, normal=1
        ;plots, tx+[0,0], ty-[rr,ychsz*2], color=color, normal=1
        arrow, tx,ty-rr, tx,ty-ychsz*2, color=color, normal=1, solid=1, hsize=hsize
    endforeach
    tmp = convert_coord(0,mean(yys), data=1, to_normal=1)
    tx = tmp[0]+xchsz*2.5
    ty = tmp[1]-ychsz*0.8
    msg = 'c) O+ outflow'
    xyouts, tx,ty,msg, normal=1, color=color, alignment=0
    msg = '1-10 keV'
    xyouts, tx,ty-ychsz*1, msg, normal=1, color=color, alignment=0, charsize=label_size

    
    ; potential drop
    pd_range = [1.5,1.7]
    color = sgcolor('gray')
    pd_ys = smkarthm(pd_range[0],pd_range[1],0.1, 'dx')
    pd_x = 0
    pd_xrange = [-1,1]*0.2
    foreach ty, pd_ys do begin
        txs = smkarthm(pd_xrange[0],pd_xrange[1],40,'n')
        tys = (txs/pd_xrange[0]*1.2)^4+ty
        index = where(tys le 2.5)
        plots, txs[index],tys[index], color=color
    endforeach
    
    text_xpos = -pd_x-0.15
    plots, text_xpos+[0,0], pd_range
    foreach ty, pd_range do begin
        tmp = convert_coord(text_xpos,ty, data=1, to_normal=1)
        txs = tmp[0]+[-1,1]*xchsz*0.15
        tys = tmp[1]+[0,0]
        plots, txs,tys, normal=1
    endforeach
    ty = mean(pd_range)
    tmp = convert_coord(text_xpos,ty, data=1, to_normal=1)
    tx = left_label_xpos;+xchsz*2
    ty = tmp[1]-ychsz*0.3
    msg = 'Potential '+tex2str('Phi')+'!D||'
    xyouts, tx,ty,msg, normal=1, alignment=left_label_alignment
    msgs = ['1-10 kV','H ~ 0.7 Re','1 Re = '+string(constant('re'),format='(I0)')+' km']
    foreach msg, msgs, ii do begin
        xyouts, tx,ty-ychsz*(ii+1), msg, normal=1, alignment=left_label_alignment, charsize=label_size
    endforeach
    
    
    
;---rbsp, and S.
    tpos = poss[*,2]
    
    xrange = [-1,1]
    yrange = [6,5]  ; Re.
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos

    ; axis.
    ;plots, axis_pos+[0,0], yrange, color=axis_color
    
    ; equator.
    color = sgcolor('black')
    equator_xs = [-1,1]*0.618*0.5
    equator_ys = yrange[0]+[0,0]
    plots, equator_xs, equator_ys, data=1, color=color
    tmp = convert_coord(min(equator_xs),equator_ys[0], data=1, to_normal=1)
    tx = left_label_xpos
    ty = tmp[1]+ychsz*2
    msg = 'Magnetic!Cequator'
    xyouts, tx,ty,msg, normal=1, color=color, alignment=left_label_alignment
    msg = 'H ~ 7.1 Re'
    xyouts, tx,ty-ychsz*2, msg, normal=1, color=color, alignment=left_label_alignment, charsize=label_size
    
    ; rbsp.
    color = sgcolor('magenta')
    rbsp_ys = yrange[1]-total(yrange*[-1,1])*0.1
    rbsp_xs = 0
    rbsp_x = rbsp_xs
    rbsp_y = rbsp_ys
    ; use panel position.
    ty = (poss[1,1]+poss[3,2])*0.5
    tmp = convert_coord(0,ty, normal=1, to_data=1)
    rbsp_x = 0
    rbsp_y = tmp[1]
    plots, rbsp_x+[-1,1]*0.25, rbsp_y+[0,0], data=1, color=color
    set_circ, fill=1
    plots, rbsp_x,rbsp_y, psym=8, symsize=symsize, color=color
    tmp = convert_coord(rbsp_x,rbsp_y, data=1, to_normal=1)
    tx = left_label_xpos;+xchsz*2
    ty = tmp[1]+ychsz*0.0
    msg = 'RBSP-A & B'
    xyouts, tx,ty,msg, normal=1, color=color, alignment=left_label_alignment
    msg = 'H ~ 6.3 Re'
    xyouts, tx,ty-ychsz*1, msg, normal=1, color=color, alignment=left_label_alignment, charsize=label_size
    
    ; alfven wave.
    color = sgcolor('red')
    wave_range = total(yrange*[-1,1])*0.6
    line_range = yrange[1]-total(yrange*[-1,1])*0.2
    yr = yrange[0]+[0,wave_range]
    wave_length = wave_range*1/1.5
    tys = smkarthm(yr[0],yr[1],30,'n')
    foreach xshift,[-1,1]*0.1 do begin
        txs = sin(tys/wave_length*2*!dpi)*0.03+xshift
        plots, txs,tys, data=1, color=color
        tx = txs[-1]
        ty = tys[-1]
        arrow, tx,ty,tx,line_range, data=1, solid=1, hsize=hsize, color=color
    endforeach
    tx = 0.1
    ty = tys[-1]
    tmp = convert_coord(0,ty, data=1, to_normal=1)
    tx = tmp[0]+xchsz*2.5
    ty = tmp[1]-ychsz*1.5
    msg = 'd) Alfven wave'
    xyouts, tx,ty,normal=1, msg, color=color
    
    return, 1
end


function fig_2015_0416_0800_overview_v01_config_panel, my_pos, event_info=event_info, $
    pos_xrange=pos_xrange, pos_yrange=pos_yrange, config_time=config_time, _extra=ex

;---SC position.
    rbsp_info = event_info.rbsp.rbspa
    probe = rbsp_info['probe']
    prefix = rbsp_info['prefix']
    tpos = my_pos
    fig_label = ' '
    label_size = 0.8

    the_time = config_time
    fline_thick = (keyword_set(test))? 2: 4
    fline_colors = sgcolor(['red','blue'])
    
    model_setting = rbsp_info['model_setting']
    external_model = event_info['external_model']
    igrf = model_setting['igrf']
    t89_par = model_setting['t89_par']
    dir = 1
    xrange = pos_xrange
    yrange = pos_yrange

    line_mlats = [60,62,64,66,68,70]
    nline = n_elements(line_mlats)
    lines = list()
    ps = geopack_recalc(the_time)
    if n_elements(h0) eq 0 then h0 = 100d
    r0 = h0/constant('re')+1
    rad = constant('rad')
    deg = constant('deg')
    
    r_gsm_var = prefix+'r_gsm'
    r_gsm = transpose(get_var_data(r_gsm_var, at=the_time))
    r_sm = cotran(r_gsm, the_time, 'gsm2sm')
    the_mlt = pseudo_mlt(r_sm)
    line_mlts = the_mlt+dblarr(nline)

    par_var = geopack_read_par(the_time+[-1,1]*600, model=external_model, t89_par=t89_par)
    par = get_var_data(par_var, at=the_time)

    model_info = geopack_resolve_model(external_model)
    t89 = model_info.t89
    t96 = model_info.t96
    t01 = model_info.t01
    t04s = model_info.ts04
    storm = model_info.storm

    for i=0,nline-1 do begin
        tmp1 = (24-line_mlts[i])*15*rad  ; angle between mlt and midnight.
        tmp2 = line_mlats[i]*rad
        v_sm = [-cos(tmp2)*cos(tmp1),cos(tmp2)*sin(tmp1),sin(tmp2)]
        v0 = cotran(v_sm, the_time, 'sm2gsm')
        geopack_trace, v0[0],v0[1],v0[2], dir, par, xf,yf,zf, fline=fline, $
            t89=t89, t96=t96, t01=t01, ts04=t04s, storm=storm, igrf=igrf, r0=r0
        ntime = n_elements(fline)/3
        fline = cotran(fline, fltarr(ntime)+the_time, 'gsm2sm')
        srotate, fline, (24-the_mlt)*15*rad, 2
        lines.add, fline
    endfor

    ; settings.
    step = 5
    if n_elements(xrange) eq 0 then begin
        xr = []
        foreach line, lines do xr = [xr, minmax(lines[*,0])]
        xr = minmax(make_bins(xr,step))
    endif else xr = xrange
    if n_elements(yrange) eq 0 then begin
        yr = []
        foreach line, lines do yr = [yr, minmax(lines[*,2])]
        yr = minmax(make_bins(yr,step))
    endif else yr = yrange

    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = ychsz*0.15*fig_size[1]
    
    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    xminor = 2
    xtickv = make_bins(xr, xminor, inner=1)
    xticks = n_elements(xtickv)
    xtickn = string(xtickv,format='(I0)')
    ;    xtitle = 'SM R!DXY!N (Re)'
    xtickn[-1] = 'SM R!DXY!N (Re)        '
    xtitle = ' '
    yminor = 2
    ytickv = make_bins(yr, yminor, inner=1)
    yticks = n_elements(ytickv)
    ytitle = 'SM Z (Re)'
    
    xstep = 4
    xgrids = make_bins(xr, xstep, inner=1)
    ystep = 2
    ygrids = make_bins(yr, ystep, inner=1)

    ; Add axis.
    plot, xr, yr, $
        xstyle=1, xrange=xr, xtitle=xtitle, xminor=xminor, xticks=xticks, xtickv=xtickv, xtickn=xtickn, $
        ystyle=1, yrange=yr, ytitle=ytitle, yminor=yminor, yticks=yticks, ytickv=ytickv, $
        nodata=1, noerase=1, position=tpos, iso=1, $
        xticklen=xticklen, yticklen=yticklen
    foreach val, xgrids do begin
        oplot, val+[0,0], yr, linestyle=1
    endforeach
    foreach val, ygrids do begin
        oplot, xr, val+[0,0], linestyle=1
    endforeach
    ; Add SC.
    the_r_gsm = get_var_data(r_gsm_var, at=the_time)
    the_r_sm = transpose(cotran(the_r_gsm, the_time, 'gsm2sm'))
    srotate, the_r_sm, (24-the_mlt)*15*rad, 2
    the_sc_pos = the_r_sm

    sc_neighbors = list()
    foreach dir, [-1,1], ii do begin
        geopack_trace, the_r_gsm[0],the_r_gsm[1],the_r_gsm[2], dir, par, xf,yf,zf, fline=fline, $
            t89=t89, t96=t96, t01=t01, ts04=t04s, storm=storm, igrf=igrf, r0=r0
        ntime = n_elements(fline)/3
        fline = cotran(fline, fltarr(ntime)+the_time, 'gsm2sm')
        srotate, fline, (24-the_mlt)*15*rad, 2
        oplot, fline[*,0], fline[*,2], linestyle=0, color=fline_colors[ii], thick=fline_thick
        
        sc_neighbors.add, reform(fline[1,*])
        
        ; Add Fpt.
        if dir eq -1 then begin
            f_gsm = [xf,yf,zf]
            f_mag = cotran(f_gsm, the_time, 'gsm2mag')
            fmlat = asin(f_mag[2]/r0)*deg
            tmp = convert_coord(fline[0,0],fline[0,2], data=1, to_normal=1)
        endif
    endforeach


    bline_color = sgcolor('silver')
    foreach line, lines do begin
        oplot, line[*,0], line[*,2], color=bline_color
    endforeach


    ; Add SC.
    sc_color = rbsp_info['sc_color']
    sc_name = rbsp_info['sc_name']
    probe = rbsp_info['probe']
    sc_p0 = the_sc_pos[[0,2]]
    plots, sc_p0[0], sc_p0[1], psym=6, color=sc_color;, symsize=label_size
    tmp = convert_coord(sc_p0[0], sc_p0[1], data=1, to_normal=1)
    tx = tmp[0]+xchsz*1
    ty = tmp[1]+ychsz*0.2
    msg = sc_name+'-'+strupcase(probe)
    xyouts, tx,ty,normal=1, msg, charsize=sc_label_size, color=sc_color

    ; Add RBSP-B.
    rbspb_info = event_info.rbsp.rbspb
    sc_color = rbspb_info['sc_color']
    sc_name = rbspb_info['sc_name']
    probe = rbspb_info['probe']
    prefix = rbspb_info['prefix']
    the_r_gsm = get_var_data(prefix+'r_gsm', at=the_time)
    the_r_sm = transpose(cotran(the_r_gsm, the_time, 'gsm2sm'))
    the_mlt = pseudo_mlt(the_r_sm)
    srotate, the_r_sm, (24-the_mlt)*15*rad, 2
    the_sc_pos = the_r_sm
    sc_p0 = the_sc_pos[[0,2]]
    plots, sc_p0[0], sc_p0[1], psym=6, color=sc_color;, symsize=label_size
    tmp = convert_coord(sc_p0[0], sc_p0[1], data=1, to_normal=1)
    tx = tmp[0]+xchsz*1
    ty = tmp[1]-ychsz*0.8
    msg = sc_name+'-'+strupcase(probe)
    ;msg = strupcase(probe)
    xyouts, tx,ty,normal=1, msg, charsize=sc_label_size, color=sc_color




    ; Add DMSP.
    dmsp_info = event_info.dmsp.dmspf19
    dmsp_dis = 1+800d/constant('re')
    dmsp_mlat = (180+64)*constant('rad')
    dmsp_color = dmsp_info['sc_color']
    tx = dmsp_dis*cos(dmsp_mlat)
    ty = dmsp_dis*sin(dmsp_mlat)
    plots, tx,ty, psym=8, symsize=0.5, color=dmsp_color
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx0 = tmp[0]
    ty0 = tmp[1]
    tx1 = tx0+xchsz*1.5
    ty1 = ty0+ychsz*0.5
    plots, [tx0,tx1],[ty0,ty1], normal=1, color=dmsp_color
    probe = dmsp_info['probe']
    dmsp_probe = dmsp_info['probe']
    sc_name = dmsp_info['sc_name']
    xyouts, tx1+xchsz*0., ty1+ychsz*0.2, normal=1, strupcase(sc_name+' '+probe), charsize=sc_label_size, color=dmsp_color


    ; Add earth.
    tmp = smkarthm(0,2*!dpi, 40, 'n')
    txs = cos(tmp)
    tys = sin(tmp)
    polyfill, txs>0, tys, color=sgcolor('white')
    polyfill, txs<0, tys, color=sgcolor('grey')
    plots, txs, tys

    ; Add labeling.
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
;    msg = time_string(the_time,tformat='YYYY-MM-DD/hh:mm')+' UT'
;    xyouts, tx,ty,normal=1, msg, charsize=label_size
;    ty = tpos[3]-ychsz*2
    msg = 'Model: '+strupcase(external_model)
    xyouts, tx,ty,normal=1, msg, charsize=label_size

    tx = tpos[0]-xchsz*8
    ty = tpos[3]-ychsz*0.6
    xyouts, tx,ty,fig_label, normal=1
end


function fig_2015_0416_0800_overview_v01_right_panel, my_pos, event_info=event_info



end


function fig_2015_0416_0800_overview_v01, event_info=event_info

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)

;---Plot file.
    test = 1
    test_plot_panels = 0
    if keyword_set(test) then plot_file = 0

;---Figure out panel size.
    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz
    margins = [6,3,10,1]
    panel_margins = [0,0,0,0]
    uniform_ticklen = -abs_ychsz*0.15

    left_ypad = 1
    xpad = 2
    pos_xrange = [2,-10]
    pos_yrange = [-1,1]*3
    pos_aspect_ratio = abs(total(pos_xrange*[-1,1])/total(pos_yrange*[-1,1]))
    pos_xpan_size = 2
    pos_ypan_size = pos_xpan_size/pos_aspect_ratio

    arc_aspect_ratio = 1
    cartoon_aspect_ratio = 1d/1.2
    cartoon_xpan_size = pos_xpan_size
    cartoon_ypan_size = cartoon_xpan_size/cartoon_aspect_ratio

    panel_xpan_size = 1.5*pos_xpan_size

    fig_ysize = cartoon_ypan_size+left_ypad*abs_ychsz+pos_ypan_size+total(margins[[1,3]])*abs_ychsz
    fig_xsize = pos_xpan_size+xpad*abs_xchsz+total(panel_margins[[0,2]])*abs_xchsz+panel_xpan_size+total(margins[[0,2]])*abs_xchsz

    fig_size = [fig_xsize,fig_ysize]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    
    all_poss = sgcalcpos(1,2, xpans=[pos_xpan_size,panel_xpan_size], margins=margins, xpad=xpad)
    left_pos = all_poss[*,0]
    right_pos = all_poss[*,1]
    left_poss = sgcalcpos(2,1, ypans=[cartoon_ypan_size,pos_ypan_size], position=left_pos, ypad=left_ypad)
    cartoon_pos = left_poss[*,0]
    cartoon_pos[0] -= xchsz*(margins[0]-1)
    config_pos = left_poss[*,1]


;---Test panels.
    if keyword_set(test_plot_panels) then begin
        panel_list = list()
        panel_list.add, cartoon_pos
        panel_list.add, config_pos
        panel_list.add, right_pos

        foreach tpos, panel_list do begin
            xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            plot, [0,1],[0,1], $
                xstyle=1, ystyle=1, $
                xtickformat='(A1)', ytickformat='(A1)', $
                xticklen=xticklen, yticklen=yticklen, $
                nodata=1, noerase=1, position=tpos
        endforeach
    endif
    

;---Panels.
    tpos = fig_2015_0416_0800_overview_v01_cartoon_panel(cartoon_pos, event_info=event_info)

    config_time = time_double('2015-04-16/08:03:30')
    tpos = fig_2015_0416_0800_overview_v01_config_panel(config_pos, event_info=event_info, $
        pos_xrange=pos_xrange, pos_yrange=pos_yrange, config_time=config_time)
    
    tpos = fig_2015_0416_0800_overview_v01_right_panel(right_pos, event_info=event_info)


;---Done.
    if keyword_set(test) then stop
    sgclose


end


print, fig_2015_0416_0800_overview_v01(event_info=event_info)
end