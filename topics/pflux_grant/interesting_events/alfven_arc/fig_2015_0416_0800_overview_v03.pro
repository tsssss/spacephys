;+
; Place panels on left and right sides of cartoon.
;-

;+
; Plot schematics, position, panels.
;-

pro set_circ, fill=fill

    tmp = smkarthm(0,2*!dpi,20,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs, circ_ys, fill=fill
    
end


function fig_2015_0416_0800_overview_v03_cartoon_panel, my_pos, event_info=event_info

    ypans = [2,5,3]
    nypan = n_elements(ypans)
    poss = sgcalcpos(nypan, ypans=ypans, position=my_pos, ypad=1)
    
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    symsize = 0.5
    label_size = 0.8
    hsize = (!d.name eq 'X')? 6: 120
    thick = (!d.name eq 'X')? 1: 2
    tab = '      '
    
    axis_pos = 0.
    axis_color = sgcolor('tan')
    
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
        arrow, txs[0],tys[0],txs[1],tys[1], normal=1, color=axis_color, solid=1, hsize=hsize, thick=thick
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
    msg = tab+'H ~ '+string(arc_ys[0],format='(I0)')+' km'
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
    symsize2 = 1.1
    
    ; electron.
    color = sgcolor('salmon')
    set_circ, fill=0
    xxs = [-0.032,0.031]
    yys = [1.2,1.4]
    rr = xchsz*symsize2
    foreach xx, xxs, pt_id do begin
        yy = yys[pt_id]
        plots, xx,yy, data=1, psym=8, color=color, symsize=symsize2
        tmp = convert_coord(xx,yy, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        tr = rr*0.3
        plots, tx+[-1,1]*tr, ty+[0,0], color=color, normal=1
        ;plots, tx+[0,0], ty+[-1,1]*tr, color=color, normal=1
        ;plots, tx+[0,0], ty-[rr,ychsz*2], color=color, normal=1
        arrow, tx,ty+rr, tx,ty+ychsz*2, color=color, normal=1, solid=1, hsize=hsize, thick=thick
    endforeach
    tmp = convert_coord(0,mean(yys), data=1, to_normal=1)
    tx = tmp[0]+xchsz*2.5
    ty = tmp[1]-ychsz*0.8
    msg = 'b) Inverted-V e-'
    xyouts, tx,ty,msg, normal=1, color=color, alignment=0
    msg = tab+'1-10 keV'
    xyouts, tx,ty-ychsz*1, msg, normal=1, color=color, alignment=0, charsize=label_size
    
    ; O+.
    color = sgcolor('blue')
    xxs = [-0.03,0.04]
    yys = [1.85,2.0]
    rr = xchsz*symsize2
    fig_size = double([!d.x_size,!d.y_size])
    foreach xx, xxs, pt_id do begin
        yy = yys[pt_id]
        plots, xx,yy, data=1, psym=8, color=color, symsize=symsize2
        tmp = convert_coord(xx,yy, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        tr = rr*0.3
        plots, tx+[-1,1]*tr, ty+[0,0], color=color, normal=1
        plots, tx+[0,0], ty+[-1,1]*tr*fig_size[0]/fig_size[1], color=color, normal=1
        ;plots, tx+[0,0], ty-[rr,ychsz*2], color=color, normal=1
        arrow, tx,ty-rr, tx,ty-ychsz*2, color=color, normal=1, solid=1, hsize=hsize, thick=thick
    endforeach
    tmp = convert_coord(0,mean(yys), data=1, to_normal=1)
    tx = tmp[0]+xchsz*2.5
    ty = tmp[1]-ychsz*0.8
    msg = 'c) O+ outflow'
    xyouts, tx,ty,msg, normal=1, color=color, alignment=0
    msg = tab+'1-10 keV'
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
    msg = 'H ~ 6.70 Re'
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
    msg = 'H ~ 5.95 Re'
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
        arrow, tx,ty,tx,line_range, data=1, solid=1, hsize=hsize, color=color, thick=thick
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


function fig_2015_0416_0800_overview_v03_config_panel, my_pos, event_info=event_info, $
    pos_xrange=pos_xrange, pos_yrange=pos_yrange, config_time=config_time, _extra=ex

;---SC position.
    rbsp_info = event_info.rbsp.rbspa
    probe = rbsp_info['probe']
    prefix = rbsp_info['prefix']
    tpos = my_pos
    fig_label = ' '
    label_size = 0.8

    the_time = config_time
    fline_thick = (!d.name eq 'X')? 4: 8
    ;fline_colors = sgcolor(['red','blue'])
    fline_colors = sgcolor(['tan','gray'])
    
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
    uniform_ticklen = -ychsz*0.15*fig_size[1]
    
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
    the_fline = []
    foreach dir, [-1,1], ii do begin
        geopack_trace, the_r_gsm[0],the_r_gsm[1],the_r_gsm[2], dir, par, xf,yf,zf, fline=fline, $
            t89=t89, t96=t96, t01=t01, ts04=t04s, storm=storm, igrf=igrf, r0=r0
        ntime = n_elements(fline)/3
        fline = cotran(fline, fltarr(ntime)+the_time, 'gsm2sm')
        srotate, fline, (24-the_mlt)*15*rad, 2
        ;oplot, fline[*,0], fline[*,2], linestyle=0, color=fline_colors[ii], thick=fline_thick
        the_fline = (dir eq -1)? [the_fline,fline]: [reverse(fline,1),the_fline]
        tmp = fline
        drs = snorm(tmp[1:-1,*]-tmp[0:-2,*])
        fline_length = total(drs)+r0-1
        print, 'RBSP-'+strupcase(probe)+' fline (Re): '+string(fline_length,format='(F5.2)')
        
        sc_neighbors.add, reform(fline[1,*])
        
        ; Add Fpt.
        if dir eq -1 then begin
            f_gsm = [xf,yf,zf]
            f_mag = cotran(f_gsm, the_time, 'gsm2mag')
            fmlat = asin(f_mag[2]/r0)*deg
            tmp = convert_coord(fline[0,0],fline[0,2], data=1, to_normal=1)
        endif
    endforeach
    
    ; Plot the south fline.
    index = where(the_fline[*,2] lt 0, complement=index2)
    oplot, the_fline[index,0], the_fline[index,2], color=fline_colors[0], thick=fline_thick
    tmp = the_fline[index,*]
    drs = snorm(tmp[1:-1,*]-tmp[0:-2,*])
    fline_length = total(drs)+r0-1
    print, 'S-hem fline (Re): '+string(fline_length,format='(F5.2)')
    
    oplot, the_fline[index2,0], the_fline[index2,2], color=fline_colors[1], thick=fline_thick
    tmp = the_fline[index2,*]
    drs = snorm(tmp[1:-1,*]-tmp[0:-2,*])
    fline_length = total(drs)+r0-1
    print, 'N-hem fline (Re): '+string(fline_length,format='(F5.2)')
    

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
    ty1 = ty0+ychsz*0.3
    plots, [tx0,tx1],[ty0,ty1], normal=1, color=dmsp_color
    probe = dmsp_info['probe']
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


function fig_2015_0416_0800_overview_v03_dmsp_panel, dmsp_poss, event_info=event_info

    label_size = 0.8
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]

    dmsp_info = event_info.dmsp.dmspf19
    prefix = dmsp_info['prefix']
    probe = dmsp_info['probe']
    dmsp_orbit_time_range = time_double(['2015-04-16/08:01','2015-04-16/08:06'])
    dmsp_plot_time_range = time_double(['2015-04-16/08:00','2015-04-16/08:07'])
    invertedv_times = time_double('2015-04-16/'+['08:02:50','08:03:30'])
    invertedv_times = time_double('2015-04-16/'+['08:02:50','08:03:22','08:03:38'])
    invertedv_color = sgcolor('salmon')
    invertedv_text = 'Inverted-V'
    dmsp_color = dmsp_info.sc_color
    ssusi_id = 'energy'
    ssusi_wavelength = strupcase(ssusi_id)

;---Load data.
    dmsp_mlt_image_var = dmsp_read_mlt_image(dmsp_plot_time_range+[-1800,0], probe=probe, id=ssusi_id)
    mlat_vars = dmsp_read_mlat_vars(dmsp_plot_time_range, probe=probe, errmsg=errmsg)
    ele_spec_var = dmsp_read_en_spec(dmsp_plot_time_range, probe=probe, species='e', errmsg=errmsg)
    ion_spec_var = dmsp_read_en_spec(dmsp_plot_time_range, probe=probe, species='p', errmsg=errmsg)
    ele_eflux_var = dmsp_read_eflux(dmsp_plot_time_range, probe=probe, species='e', errmsg=errmsg)
    db_xyz_var = dmsp_read_bfield_madrigal(dmsp_plot_time_range, probe=probe)
    r_var = dmsp_read_orbit(dmsp_plot_time_range, probe=probe)

    mlat_var = mlat_vars[0]
    mlt_var = mlat_vars[1]
    mlon_var = mlat_vars[2]


    ; ssusi eflux along sc track.
    min_mlat = 50.
    mlats = get_var_data(mlat_var, in=dmsp_orbit_time_range, times=the_times)
    mlons = get_var_data(mlon_var, at=the_times)
    ; use aacgm mlon.
    mlts = aacgm_mlon2mlt(mlons, the_times)
    sc_tt = (mlts*15-90)*constant('rad')
    sc_rr = (90-abs(mlats))/(90-min_mlat)
    sc_xx = sc_rr*cos(sc_tt)
    sc_yy = sc_rr*sin(sc_tt)
    ntime = n_elements(the_times)
    ssusi_eflux = fltarr(ntime)
    
    ; mlt image.
    get_data, dmsp_mlt_image_var, times, mlt_images, limits=lim
    tmp = min(times-mean(dmsp_plot_time_range), abs=1, time_id)
    mlt_image = reform(mlt_images[time_id,*,*])
    pixel_mlat = lim.pixel_mlat
    pixel_mlt = lim.pixel_mlt
    
;    ; use ssusi time per pixel.
;    files = dmsp_load_ssusi(dmsp_plot_time_range+[-120d,0]*60, probe=probe, errmsg=errmsg)
;    file = files[0]
;    time_var = 'UT_S'
;    date = dmsp_plot_time_range[0]
;    date = date-(date mod 86400d)
;    pixel_time = netcdf_read_var(time_var, filename=files[0])*3600+date
;    pixel_time = pixel_time[*]
;    mlt_image = mlt_image[*]
;    for ii=0,ntime-1 do begin
;        tmp = min(pixel_time-the_times[ii], abs=1, index)
;        ssusi_eflux[ii] = mlt_image[index]
;    endfor
;    stop


    ; use pixel location.
    min_mlat = 50
    pixel_tt = (pixel_mlt*15-90)*constant('rad')
    pixel_rr = (90-pixel_mlat)/(90-min_mlat)
    pixel_xx = pixel_rr*cos(pixel_tt)
    pixel_yy = pixel_rr*sin(pixel_tt)

    mlt_image = mlt_image[*]
    for ii=0,ntime-1 do begin
        dis = sqrt((pixel_xx-sc_xx[ii])^2+(pixel_yy-sc_yy[ii])^2)
        tmp = min(dis[*], abs=1, index)
        ssusi_eflux[ii] = mlt_image[index]
    endfor
    ssusi_eflux_var = prefix+'ssusi_eflux'
    index = where(ssusi_eflux[1:-1]-ssusi_eflux[0:-2] ne 0)
    store_data, ssusi_eflux_var, the_times[index+1], ssusi_eflux[index+1]
    dmsp_eflux_yrange = [1e-1,1e2]
    add_setting, ssusi_eflux_var, smart=1, dictionary($
        'ylog', 1, $
        'yrange', dmsp_eflux_yrange, $
        'display_type', 'scalar', $
        'short_name', 'Aurora eflux', $
        'unit', 'mW/m!U2!N' )


    ; map eflux and convert unit.
    get_data, prefix+'e_eflux', times, eflux, limits=lim
    var = prefix+'e_eflux_map'
    cmap = 1.4  ; this is for 800 km.
    theta_loss_cones = [45,60]
    ndim = n_elements(theta_loss_cones)
    ntime = n_elements(times)
    eflux_map = fltarr(ntime,ndim+1)
    foreach theta_loss_cone, theta_loss_cones, lc_id do begin
        sr = !dpi*sin(theta_loss_cone*constant('rad'))^2
        eflux_map[*,lc_id+1] = eflux*cmap*sr*1.6e-19*1e4*1e3  ; convert from eV/cm^2-s-sr to mW/m^2.
    endforeach
    store_data, var, times, eflux_map
    yrange = [2e-1,9e2]
    constant = 10d^[0,1,2]
    add_setting, var, smart=1, dictionary($
        'ylog', 1, $
        'yrange', yrange, $
        'constant', constant, $
        'display_type', 'stack', $
        'labels', ['Aurora',tex2str('theta')+'!DLC!N='+string(theta_loss_cones,format='(I0)')+'!U'+tex2str('circ')], $
        'colors', sgcolor(['red','green','blue']), $
        'ytitle', '(mW/m!U2!N)' )
    

    ; convert dB from xyz to fac.
    db_xyz = get_var_data(db_xyz_var, times=times)
    db_fac = db_xyz
    fmlt = get_var_data(mlt_var, at=times)
    fmlat = get_var_data(mlat_var, at=times)
    theta = (fmlt*15-90)*constant('rad')
    n_hat = -[[cos(theta)],[sin(theta)]]
    w_hat = [[cos(theta-0.5*!dpi)],[sin(theta-0.5*!dpi)]]
    db_fac[*,0] = n_hat[*,0]*db_xyz[*,0]+n_hat[*,1]*db_xyz[*,1]
    db_fac[*,1] = w_hat[*,0]*db_xyz[*,0]+w_hat[*,1]*db_xyz[*,1]
    db_var = prefix+'db_fac'
    store_data, db_var, times, db_fac
    add_setting, db_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'dB', $
        'unit', 'nT', $
        'coord', '', $
        'coord_labels', fac_labels )
    var = db_var
    yrange = [-1,1]*300
    ytickv = [-1,0,1]*200
    yticks = n_elements(ytickv)-1
    yminor = 2
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor


;---Plot
    prefix = dmsp_info['prefix']
    dmsp_vars = prefix+['e_en_spec','e_eflux_map']
    dmsp_labels = 'b-'+['1) e-','2) KEflux']
    ndmsp_var = n_elements(dmsp_vars)
    
    
    var = prefix+'e_en_spec'
    get_data, var, times, data, vals
    index = where(data eq 0 or finite(data,nan=1), count)
    if count ne 0 then begin
        data[index] = 0.01
        store_data, var, times, data, vals
    endif
    
    ytickv = 10d^[2,3,4]
    ytickn = '10!U'+string([2,3,4],format='(I0)')
    zrange = [1e5,5e8]
    ztickv = 10d^[5,6,7,8]
    zticks = n_elements(ztickv)-1
    zminor = 9
    ztickn = '10!U'+string([5,6,7,8],format='(I0)')
    ztickn[0:*:2] = ' '
    zticklen = uniform_ticklen/xchsz*1/fig_size[0]
    options, var, 'constant', ytickv
    options, var, 'ytickv', ytickv
    options, var, 'yticks', n_elements(ytickv)-1
    options, var, 'yminor', 10
    options, var, 'ytickname', ytickn
    options, var, 'ytitle', 'Energy!C(eV)'
    options, var, 'zcharsize', label_size
    options, var, 'zrange', zrange
    options, var, 'ztickname', ztickn
    options, var, 'ztickv', ztickv
    options, var, 'zticks', zticks
    options, var, 'zminor', zminor
    options, var, 'zticklen', zticklen
    options, var, 'color_table', 65
    
    the_poss = dmsp_poss
    for pid=0,ndmsp_var-1 do begin
        tpos = the_poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, dmsp_vars[pid], 'xticklen', xticklen
        options, dmsp_vars[pid], 'yticklen', yticklen
    endfor
    
    tplot_options, 'version', 2
    tplot, dmsp_vars, trange=dmsp_plot_time_range, position=the_poss, noerase=1, single_line_uttick=1
    timebar, invertedv_times, linestyle=1, color=invertedv_color
    
    ; Add inverted V labels.
    pid = where(dmsp_vars eq prefix+'e_en_spec')
    tpos = the_poss[*,pid]
    plot, dmsp_plot_time_range, [0,1], $
        xstyle=5, ystyle=5, position=tpos, nodata=1, noerase=1
    foreach invertedv_time, invertedv_times do begin
        tmp = convert_coord(invertedv_time,1, data=1, to_normal=1)
        tx = tmp[0]
        plots, tx,tpos[3],normal=1, psym=8, symsize=0.5, color=invertedv_color
    endforeach
    invertedv_time = mean(invertedv_times)
    tmp = convert_coord(invertedv_time,1, data=1, to_normal=1)
    tx = tmp[0]
    ty = tpos[3]+ychsz*0.4
    msg = invertedv_text
    xyouts, tx,ty,normal=1, msg, charsize=sc_label_size, alignment=0.5, color=invertedv_color
    
    for pid=0,ndmsp_var-1 do begin
        tpos = the_poss[*,pid]
        tx = tpos[0]-xchsz*8
        ty = tpos[3]-ychsz*0.5
        msg = dmsp_labels[pid]
        xyouts, tx,ty,normal=1, msg
    endfor
    
    
    ; Add SSUSI eflux.
    e_eflux_var = prefix+'e_eflux_map'
    invertedv_time_range = time_double('2015-04-16/'+['08:02:20','08:04:00'])
    pid = where(dmsp_vars eq e_eflux_var, count)
    if count ne 0 then begin
        xrange = dmsp_plot_time_range
        yrange = get_setting(e_eflux_var,'yrange')
        ylog = get_setting(e_eflux_var,'ylog', exist)
        if ~exist then ylog=0
        tpos = the_poss[*,pid]
        plot, xrange, yrange, ylog=ylog, xlog=0, $
            xstyle=5, ystyle=5, position=tpos, nodata=1, noerase=1
        ssusi_eflux = get_var_data(prefix+'ssusi_eflux', times=times, in=invertedv_time_range)
        oplot, times, ssusi_eflux, psym=6, symsize=0.25, color=sgcolor('red')
;        foreach time, invertedv_times do oplot, time+[0,0], yrange, color=invertedv_color, linestyle=1
    endif


end


function fig_2015_0416_0800_overview_v03_asi_panel, my_poss, event_info=event_info

    label_size = 0.8
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]

    asi_time = time_double('2015-04-16/08:03:30')
    asi_setting = (event_info.ground)['asi_setting']
    asi_sites = sort_uniq(asi_setting.sites)

    dmsp_info = event_info.dmsp.dmspf19
    prefix = dmsp_info['prefix']
    probe = dmsp_info['probe']
    dmsp_orbit_time_range = time_double(['2015-04-16/08:01','2015-04-16/08:06'])
    dmsp_plot_time_range = time_double(['2015-04-16/08:00','2015-04-16/08:07'])
    invertedv_times = time_double('2015-04-16/'+['08:02:50','08:03:22','08:03:38'])
    invertedv_color = sgcolor('salmon')
    invertedv_text = 'Inverted-V'
    dmsp_color = dmsp_info.sc_color
    ssusi_id = 'energy'
    ssusi_wavelength = strupcase(ssusi_id)

    dmsp_mlt_image_var = dmsp_read_mlt_image(dmsp_plot_time_range+[-1800,0], probe=probe, id=ssusi_id)
    mlat_vars = dmsp_read_mlat_vars(dmsp_plot_time_range, probe=probe, errmsg=errmsg)
    mlat_var = mlat_vars[0]
    mlt_var = mlat_vars[1]
    mlon_var = mlat_vars[2]


    ; Settings.
    mlt_range = [-1d,0]*6
    min_mlat = 50d
    mlat_range = [min_mlat,90]
    ct_ssusi = 70
    ssusi_zrange = [-1,1]*100
    ct_asi = 70
    asi_zlog = 0
    asi_zrange = [-1,1]*800
    color_top = 254
    nangle = 50
    dangle = !dpi/nangle
    angles = make_bins([-!dpi,0], dangle)


    ; xlabels.
    xtickn_pos = (90.-min_mlat)/(90.-min_mlat+6)
    xrange = mlt_range
    xstep = 6d
    xtickv = smkarthm(xrange[0]+1,xrange[1]-1, 2, 'dx')
    xtick_minor = make_bins(xrange, 1, inner=1)
    xtickn = xtickv
    index = where(xtickv lt 0, count)
    if count ne 0 then xtickn[index] += 24
    xtickn = string(xtickn,format='(I02)')
    
    ; ylabels.
    ytick_pos = (-1+1d/12)*!dpi
    ytick_pos = (-1)*!dpi
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange, 10, inner=1)
    ytick_minor = make_bins(yrange, ystep, inner=1)
    ytickn = string(ytickv,format='(I0)')
    
    aurora_info = dictionary($
        'asi', dictionary(), $
        'ssusi', dictionary() )
    
    ; ASI.
    asi_mlt_image_var = 'thg_asf_mlt_image'
    mlt_image = get_var_data(asi_mlt_image_var, at=asi_time)
    npx = n_elements(mlt_image[0,*])
    mlt_image = mlt_image[0:npx*0.5-1,0:npx*0.5-1]
;    mlt_image = mlt_image[*,0:npx*0.5-1]
    if asi_zlog eq 1 then begin
        asi_log_zrange = alog10(asi_zrange)
        asi_zzs = bytscl(alog10(mlt_image), min=asi_log_zrange[0],max=asi_log_zrange[1], top=color_top)
        asi_log_ztickv = make_bins(asi_log_zrange,1, inner=1)
        asi_ztickv = 10^asi_log_ztickv
        asi_zticks = n_elements(asi_ztickv)-1
        asi_zminor = 9
        asi_ztickn = '10!U'+string(asi_log_ztickv,format='(I0)')
    endif else begin
        asi_zzs = bytscl((mlt_image), min=asi_zrange[0],max=asi_zrange[1], top=color_top)
        asi_ztickv = sort_uniq([asi_zrange,0])
        asi_zticks = n_elements(asi_ztickv)-1
        asi_zminor = 10
        asi_ztickn = string(asi_ztickv,format='(I0)')
    endelse
    ;asi_zrange = [0,5e3]
    ;asi_zzs = bytscl((mlt_image), min=asi_zrange[0],max=asi_zrange[1], top=color_top)
    
    mlt_images = get_var_data(dmsp_mlt_image_var, times=times, limits=lim)
    tmp = min(times-mean(dmsp_plot_time_range), abs=1, time_id)
    npx = n_elements(mlt_images[0,0,*])
;    mlt_image = reform(mlt_images[time_id,*,0:npx*0.5-1])
    mlt_image = reform(mlt_images[time_id,0:npx*0.5-1,0:npx*0.5-1])
    npx = n_elements(mlt_image[0,*])
    ssusi_zzs = bytscl(mlt_image, min=ssusi_zrange[0], max=ssusi_zrange[1], top=color_top)
    ssusi_time_range = reform(lim.time_range[time_id,*])
    ssusi_unit = lim.unit

    aurora_info.asi = dictionary($
        'msg', ['a-1) North | white light',strjoin(strupcase(asi_sites),' ')+' | '+$
        time_string(asi_time,tformat='hh:mm:ss')+' UT'], $
        'hemisphere', 'north', $
        'position', my_poss[*,0], $
        'zzs', asi_zzs, $
        'ct', ct_asi )
    dmsp_name = dmsp_info['sc_name']
    aurora_info.ssusi = dictionary($
        'msg', ['a-2) South | '+ssusi_wavelength,strupcase(dmsp_name+' '+probe)+' | '+$
            strjoin(time_string(ssusi_time_range,tformat='hh:mm'),'-')+' UT'], $
        'hemisphere', 'south', $
        'position', my_poss[*,1], $
        'zzs', ssusi_zzs, $
        'ct', ct_ssusi )

    foreach the_info, aurora_info do begin
        tpos = the_info.position

        ; Draw data.
        zzs = the_info.zzs
        ct = the_info.ct
        sgtv, zzs, ct=ct, position=tpos

        ; Add labels, etc.
        ; Draw axes.
        plot, [-1,0], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos
        plots, [0,-1], [0,0], color=sgcolor('silver')
        plots, [0,0], [0,-1], color=sgcolor('silver')
        ;plots, [0,1], [0,0], color=sgcolor('silver')
        

        ; circles for ytickv.
        foreach yminor, ytick_minor, val_id do begin
            rr = (yminor-min_mlat)/(90-min_mlat)
            txs = rr*cos(angles)
            tys = rr*sin(angles)
            linestyle = 1
            index = where(ytickv eq yminor, count)
            if count ne 0 then linestyle = 0
            oplot, txs,tys, linestyle=linestyle, color=sgcolor('silver')
        endforeach


        ; lines for xickv.
        foreach xminor, xtick_minor, val_id do begin
            linestyle = 1
            index = where(xtickv eq xminor, count)
            if count ne 0 then linestyle = 0
            
            tt = (xminor*15-90)*constant('rad')
            txs = [0,1]*cos(tt)
            tys = [0,1]*sin(tt)

            plots, txs,tys, data=1, linestyle=linestyle, color=sgcolor('silver')
        endforeach
        
        ; add yticknames.
        foreach yminor, ytickv, val_id do begin
            rr = 1-(yminor-min_mlat)/(90-min_mlat)
            tt = ytick_pos
            tx = rr*cos(tt)
            ty = rr*sin(tt)
            msg = ytickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*label_size*0.35
            xyouts, tx,ty,normal=1, alignment=0.5, msg, charsize=label_size
        endforeach

        ; add xticknames.
        foreach xminor, xtickv, val_id do begin
            tmp = (xminor*15-90)*constant('rad')
            rr = xtickn_pos
            tx = rr*cos(tmp)
            ty = rr*sin(tmp)
            msg = xtickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            xyouts, tx,ty-ychsz*0.3,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
        endforeach
        
        ; add label of arc.
        label_x = -4
        label_y = 65
        if the_info.hemisphere eq 'south' then label_y = 64
        rr = (90-label_y)/(90-min_mlat)
        tt = (label_x*15-90)*constant('rad')
        tx = rr*cos(tt)
        ty = rr*sin(tt)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        tx1 = tx+xchsz*0.5
        ty1 = ty+ychsz*2
        plots, [tx,tx1],[ty,ty1], normal=1
        plots, tx,ty, psym=8, normal=1, symsize=0.5
        msg = 'Arc '+the_info.hemisphere
        xyouts, tx1+xchsz*0.5, ty1+ychsz*0.3, msg, normal=1, alignment=0.5
        
        ; Add panel label.
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.2
        msgs = the_info.msg
        xyouts, tx,ty,normal=1, msgs[1], charsize=label_size
        xyouts, tx,ty+ychsz*1,normal=1, msgs[0]
    endforeach


    ; Colorbar.
    asi_cbpos = aurora_info.asi.position
    asi_cbpos[1] = asi_cbpos[3]+ychsz*0.4
    asi_cbpos[3] = asi_cbpos[1]+ychsz*0.4
    asi_cbpos[0] += xchsz*1
    asi_cbpos[2] -= xchsz*1
    asi_ztitle = 'N-hem ASI (#)'
    asi_linestyle = 1
    zticklen = uniform_ticklen/(asi_cbpos[3]-asi_cbpos[1])/fig_size[1]
    sgcolorbar, findgen(color_top), horizontal=1, $
        ztitle=asi_ztitle, zrange=asi_zrange, ct=ct_asi, position=asi_cbpos, $
        ztickv=asi_ztickv, zticks=asi_zticks, ztickname=asi_ztickn, zminor=asi_zminor, zticklen=zticklen, log=asi_zlog


    ssusi_cbpos = aurora_info.ssusi.position
    ssusi_cbpos[1] = ssusi_cbpos[3]+ychsz*0.4
    ssusi_cbpos[3] = ssusi_cbpos[1]+ychsz*0.4
    ssusi_cbpos[0] += xchsz*1
    ssusi_cbpos[2] -= xchsz*1
    ztitle = 'S-hem SSUSI '+strupcase(ssusi_id)+' ('+ssusi_unit+')'
    zrange = ssusi_zrange
    sgcolorbar, findgen(color_top), horizontal=1, $
        ztitle=ztitle, zrange=zrange, ct=ct_ssusi, position=ssusi_cbpos, zticklen=zticklen, zminor=5


;---DMSP SSUSI.
    line_color = sgcolor('silver')
    probe = dmsp_info['probe']
    sc_name = dmsp_info['sc_name']
    
    tpos = aurora_info.ssusi.position
    plot, [-1,0], [-1,0], /nodata, /noerase, $
        xstyle=5, ystyle=5, position=tpos
    
    ; Add SC track.
    mlts = get_var_data(mlt_var, in=dmsp_orbit_time_range, times=the_times)
    mlats = get_var_data(mlat_var, in=dmsp_orbit_time_range, times=the_times)
    mlons = get_var_data(mlon_var, at=the_times)
    ; use aacgm mlon.
    mlts = aacgm_mlon2mlt(mlons, the_times)
    
    rad = constant('rad')
    tts = (mlts*15-90)*rad
    rrs = abs((90-abs(mlats))/total(mlat_range*[-1,1]))
    txs = rrs*cos(tts)
    tys = rrs*sin(tts)
    oplot, txs, tys, color=line_color
    
    tx = txs[-1]
    ty = tys[-1]
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]
    msg = strupcase(sc_name+' '+probe)
    xyouts, tx-xchsz*1,ty-ychsz*1.2,normal=1, msg, alignment=0.7, color=dmsp_color

    minor_times = make_bins(dmsp_orbit_time_range, 60, inner=1)
    minor_tts = interpol(tts, the_times, minor_times)
    minor_rrs = interpol(rrs, the_times, minor_times)
    plots, minor_rrs*cos(minor_tts), minor_rrs*sin(minor_tts), psym=8, symsize=0.5, color=line_color

    ; The times.
    major_times = smkarthm(dmsp_orbit_time_range[0],dmsp_orbit_time_range[1], 300, 'dx')
    major_tts = interpol(tts, the_times, major_times)
    major_rrs = interpol(rrs, the_times, major_times)
    major_tickns = time_string(major_times,tformat='hh:mm')
    major_xxs = major_rrs*cos(major_tts)
    major_yys = major_rrs*sin(major_tts)
    foreach msg, major_tickns, ii do begin
        tmp = convert_coord(major_xxs[ii],major_yys[ii], data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]+ychsz*0.4
        plots, tmp[0], tmp[1], normal=1, psym=8, symsize=0.5, color=dmsp_color
        xyouts, tx,ty,normal=1, msg, alignment=0.5, color=dmsp_color, charsize=label_size
    endforeach
    
    
    ; inverted V label.
    foreach invertedv_time, invertedv_times do begin
        invertedv_tt = interpol(tts, the_times, invertedv_time)
        invertedv_rr = interpol(rrs, the_times, invertedv_time)
        tx = invertedv_rr*cos(invertedv_tt)
        ty = invertedv_rr*sin(invertedv_tt)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        plots, tx,ty, normal=1, psym=8, symsize=0.5, color=invertedv_color
        label_mlt = -3.5
        label_mlat = 72
        tr = (90-label_mlat)/(90-min_mlat)
        tt = (label_mlt*15-90)*!dtor
        tx1 = tr*cos(tt)
        ty1 = tr*sin(tt)
        tmp = convert_coord(tx1,ty1, data=1, to_normal=1)
        tx1 = tmp[0]
        ty1 = tmp[1]
        plots, [tx,tx1],[ty,ty1], normal=1, color=invertedv_color
        msg = invertedv_text
        xyouts, tx1+xchsz*0.5, ty1+ychsz*0.3, msg, normal=1, alignment=0.2, color=invertedv_color
    endforeach

;---SC.
    rbsp_info = event_info.rbsp.rbspa
    model_setting = rbsp_info['model_setting']
    internal_model = event_info['internal_model']
    external_model = event_info['external_model']
    ;external_model = 't04s'
    models = model_setting.models
    model_index = where(models eq external_model)
    foreach the_info, aurora_info do begin
        tpos = the_info.position
        plot, [-1,0], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos

        foreach sc_info, event_info.rbsp do begin
            prefix = sc_info['prefix']
            probe = sc_info['probe']
            sc_name = sc_info['sc_name']
            sc_color = sc_info['sc_color']

            suffix = '_'+internal_model+'_'+the_info.hemisphere
            fmlts = get_var_data(prefix+'fmlt'+suffix, at=asi_time)+24
            fmlats = get_var_data(prefix+'fmlat'+suffix, at=asi_time)
            fmlt = fmlts[model_index]
            fmlat = abs(fmlats[model_index])

            tr = (90-fmlat)/(90-min_mlat)
            tt = (fmlt*15-90)*!dtor
            tx = tr*cos(tt)
            ty = tr*sin(tt)
            tmp = convert_coord(tx, ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            plots, tx,ty,normal=1, psym=6, symsize=label_size, color=sc_color
            msg = strupcase(sc_name)+'-'+strupcase(probe)
            if the_info.hemisphere eq 'south' then msg = strupcase(probe)
            xyouts, tx-xchsz*0.5,ty-ychsz*0.5, alignment=1,normal=1, $
                msg, color=sc_color, charsize=sc_label_size
        endforeach
    endforeach

    

end


function fig_2015_0416_0800_overview_v03_rbsp_panel, my_pos, event_info=event_info


    time_range = time_double(['2015-04-16/08:00','2015-04-16/08:25'])
    ;bar_times = !null
    event_info['dp_times'] = time_double('2015-04-16/'+['08:07:42','08:10:00'])
    dp_times = event_info['dp_times']
    bar_times = dp_times[0]+[0,1.5]*60
    ;bar_times = time_double('2015-04-16/'+['08:07:50','08:08:50'])
    rbsp_info = event_info.rbsp.rbspa
    model_setting = rbsp_info['model_setting']
    all_models = model_setting['models']
    internal_model = event_info.internal_model
    external_model = event_info.external_model
    model_index = where(all_models eq external_model)


    ; For O+ tracing.
    obs_info = dictionary()
    break_time = time_double('2015-04-16/08:11:46')
    break_energy = 3200d
    obs_info['o1'] = dictionary($
        'species', 'o', $
        'spec_var', 'rbspa_o_en_spec_para', $
        'nanti_appearance', 1, $
        'npara_appearance', 1, $
        'trs', [time_double('2015-04-16/08:10:38'),break_time], $
        'ens', [9000d,break_energy] )
    obs_info['o2'] = dictionary($
        'species', 'o', $
        'spec_var', 'rbspa_o_en_spec_para', $
        'nanti_appearance', 1, $
        'npara_appearance', 1, $
        'trs', [break_time,time_double('2015-04-16/08:13:41')], $
        'ens', [break_energy,1200d] )
    

    ; Plot settings.
    label_size = 0.8
    tmp = smkarthm(0,2*!dpi,50,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1
    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    

;---ewogram.
    ewo_var = 'thg_asf_ewo'
    mlat_range = [62,64]
    mlt_range = [-4,-2]
    ewo_zrange = [5e2,1e4]

    
    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
    get_data, mlt_image_rect_var, times, data, limits=lim
    ntime = n_elements(times)
    mlt_bins = lim.mlt_bins
    mlt_index = where_pro(mlt_bins, '[]', mlt_range, count=nmlt_bin)
    mlt_bins = mlt_bins[mlt_index]
    index = where(mlt_bins le 0, count)
    if count ne 0 then mlt_bins[index] += 24
    mlat_bins = lim.mlat_bins
    mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nmlat_bin)
    mlat_bins = mlat_bins[mlat_index]
    ewo = total(data[*,mlt_index,mlat_index],3)/nmlat_bin
    store_data, ewo_var, times, ewo, mlt_bins
    yrange = mlt_range+24
    ystep = 1
    ytickv = make_bins(yrange,ystep,inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 4
    add_setting, ewo_var, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'ASI Count', $
        'unit', '#', $
        'ytitle', 'MLT!C(h)', $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'zrange', ewo_zrange, $
        'zlog', 1, $
        'color_table', 49 )
    
;---keogram.
    keo_var = 'thg_asf_keo'
    mlat_range = [60,67]
    mlt_range = [-3.5,-3]
    keo_zrange = ewo_zrange

    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
    get_data, mlt_image_rect_var, times, data, limits=lim
    ntime = n_elements(times)
    mlt_bins = lim.mlt_bins
    mlt_index = where_pro(mlt_bins, '[]', mlt_range, count=nmlt_bin)
    mlt_bins = mlt_bins[mlt_index]
    index = where(mlt_bins le 0, count)
    if count ne 0 then mlt_bins[index] += 24
    mlat_bins = lim.mlat_bins
    mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nmlat_bin)
    mlat_bins = mlat_bins[mlat_index]
    keo = total(data[*,mlt_index,mlat_index],2)/nmlt_bin
    store_data, keo_var, times, keo, mlat_bins
    add_setting, keo_var, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'ASI Count', $
        'unit', '#', $
        'ytitle', 'MLat!C(deg)', $
        'zrange', keo_zrange, $
        'zlog', 1, $
        'color_table', 49 )

    
    ; keogram around sc mlt.
    dmlt = 0.1
    mlat_range = [60,68]

    foreach info, event_info.rbsp do begin
        prefix = info['prefix']
        keo_var = prefix+'keo'
        
        mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
        get_data, mlt_image_rect_var, times, data, limits=lim
        ntime = n_elements(times)
        mlt_bins = lim.mlt_bins
        mlat_bins = lim.mlat_bins
        mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nmlat_bin)
        mlat_bins = mlat_bins[mlat_index]
        data = data[*,*,mlat_index]
        
        
        fmlt_var = prefix+'fmlt_'+internal_model+'_north'
        fmlts = (get_var_data(fmlt_var, at=times))[*,model_index]
        hehe = fltarr(ntime,nmlat_bin)
        foreach time, times, time_id do begin
            the_mlt_range = fmlts[time_id]+[-1,1]*dmlt
            mlt_index = where_pro(mlt_bins, '[]', the_mlt_range, count=nmlt_bin)
            hehe[time_id,*] = total(data[time_id,mlt_index,*],2)/nmlt_bin
        endforeach
        yrange = mlat_range
        ytitle = 'MLat (deg)'

        ystep = 5
        ytickv = make_bins(yrange, ystep, inner=1)
        yticks = n_elements(ytickv)-1
        yminor = ystep
        store_data, keo_var, times, hehe, mlat_bins
        add_setting, keo_var, smart=1, dictionary($
            'display_type', 'spec', $
            'short_name', 'ASI Count', $
            'unit', '#', $
            'ytitle', 'MLat!C(deg)', $
            'zrange', keo_zrange, $
            'zlog', 1, $
            'yrange', yrange, $
            'ytickv', ytickv, $
            'yticks', yticks, $
            'yminor', yminor, $
            'color_table', 49 )
    endforeach
    
    vars = ['thg_asf_'+['ewo','keo'],'rbsp'+['a','b']+'_keo']
    zrange = [2e2,1e4]
    log_zrange = alog10(zrange)
    log_ztickv = make_bins(log_zrange,1,inner=1)
    ztickv = 10d^log_ztickv
    zticks = n_elements(ztickv)-1
    zminor = 9
    ztickn = '10!U'+string(log_ztickv,format='(I0)')
    
    options, vars, 'zrange', zrange
    options, vars, 'ztickv', ztickv
    options, vars, 'zticks', zticks
    options, vars, 'zminor', zminor
    options, vars, 'ztickname', ztickn


;---rbsp particle and B1 spec.
    foreach info, event_info.rbsp do begin
        prefix = info['prefix']
        probe = info['probe']
    
        ; e- spec.
        e_spec_var = prefix+'e_en_spec'+['','_'+['anti','perp','para']]
        ct_electron = 65
        zrange = [1e5,1e10]
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,2,inner=1)
        ztickv = 10.^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        
        var = e_spec_var
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'zminor', zminor
        options, var, 'zcharsize', label_size
        options, var, 'color_table', ct_electron

        yrange = [15,5e4]
        log_yrange = alog10(yrange)
        log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
        log_ytickv = make_bins(log_yrange,1,inner=1)
        ytickv = 10d^log_ytickv
        yticks = n_elements(ytickv)-1
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        yminor = 10
        foreach tx, ytickv, ii do begin
            if tx eq 1 then begin
                ytickn[ii] = '1'
            endif else if tx eq 10 then begin
                ytickn[ii] = '10'
            endif
        endforeach
        options, var, 'yrange', yrange
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'yminor', yminor
        options, var, 'ytickname', ytickn
        options, var, 'constant', ytickv
        options, var, 'ytitle', 'Energy!C(eV)'

        foreach tvar, var do begin
            get_data, tvar, times, data, val, limits=lim
            index = where(finite(data,nan=1) or data eq 0, count)
            if count ne 0 then begin
                data[index] = 0.001
                store_data, tvar, times, data, val, limits=lim
            endif
        endforeach


        ; H+ spec.
        p_spec_var = prefix+'p_en_spec'+['','_'+['anti','perp','para']]
        ct_proton = 63
        zrange_proton = [1e4,1e6]

        zrange = zrange_proton
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,1,inner=1)
        ztickv = 10.^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9

        var = p_spec_var
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'zminor', zminor
        options, var, 'zcharsize', label_size
        options, var, 'color_table', ct_proton

        yrange = [15,5e4]
        log_yrange = alog10(yrange)
        log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
        log_ytickv = make_bins(log_yrange,1,inner=1)
        ytickv = 10d^log_ytickv
        yticks = n_elements(ytickv)-1
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        yminor = 10
        foreach tx, ytickv, ii do begin
            if tx eq 1 then begin
                ytickn[ii] = '1'
            endif else if tx eq 10 then begin
                ytickn[ii] = '10'
            endif
        endforeach
        options, var, 'yrange', yrange
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'yminor', yminor
        options, var, 'ytickname', ytickn
        options, var, 'constant', ytickv
        options, var, 'ytitle', 'Energy!C(eV)'

        foreach tvar, var do begin
            get_data, tvar, times, data, val, limits=lim
            index = where(finite(data,nan=1) or data eq 0, count)
            if count ne 0 then begin
                data[index] = 0.001
                store_data, tvar, times, data, val, limits=lim
            endif
        endforeach


        ; O+ spec.
        o_spec_var = prefix+'o_en_spec'+['','_'+['anti','perp','para']]
        ct_oxygen = 64
        zrange_proton = [1e4,1e6]

        zrange = zrange_proton
        log_zrange = alog10(zrange)
        log_ztickv = make_bins(log_zrange,1,inner=1)
        ztickv = 10.^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        
        var = o_spec_var
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'zminor', zminor
        options, var, 'zcharsize', label_size
        options, var, 'color_table', ct_oxygen

        yrange = [15,5e4]
        log_yrange = alog10(yrange)
        log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
        log_ytickv = make_bins(log_yrange,1,inner=1)
        ytickv = 10d^log_ytickv
        yticks = n_elements(ytickv)-1
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        yminor = 10
        foreach tx, ytickv, ii do begin
            if tx eq 1 then begin
                ytickn[ii] = '1'
            endif else if tx eq 10 then begin
                ytickn[ii] = '10'
            endif
        endforeach
        options, var, 'yrange', yrange
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'yminor', yminor
        options, var, 'ytickname', ytickn
        options, var, 'constant', ytickv
        options, var, 'ytitle', 'Energy!C(eV)'

        foreach tvar, var do begin
            get_data, tvar, times, data, val, limits=lim
            index = where(finite(data,nan=1) or data eq 0, count)
            if count ne 0 then begin
                data[index] = 0.001
                store_data, tvar, times, data, val, limits=lim
            endif
        endforeach
        
        ; density.
        var = prefix+'density_hope'
        options, var, 'yrange', [0.1,1]
        ;options, var, 'labels', 'e- Density!C  >200 eV'
        options, var, 'labels', ' '
        options, var, 'ytitle', '(cm!U-3!N)'
        options, var, 'ystyle', 1

        ; B var.
        var = prefix+'b_gsm'
        get_data, var, times, b_gsm
        b_sm = cotran(b_gsm, times, 'gsm2sm')
        b_sm_var = prefix+'b_sm'
        store_data, b_sm_var, times, b_sm
        add_setting, b_sm_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'B', $
            'unit', 'nT', $
            'coord', 'SM' )
        b_tilt = asin(b_sm[*,2]/snorm(b_sm))*constant('deg')
        b_tilt_var = prefix+'b_tilt'
        store_data, b_tilt_var, times, b_tilt
        add_setting, b_tilt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'arcsin(B!Dz!N/|B|)', $
            'unit', 'deg', $
            'yrange', [0,50], $
            'ytickv', [20,40], $
            'yticks', 1, $
            'yminor', 4, $
            'ystyle', 1 )
        
        ; B spec.
        b_fac_var = prefix+'b1_fac'
        b_vars = stplot_split(b_fac_var)
        the_b_var = b_vars[1]
        b_mor_var = stplot_mor_new(the_b_var, scale_info=scale_info)
        var = b_mor_var
        zrange = [1,1e6]
        log_ztickv = make_bins(alog10(zrange),1,inner=1)
        ztickv = 10d^log_ztickv
        zticks = n_elements(ztickv)-1
        zminor = 9
        ztickn = '10!U'+string(log_ztickv,format='(I0)')
        index = where(log_ztickv eq 0, count)
        if count ne 0 then ztickn[index] = '1'
        index = where(log_ztickv eq 1, count)
        if count ne 0 then ztickn[index] = '10'
        ztickn[0:*:2] = ' '
        
        time_step = info['field_time_step']
        b0_window = info['b0_window']
        yrange = minmax(1d3/[b0_window*0.5,time_step*4])
        log_ytickv = [1,2,3]
        ytickv = 10d^log_ytickv
        yticks = n_elements(ytickv)-1
        yminor = 9
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        index = where(log_ytickv eq 0, count)
        if count ne 0 then ytickn[index] = '1'
        index = where(log_ytickv eq 1, count)
        if count ne 0 then ytickn[index] = '10'
        options, var, 'display_type', 'spec'
        options, var, 'ztitle', 'B!D'+fac_labels[1]+'!N (nT)!U2!N'
        options, var, 'zrange', zrange
        options, var, 'ztickv', ztickv
        options, var, 'zticks', zticks
        options, var, 'zminor', zminor
        options, var, 'ztickname', ztickn
        get_data, var, times, mor, fs, limits=lim
        store_data, var, times, mor, fs*1e3
        options, var, 'yrange', yrange
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'ytickname', ytickn
        options, var, 'yminor', yminor
        options, var, 'ytitle', 'Freq!C(mHz)'
        options, var, 'constant', ytickv
        options, var, 'zcharsize', label_size
    endforeach
    
    tilt_var = 'b_tilt_combo'
    rbsp_colors = list()
    rbsp_probes = list()
    foreach info, event_info.rbsp do begin
        rbsp_colors.add, info['sc_color']
        rbsp_probes.add, info['probe']
    endforeach
    rbsp_colors = rbsp_colors.toarray()
    rbsp_probes = rbsp_probes.toarray()
    index = sort(rbsp_probes)
    rbsp_probes = rbsp_probes[index]
    rbsp_colors = rbsp_colors[index]
    tmp = stplot_merge('rbsp'+rbsp_probes+'_b_tilt', output=tilt_var)
    add_setting, tilt_var, smart=1, dictionary($
        'display_type', 'stack', $
        'ytitle', '(deg)', $
        'yrange', [5,55], $
        'ytickv', [20,40], $
        'yticks', 1, $
        'yminor', 4, $
        'colors', rbsp_colors, $
        'labels', strupcase('rbsp-'+rbsp_probes) )
    ; calculate correlation.
;    ; time to be compared for RBSP-A.
;    the_tr = time_double('2015-04-16/'+['08:05','08:10'])
;    ; time shift to be applied to RBSP-B.
;    dt = 1
;    time_shifts = make_bins(-120+[-1,1]*30, dt)
;    ntime_shift = n_elements(time_shifts)
;    corr = fltarr(ntime_shift)
;    a_tilt = get_var_data('rbspa_b_tilt', in=the_tr, times=times)
;    foreach time_shift, time_shifts, ts_id do begin
;        b_tilt = get_var_data('rbspb_b_tilt', at=times)
;        corr[ts_id] = c_correlate(a_tilt, b_tilt, 0)
;    endforeach
;    max_corr = max(corr, index)
;    best_time_shift = time_shifts[index]
;    print, 'Best time shift (min): ', best_time_shift/60
;    mlt_a = get_var_data(rbsp_read_mlt(time_range, probe='a'), at=the_tr[0])
;    mlt_b = get_var_data(rbsp_read_mlt(time_range, probe='b'), at=the_tr[0]+abs(best_time_shift))
;    omega_b = abs(mlt_a-mlt_b)*15/(best_time_shift/60)
;    print, 'omega (deg/min): ', omega_b
;    stop
    
    
;---E and Poynting flux.
    edot0_fac_var = prefix+'edot0_fac_spinfit_phasef'
    e_vars = stplot_split(edot0_fac_var)
    e_var = e_vars[2]
    var = e_var
    add_setting, var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'E', $
        'unit', 'mV/m' )
    options, var, 'yrange', [-1,1]*120
    options, var, 'ytickv', [-1,1]*100
    options, var, 'yticks', 2
    options, var, 'yminor', 5
    options, var, 'labels', 'E!D'+tex2str('perp')+',out'
        
    
    ; Poynting flux.
    pf_var = 'rbspb_pfdot0_fac_map'
    var = pf_var
    options, var, 'yrange', [-400,100]
    options, var, 'ytickv', [-400,-200,0]
    options, var, 'yticks', 2
    options, var, 'yminor', 4
    options, var, 'constant', [0]
    
    
    pf_var = 'rbspa_pfdot0_fac_spinfit_phasef_map'
    var = pf_var
    options, var, 'yrange', [-60,20]
    options, var, 'ytickv', [-60,-30,0]
    options, var, 'yticks', 2
    options, var, 'yminor', 3
    options, var, 'constant', [0]
    options, var, 'labels', 'S!D'+fac_labels+'!N'

   

;---Set plot_vars.
    big_ypad = 1
    plot_info = orderedhash()
    plot_info['rbspa_o_en_spec_para'] = dictionary($
        'fig_label', 'c-1) S-hem O+' )
    plot_info['rbspa_o_en_spec_anti'] = dictionary($
        'fig_label', 'c-2) N-hem O+', $
        'ypad', big_ypad )    
    plot_info['rbspa_b1_fac_comp2_mor'] = dictionary($
        'fig_label', 'd-1) Wave')
    plot_info['rbspa_pfdot0_fac_spinfit_phasef_map'] = dictionary($
        'fig_label', 'd-2) S FAC')
        ;'ypan', 1, $
        ;'ypad', big_ypad )

    plot_vars = (plot_info.keys()).toarray()
    nvar = n_elements(plot_vars)
    fig_labels = strarr(nvar)
    ypans = fltarr(nvar)
    ypads = fltarr(nvar)
    foreach key, plot_info.keys(), var_id do begin
        info = plot_info[key]
        fig_labels[var_id] = info['fig_label']
        if ~info.haskey('ypan') then info['ypan'] = 1.
        if ~info.haskey('ypad') then info['ypad'] = 0.4
        ypans[var_id] = info['ypan']
        ypads[var_id] = info['ypad']
    endforeach
    ypads = ypads[0:nvar-2]


;---Plot.
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]
    poss = sgcalcpos(nvar, margins=margins, ypans=ypans, ypad=ypads, position=my_pos)
    
    
    ; ticklen.
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, plot_vars[pid], 'xticklen', xticklen
        options, plot_vars[pid], 'yticklen', yticklen
        var = plot_vars[pid]
        display_type = get_setting(var, 'display_type')
        if display_type eq 'spec' then begin
            zticklen = uniform_ticklen/xchsz*1/fig_size[0]
            options, var, 'zticklen', zticklen
        endif
    endfor
    symsize = 0.4
    
    ; tplot.
    tplot_options, 'version', 2
    tplot, plot_vars, position=poss, trange=time_range, noerase=1, single_line_uttick=1
    label_yshift = -ychsz*0.7
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]+label_yshift
        xyouts, tx,ty, normal=1, fig_labels[pid]
    endfor
    if n_elements(bar_times) ne 0 then timebar, bar_times, color=sgcolor('red'), linestyle=1


;---Add labelings.
    ; Add oxygen labels.
    var = 'rbspa_o_en_spec_para'
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.3
        msg = 'Away from Earth in S-hem, PA [0,45] deg'
        xyouts, tx,ty,msg, normal=1, charsize=label_size
    endif
    
    var = 'rbspa_o_en_spec_anti'
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos

        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.3
        msg = 'Toward Earth in S-hem, PA [135,180] deg'
        xyouts, tx,ty,msg, normal=1, charsize=label_size
    endif
    
    
    
    
    ; Add labels for Pflux.
    var = 'rbspa_pfdot0_fac_spinfit_phasef_map'
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=0, $
            nodata=1, noerase=1, position=tpos

        filter = (rbsp_info['pflux_setting']).filter
        spin_period = 11d
        filter[0] = spin_period

        tx = tpos[2]-xchsz*0.5
        ty = tpos[1]+ychsz*0.3
        msg = 'Normalized to 100 km altitude'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size
        ty = tpos[1]+ychsz*1.3
        msg = 'Filtered in '+string(1d3/filter[1],format='(F3.1)')+'mHz-'+string(1d/filter[0],format='(F3.1)')+'Hz'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size

        tmp = convert_coord(xrange[0],0, data=1, to_normal=1)
        tx = tpos[2]-xchsz*0.5
        ty = tmp[1]+ychsz*0.3
        xyouts, tx,ty,normal=1, alignment=1, 'Away from Earth', charsize=label_size;, color=sgcolor('red')
        ty = tmp[1]-ychsz*0.9
        xyouts, tx,ty,normal=1, alignment=1, 'Toward Earth', charsize=label_size;, color=sgcolor('red')
    endif
    
    var = 'rbspb_pfdot0_fac_map'
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        
        xrange = time_range
        yrange = get_setting(var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=0, $
            nodata=1, noerase=1, position=tpos
        
        filter = (rbsp_info['pflux_setting']).filter
        

        tx = tpos[2]-xchsz*0.5
        ty = tpos[1]+ychsz*0.3
        msg = 'Normalized to 100 km altitude'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size
        ty = tpos[1]+ychsz*1.3
        msg = 'Filtered in '+string(1d3/filter[1],format='(F3.1)')+'mHz-'+string(1d/filter[0],format='(I0)')+'Hz'
        xyouts, tx,ty,normal=1, alignment=1, msg, charsize=label_size
        
        tmp = convert_coord(xrange[0],0, data=1, to_normal=1)
        tx = tpos[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.3
        xyouts, tx,ty,normal=1, alignment=0, 'Away from Earth in S-hem', charsize=label_size;, color=sgcolor('red')
        ty = tmp[1]-ychsz*0.9
        xyouts, tx,ty,normal=1, alignment=0, 'Toward Earth in S-hem', charsize=label_size;, color=sgcolor('red')
    endif


;---O+ tracing.
    foreach the_info, obs_info do begin
        spec_var = the_info['spec_var']
        prefix = get_prefix(spec_var)

        get_data, spec_var, times, data, vals, limits=lim
        time_index = where_pro(times, '[]', the_info['trs'], count=ntest_time)
        if ntest_time eq 0 then message, 'Inconsistency ...'

        test_times = times[time_index]
        test_energys = 10.^interpol(alog10(the_info['ens']), the_info['trs'], test_times)
        the_info['test_times'] = test_times
        the_info['test_energys'] = test_energys
        
;        energy_range = minmax(the_info['ens'])
;        test_energys = fltarr(ntest_time)
;        foreach ii, time_index, test_id do begin
;            energy_index = where_pro(reform(vals[ii,*]), '[]', energy_range, count=count)
;            if count eq 0 then message, 'Inconsistency ...'
;            max_flux = max(data[ii,energy_index], index)
;            test_energys[test_id] = vals[ii,energy_index[index]]
;        endforeach

        pid = where(plot_vars eq spec_var)
        tpos = poss[*,pid]
        yrange = lim.yrange
        xrange = time_range
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach test_time, test_times, ii do begin
            plots, test_times[ii], test_energys[ii], psym=1, symsize=symsize, color=sgcolor('red')
        endforeach
    endforeach
    
    
;---Solve for beam generation altitude and time.
    r_gsm_var = prefix+'r_gsm'
    foreach the_info, obs_info do begin
        model_time = (the_info['trs'])[0]
        r_gsm = get_var_data(r_gsm_var, at=model_time)
        par_var = external_model+'_par'
        if check_if_update(par_var, time_range) then par_var = geopack_read_par(time_range, model=external_model)
        par = get_var_data(par_var, at=model_time)

    ;---Get the model field line.
        time = model_time
        model = external_model
        fline_info = dictionary()
        foreach trace_dir, [-1,1] do begin
            msg = (trace_dir eq -1)? 'north': 'south'

            ps = geopack_recalc(time)
            xp = r_gsm[0]
            yp = r_gsm[1]
            zp = r_gsm[2]

            geopack_trace, xp,yp,zp, trace_dir, par, $
                xf,yf,zf, r0=r0, refine=1, ionosphere=1, $
                t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm, igrf=igrf, $
                fline=fline

            ntrace = n_elements(fline[*,0])
            ndim = 3
            bline = fltarr(ntrace,ndim)
            for ii=0,ntrace-1 do begin
                rx = fline[ii,0]
                ry = fline[ii,1]
                rz = fline[ii,2]

                if keyword_set(igrf) then begin
                    geopack_igrf_gsm, rx,ry,rz, bx,by,bz
                endif else begin
                    geopack_dip, rx,ry,rz, bx,by,bz
                endelse

                if model eq 'dip' or model eq 'dipole' or model eq 'igrf' then begin
                    dbx = -1
                    dby = -1
                    dbz = -1
                endif else begin
                    routine = 'geopack_'+model
                    if model eq 't04s' then routine = 'geopack_ts04'
                    call_procedure, routine, par, rx,ry,rz, dbx,dby,dbz
                endelse

                bline[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
            endfor

            fline_info[msg] = dictionary($
                'fline', fline, $
                'bline', bline )
        endforeach
        
    ;---Solve each energy bin.
        species = the_info.species
        case species of
            'e': mass0 = 1d/1836
            'p': mass0 = 1d
            'o': mass0 = 16d
            'he': mass0 = 4d
        endcase
        mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.

        test_energys = the_info.test_energys
        test_times = the_info.test_times
        ntest_energy = n_elements(test_energys)

    ;---Solve the closer hemisphere.
        flines = fline_info.south.fline
        blines = fline_info.south.bline
        bmags = snorm(blines)
        nfline = n_elements(blines[*,0])
        trace_times = dblarr(nfline,ntest_energy)

        foreach energy, test_energys, test_id do begin
            vmag = sqrt(2*energy/mass0)*1e-3
            drs = snorm(flines[1:nfline-1,*]-flines[0:nfline-2,*])*constant('re')
            dts = drs/vmag
            trace_times[0,test_id] = test_times[test_id]
            for ii=1,nfline-1 do trace_times[ii,test_id] = trace_times[ii-1,test_id]-dts[ii-1]
        endforeach

        fline_diss = snorm(flines)
        dis_step = 0.1
        uniform_diss = make_bins(minmax(fline_diss), dis_step)
        nuniform_dis = n_elements(uniform_diss)
        uniform_times = sinterpol(trace_times, fline_diss, uniform_diss)
        time_diffs = dblarr(nuniform_dis)
        for ii=0,nuniform_dis-1 do begin
            time_diffs[ii] = stddev(uniform_times[ii,*])
        endfor
        beam_time_error = min(time_diffs, index)
        beam_dis = uniform_diss[index]
        beam_time = mean(uniform_times[index,*])
        the_info['beam_time'] = beam_time
        the_info['beam_dis'] = beam_dis
        the_info['beam_time_error'] = beam_time_error
        
        ; O+ tracing results.
        ;plots, uniform_times[index,*], test_energys, psym=1, symsize=symsize, color=sgcolor('blue')
        print, 'Beam dis (Re): ', beam_dis
        print, time_string(beam_time)
        print, beam_time_error
        tmp = convert_coord(beam_time, min(test_energys), data=1, to_normal=1)
        tx = tmp[0]
        ty = tpos[1]+ychsz*3
;        msg = 'time of flight'
;        xyouts, tx,ty,msg, normal=1, charsize=label_size, color=sgcolor('blue')
;        stop
        hsize = (!d.name eq 'X')? 6: 120
        thick = (!d.name eq 'X')? 1: 2
        arrow, tx,ty, tx+xchsz*2.5,ty, normal=1, hsize=hsize, solid=1, thick=thick, color=sgcolor('blue')
        msg = 'time of flight'
        xyouts, tx-xchsz*1.5,ty-ychsz*1, msg, normal=1, charsize=label_size, color=sgcolor('blue')


    ;---Basic travel time.
        travel_info = dictionary()

        ; Closer.
        test_times = the_info.test_times
        beam_time = the_info.beam_time
        travel_times = abs(test_times-beam_time)
        travel_info['south'] = travel_times

        ; Farther.
        flines = fline_info.north.fline
        blines = fline_info.north.bline
        bmags = snorm(blines)
        nfline = n_elements(blines[*,0])
        trace_times = dblarr(nfline,ntest_energy)

        foreach energy, test_energys, test_id do begin
            vmag = sqrt(2*energy/mass0)*1e-3
            drs = snorm(flines[1:nfline-1,*]-flines[0:nfline-2,*])*constant('re')
            dts = drs/vmag
            trace_times[0,test_id] = test_times[test_id]
            for ii=1,nfline-1 do trace_times[ii,test_id] = trace_times[ii-1,test_id]-dts[ii-1]
        endforeach

        fline_diss = snorm(flines)
        beam_times = dblarr(ntest_energy)
        for ii=0,ntest_energy-1 do beam_times[ii] = sinterpol(trace_times[*,ii], fline_diss, beam_dis)
        travel_times = abs(test_times-beam_times)
        travel_info['north'] = travel_times
        the_info['travel_info'] = travel_info
    endforeach



;---Trace multiple times.
    south_color = sgcolor('magenta')
    north_color = sgcolor('red')
    coef = 1.0
    L_phi = 0.5
    dphi = 5
    dt_phi = 117d*L_phi/sqrt(dphi)    ; time to go through potential drop one time, t = 117*L/sqrt(phi), L in Re, phi in kV.
    dt_phi = 0

    pid = where(plot_vars eq prefix+'o_en_spec_anti', count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        yrange = lim.yrange
        xrange = time_range
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach the_info, obs_info do begin
            beam_time = the_info.beam_time
            test_energys = the_info.test_energys
            travel_info = the_info.travel_info
            nanti_appearance = the_info['nanti_appearance']
            
            ; North.
            hem_color = north_color
            t0 = beam_time
            ; 1st appearance.
            if nanti_appearance ge 1 then begin            
                t1 = t0+travel_info.north
                plots, t1, test_energys, psym=1, symsize=symsize, color=hem_color
            endif
;            ; 2nd appearnce.
;            if nanti_appearance ge 2 then begin
;                t2 = t1+travel_info.south+dt_phi+$ ; before reaching southern hemisphere.
;                    (travel_info.south+travel_info.north*2+dt_phi*3)/sqrt(coef)  ; after reaching southern hemisphere.
;                plots, t2, test_energys*coef, psym=1, symsize=symsize, color=hem_color
;            endif
            
;            ; South.
;            hem_color = south_color
;            t0 = beam_time
;            ; 1st appearance.
;            t1 = t0+(travel_info.south+travel_info.north*2+dt_phi*2)/sqrt(coef)
;            plots, t1, test_energys*coef, psym=1, symsize=symsize, color=hem_color
        endforeach
    endif


;---Predict H+.
    pid = where(plot_vars eq prefix+'p_en_spec_anti', count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        yrange = lim.yrange
        xrange = time_range
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach the_info, obs_info do begin
            beam_time = the_info.beam_time
            test_energys = the_info.test_energys
            travel_info = the_info.travel_info
            nanti_appearance = the_info['nanti_appearance']
    
            ; North.
            hem_color = north_color
            t0 = beam_time
            ; 1st appearance.
            if nanti_appearance ge 1 then begin
                t1 = t0+(travel_info.north)/4
                plots, t1, test_energys, psym=1, symsize=symsize, color=hem_color
            endif
            ; 2nd appearnce.
            if nanti_appearance ge 2 then begin
                t2 = t1+(travel_info.south+$ ; before reaching southern hemisphere.
                    (travel_info.south+travel_info.north*2+dt_phi*3)/sqrt(coef))/4  ; after reaching southern hemisphere.
                plots, t2, test_energys*coef, psym=1, symsize=symsize, color=hem_color
            endif
            
            ; South.
            hem_color = south_color
            t0 = beam_time
            ; 1st appearance.
            t1 = t0+(travel_info.south+travel_info.north*2+dt_phi*2)/sqrt(coef)
            plots, t1, test_energys*coef, psym=1, symsize=symsize, color=hem_color
        endforeach
    endif
    
    pid = where(plot_vars eq prefix+'p_en_spec_para', count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        yrange = lim.yrange
        xrange = time_range
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        foreach the_info, obs_info do begin
            beam_time = the_info.beam_time
            test_energys = the_info.test_energys
            travel_info = the_info.travel_info
            n_appearance = the_info['npara_appearance']

            ; South.
            hem_color = south_color
            t0 = beam_time
            ; 1st appearance.
            if n_appearance ge 1 then begin
                t1 = t0+(travel_info.south)/4
                plots, t1, test_energys, psym=1, symsize=symsize, color=hem_color
            endif
        endforeach
    endif

end





function fig_2015_0416_0800_overview_v03, event_info=event_info, test=test

;---Load data and settings.
    version = 'v03'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)

;---Plot file.
    test_plot_panels = 0
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,$
        strlowcase('fig_'+event_info.id+'_overview_'+version+'.pdf')])
    if keyword_set(test) then plot_file = 0

;---Figure out panel size.
    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz
    margins = [2,3,6,1]
    panel_margins = [0,0,0,0]
    uniform_ticklen = -abs_ychsz*0.15

    default_ypad = 0.4
    middle_ypad = 1
    xpads = [6,11]

    ; middle column.
    pos_xrange = [2,-10]
    pos_yrange = [-1,1]*3
    pos_aspect_ratio = abs(total(pos_xrange*[-1,1])/total(pos_yrange*[-1,1]))
    pos_xpan_size = 2
    pos_ypan_size = pos_xpan_size/pos_aspect_ratio

    cartoon_aspect_ratio = 1d/1.2
    cartoon_xpan_size = pos_xpan_size
    cartoon_ypan_size = cartoon_xpan_size/cartoon_aspect_ratio

    ; left column.
    left_ypad = 1.5
    asi_xpad = 1.5
    nasi_panel = 2
    asi_cb = 2
    arc_aspect_ratio = 1
    left_ysize = cartoon_ypan_size+middle_ypad*abs_ychsz+pos_ypan_size
    nleft_panel = 2
    left_ypan_size = pos_ypan_size*0.7
    asi_ypan_size = left_ysize-left_ypan_size*nleft_panel-left_ypad*abs_ychsz-default_ypad*(nleft_panel-1)*abs_ychsz-asi_cb*abs_ychsz
    asi_xpan_size = asi_ypan_size*arc_aspect_ratio
    left_xsize = asi_xpan_size*nasi_panel+asi_xpad*abs_xchsz

    ; right column
    right_xpan_size = 1*pos_xpan_size

    fig_ysize = cartoon_ypan_size+middle_ypad*abs_ychsz+pos_ypan_size+total(margins[[1,3]])*abs_ychsz
    fig_xsize = left_xsize+pos_xpan_size+total(xpads)*abs_xchsz+total(panel_margins[[0,2]])*abs_xchsz+right_xpan_size+total(margins[[0,2]])*abs_xchsz

    fig_size = [fig_xsize,fig_ysize]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    
    all_poss = sgcalcpos(1,3, xpans=[left_xsize,pos_xpan_size,right_xpan_size], margins=margins, xpad=xpads)
    left_pos = all_poss[*,0]
    if nleft_panel gt 1 then begin
        left_poss = sgcalcpos(1+nleft_panel, ypans=[asi_ypan_size,left_ypan_size+fltarr(nleft_panel)], ypad=[left_ypad,default_ypad+fltarr(nleft_panel-1)], region=left_pos, margins=[0,0,0,asi_cb])
    endif else begin
        left_poss = sgcalcpos(1+nleft_panel, ypans=[asi_ypan_size,left_ypan_size+fltarr(nleft_panel)], ypad=[left_ypad], region=left_pos, margins=[0,0,0,asi_cb])
    endelse
    asi_poss = sgcalcpos(1,nasi_panel, position=left_poss[*,0], xpad=asi_xpad)
    dmsp_poss = left_poss[*,1:nleft_panel]
    dmsp_margins = [8,0,8,0]
    dmsp_poss[0,*] += xchsz*dmsp_margins[0]
    dmsp_poss[2,*] -= xchsz*dmsp_margins[2]
    middle_pos = all_poss[*,1]
    right_pos = all_poss[*,2]
    middle_poss = sgcalcpos(2,1, ypans=[cartoon_ypan_size,pos_ypan_size], position=middle_pos, ypad=middle_ypad)
    cartoon_pos = middle_poss[*,0]
    cartoon_pos[0] -= xchsz*5
    config_pos = middle_poss[*,1]


;---Test panels.
    if keyword_set(test_plot_panels) then begin
        panel_list = list()
        panel_list.add, cartoon_pos
        panel_list.add, config_pos
        panel_list.add, right_pos
        ;panel_list.add, left_pos
        ;panel_list.add, middle_pos
        for ii=0,nasi_panel-1 do panel_list.add, reform(asi_poss[*,ii])
        for ii=0,nleft_panel-1 do panel_list.add, reform(left_poss[*,ii+1])

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
    tpos = fig_2015_0416_0800_overview_v03_cartoon_panel(cartoon_pos, event_info=event_info)

    config_time = time_double('2015-04-16/08:03:30')
    tpos = fig_2015_0416_0800_overview_v03_config_panel(config_pos, event_info=event_info, $
        pos_xrange=pos_xrange, pos_yrange=pos_yrange, config_time=config_time)
    
    tpos = fig_2015_0416_0800_overview_v03_asi_panel(asi_poss, event_info=event_info)
    tpos = fig_2015_0416_0800_overview_v03_dmsp_panel(dmsp_poss, event_info=event_info)

    tpos = fig_2015_0416_0800_overview_v03_rbsp_panel(right_pos, event_info=event_info)

;---Done.
    if keyword_set(test) then stop
    sgclose
    
    return, plot_file

end

test = 1
print, fig_2015_0416_0800_overview_v03(event_info=event_info, test=test)
end