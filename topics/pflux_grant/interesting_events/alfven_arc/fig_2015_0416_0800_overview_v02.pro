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


function fig_2015_0416_0800_overview_v02_config_panel, my_pos, event_info=event_info, $
    pos_xrange=pos_xrange, pos_yrange=pos_yrange, config_time=config_time, _extra=ex

;---SC position.
    rbsp_info = event_info.rbsp.rbspa
    probe = rbsp_info['probe']
    prefix = rbsp_info['prefix']
    tpos = my_pos
    fig_label = 'c)'
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
        xticklen=xticklen, yticklen=yticklen, xtickformat='(A1)'
    foreach msg, xtickn, id do begin
        tx = xtickv[id]
        tmp = convert_coord(tx,0, data=1, to_normal=1)
        tx = tmp[0]
        ty = tpos[1]-ychsz*1.1
        xyouts, tx,ty,normal=1, msg, alignment=0.5
    endforeach
    
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

    tx = tpos[0]-xchsz*4
    ty = tpos[3]-ychsz*0.6
    xyouts, tx,ty,fig_label, normal=1
end


function fig_2015_0416_0800_overview_v02_asi_panel, my_poss, event_info=event_info

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
        'msg', ['a) North | white light',strjoin(strupcase(asi_sites),' ')+' | '+$
        time_string(asi_time,tformat='hh:mm:ss')+' UT'], $
        'hemisphere', 'north', $
        'position', my_poss[*,0], $
        'zzs', asi_zzs, $
        'ct', ct_asi )
    dmsp_name = dmsp_info['sc_name']
    aurora_info.ssusi = dictionary($
        'msg', ['b) South | '+ssusi_wavelength,strupcase(dmsp_name+' '+probe)+' | '+$
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


function fig_2015_0416_0800_overview_v02_rbsp_panel, my_pos, event_info=event_info


    time_range = time_double(['2015-04-16/08:00','2015-04-16/08:25'])
    ;bar_times = !null
    event_info['dp_times'] = time_double('2015-04-16/'+['08:07:42','08:10:00'])
    dp_times = event_info['dp_times']
    bar_times = dp_times[0]+[0,1.5]*60
    ;bar_times = time_double('2015-04-16/'+['08:07:50','08:08:50'])
    

    ; Plot settings.
    label_size = 0.8
    fig_labels = 'd)'

;---Plot.
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]
    
    tilt_var = 'b_tilt_combo'
    options, tilt_var, 'labels', [' ',' ']
    
    plot_vars = tilt_var
    nvar = n_elements(plot_vars)
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
    tplot, plot_vars, position=poss, trange=time_range, noerase=1, single_line_uttick=1, nouttick=1, novtitle=1
    
    ; manual xticknames.
    tpos = poss[*,-1]
    xstep = 5*60d
    xrange = time_range
    xtickv = make_bins(xrange, xstep, inner=1)
    xtickn = time_string(xtickv, tformat='hhmm')
    xtickn[0] = time_string(xtickv[0], tformat='YYYY MTH DD        ')
    
    plot, xrange, [0,1], $
        xstyle=5, ystyle=5, $
        nodata=1, noerase=1, position=tpos
    foreach msg, xtickn, id do begin
        tx = xtickv[id]
        tmp = convert_coord(tx,0, data=1, to_normal=1)
        tx = tmp[0]
        ty = tpos[1]-ychsz*1.1
        xyouts, tx,ty,normal=1, msg, alignment=0.5
    endforeach
    
    
    label_yshift = -ychsz*0.6
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*4
        ty = tpos[3]+label_yshift
        xyouts, tx,ty, normal=1, fig_labels[pid]
    endfor
;    if n_elements(bar_times) ne 0 then begin
;        get_data, plot_vars, limits=lim
;        foreach dp_time, dp_times, dp_id do begin
;            timebar, dp_time, color=lim.colors[dp_id], linestyle=1
;        endforeach
;    endif
    
    pid = where(plot_vars eq tilt_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        get_data, tilt_var, limits=lim
        plot, time_range, lim.yrange, $
            nodata=1, noerase=1, $
            xstyle=5, ystyle=5, position=tpos
        foreach msg, strupcase('rbsp-'+['a','b']), probe_id do begin
            ty = tpos[1]+ychsz*0.3+(1-probe_id)*ychsz*0.8
            tx = tpos[2]-xchsz*6
            xyouts, tx,ty,msg, normal=1, color=lim.colors[probe_id]
            plots, dp_times[probe_id]+[0,0], lim.yrange, data=1, color=lim.colors[probe_id], linestyle=1
        endforeach
    endif


end





function fig_2015_0416_0800_overview_v02, event_info=event_info, test=test

;---Load data and settings.
    version = 'v02'
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
    margins = [1,1,2,1]
    panel_margins = [0,0,0,0]
    uniform_ticklen = -abs_ychsz*0.15

    ; left column.
    asi_xpad = 1.5
    nasi_panel = 2
    asi_cb = 2
    arc_aspect_ratio = 1
    asi_ypan_size = 1.65
    asi_xpan_size = asi_ypan_size*arc_aspect_ratio
    left_xsize = asi_xpan_size*nasi_panel+asi_xpad*(nasi_panel-1)*abs_xchsz
    left_ysize = asi_ypan_size+asi_cb*abs_ychsz

    ; right column.
    right_ypad = 2
    pos_xrange = [2.,-10]
    pos_yrange = [-1.,1]*3
    pos_aspect_ratio = abs(total(pos_xrange*[-1,1])/total(pos_yrange*[-1,1]))
    pos_xpan_size = 2.
    pos_ypan_size = pos_xpan_size/pos_aspect_ratio
    rbsp_ypan_size = left_ysize-pos_ypan_size-right_ypad*abs_ychsz
    
    xpads = 8
    fig_ysize = left_ysize+total(margins[[1,3]])*abs_ychsz
    fig_xsize = left_xsize+pos_xpan_size+total(xpads)*abs_xchsz+total(margins[[0,2]])*abs_xchsz

    fig_size = [fig_xsize,fig_ysize]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    all_poss = sgcalcpos(1,2, xpans=[left_xsize,pos_xpan_size], margins=margins, xpad=xpads)
    left_pos = all_poss[*,0]
    asi_poss = sgcalcpos(1,nasi_panel, region=left_pos, xpad=asi_xpad, margins=[0,0,0,asi_cb])
    middle_pos = all_poss[*,1]
    middle_pos[[1,3]] += (margins[3]-0.5)*ychsz
    middle_poss = sgcalcpos(2,1, ypans=[pos_ypan_size,rbsp_ypan_size], position=middle_pos, ypad=right_ypad)
    config_pos = middle_poss[*,0]
    rbsp_pos = middle_poss[*,1]
    ;config_pos[[1,3]] += (margins[3]-0.5)*abs_ychsz/fig_ysize


;---Test panels.
    if keyword_set(test_plot_panels) then begin
        panel_list = list()
        panel_list.add, cartoon_pos
        panel_list.add, config_pos
        panel_list.add, rbsp_pos
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
    config_time = time_double('2015-04-16/08:03:30')
    tpos = fig_2015_0416_0800_overview_v02_config_panel(config_pos, event_info=event_info, $
        pos_xrange=pos_xrange, pos_yrange=pos_yrange, config_time=config_time)
    
    tpos = fig_2015_0416_0800_overview_v02_asi_panel(asi_poss, event_info=event_info)
    tpos = fig_2015_0416_0800_overview_v02_rbsp_panel(rbsp_pos, event_info=event_info)

;---Done.    if keyword_set(test) then stop
    sgclose
    
    return, plot_file

end

test = 0
print, fig_2015_0416_0800_overview_v02(event_info=event_info, test=test)
end