;+
; Show SAPS E field and the association of auroral boundary.
;-

function fig_saps_v02_rbsp_panel, rbsp_pos, event_info=event_info

    label_size = 0.8
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]

    saps_times = orderedhash()
    saps_times['rbspa'] = ['2015-03-17/15:30','2015-03-17/16:00']
    saps_times['rbspb'] = ['2015-03-17/21:05','2015-03-17/21:35']
    
    pp_times = ((event_info['rbsp'])['rbspa'])['pp_times']
    pp_time = pp_times[where_pro(pp_times,'[]',time_double(saps_times['rbspa']))]
    pp_color = sgcolor('red')
    colors = sgcolor(['purple','red'])

    margins = [10,3,0,0]
    poss = sgcalcpos(3, region=rbsp_pos, margins=margins)

;---L shell as x-axis.
    xrange = [2,3.5]
    xstep = 0.3
    xtickv = smkarthm(xrange[0]+0.2,xrange[1],xstep,'dx')
    xticks = n_elements(xtickv)-1
    xminor = 3
    xtitle = 'L-shell (#)'
    xtickn = strarr(xticks+1)+' '
    

;;---Density panel.
;    tpos = poss[*,0]
;    
;    yrange = [30,3e3]
;    ylog = 1
;    ytitle = 'Density (cm!U-3!N)'
;    
;    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
;    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
;    
;    plot, xrange, yrange, $
;        xstyle=5, ystyle=5, $
;        nodata=1, noerase=1, position=tpos, ylog=ylog
;    foreach mission_probe, saps_times.keys(), probe_id do begin
;        probe_info = resolve_probe(mission_probe)
;        prefix = probe_info['prefix']
;        probe = probe_info['probe']
;        lshell_var = prefix+'lshell'
;        e_var = prefix+'edot0_fac'
;        tr = time_double(saps_times[mission_probe])
;        lshell = get_var_data(lshell_var, in=tr, times=times)
;
;        dens_var = rbsp_read_density_emfisis(tr+[-1,1]*600, probe=probe)
;        dens = get_var_data(dens_var, at=times)
;        
;        plots, lshell, dens, color=colors[probe_id]
;        
;;        tx = tpos[0]+xchsz*0.5
;;        ty = tpos[1]+ychsz*(1.5-probe_id)
;;        msg = time_string(tr[0],tformat='YYYY-MM-DD ')+strjoin(time_string(tr,tformat='hh:mm'),'-')+' UT'
;;        xyouts, tx,ty,msg, normal=1, color=colors[probe_id]
;    endforeach
;    
;    plot, xrange, yrange, $
;        xstyle=1, xtitle='', xtickv=xtickv, xminor=xminor, xticks=xticks, $
;        ystyle=1, ytitle=ytitle, ytickv=ytickv, yminor=yminor, yticks=yticks, $
;        nodata=1, noerase=1, position=tpos, $
;        xticklen=xticklen, yticklen=yticklen, xtickn=xtickn, ylog=ylog
;        
;    tx = tpos[0]-xchsz*8
;    ty = tpos[3]-ychsz*0.7
;    msg = 'a)'
;    xyouts, tx,ty,msg, normal=1
;

;---H+ panel.
    
    foreach mission_probe, saps_times.keys(), probe_id do begin
        probe_info = resolve_probe(mission_probe)
        prefix = probe_info['prefix']
        probe = probe_info['probe']
        if probe eq 'b' then continue
        spec_vars = prefix+['p','e','o']+'_en_spec'
        lshell_var = prefix+'lshell'
        tr = time_double(saps_times[mission_probe])
        lshell = get_var_data(lshell_var, in=tr, times=times)
        dlshell = min(lshell[1:-1]-lshell[0:-2])
        the_lshells = make_bins(lshell, dlshell, inner=1)
        the_times = interpol(times,lshell,the_lshells)

        foreach spec_var, spec_vars do begin
            new_var = spec_var+'_lshell'
            copy_data, spec_var, new_var
            interp_time, new_var, the_times
        endforeach

        tposs = poss[*,0:-2]
        vars = prefix+['e','p']+'_en_spec'
        options, vars, 'xstyle', 5
        options, vars, 'ytitle', 'Energy!C(eV)'
        options, vars, 'zticklen', -0.5
        options, vars, 'zminor', 9
        

        the_vars = vars[0] 
        ztickn = '10!U'+['4','5','6','7','8','9','10']
        ztickn[1:*:2] = ' '       
        options, the_vars, ytickname=['10!U2','10!U3','10!U4'], $
            ytickv=[100,1000,10000], yticks=2, yminor=9, $
            ztickv=10d^[4,5,6,7,8,9,10], zticks=6, ztickname=ztickn

        the_vars = vars[1]
        options, the_vars, ytickname=['10!U2','10!U3','10!U4'], $
            ytickv=[100,1000,10000], yticks=2, yminor=9, $
            yrange=[20,3e4]

        
        tplot_options, 'tickinterval', 60
        tplot, vars, trange=tr, position=tposs, noerase=1, nouttick=1
        timebar, pp_time, linestyle=2, color=pp_color

        ; Add label for plasmapause
        tpos = tposs[*,0]
        plot, tr, [0,1], xstyle=5, ystyle=5, nodata=1, noerase=1, position=tpos
        tmp = convert_coord(pp_time,0,data=1,to_normal=1)
        tx = tmp[0]-xchsz*0.5
        ty = tpos[3]-ychsz
        msg = 'Plasmapause'
        xyouts, tx,ty,msg, normal=1, color=pp_color, charsize=label_size, alignment=1
        msg = 'PS'
        tx = tmp[0]+xchsz*0.5
        xyouts, tx,ty,msg, normal=1, color=pp_color, charsize=label_size, alignment=0

        ; Add label for ion nose
        tpos = tposs[*,1]
        plot, tr, [0,1], xstyle=5, ystyle=5, nodata=1, noerase=1, position=tpos
        tmp = convert_coord(pp_time,0,data=1,to_normal=1)
        tx = tmp[0]-xchsz*5
        ty = tpos[3]-ychsz*2
        msg = 'Ion nose'
        xyouts, tx,ty,msg, normal=1, color=sgcolor('yellow'), charsize=label_size, alignment=0.5
        txs = tx+[0,0.7]*xchsz
        tys = ty+[0.7,1.5]*ychsz
        plots, txs, tys, normal=1, color=sgcolor('yellow')
        tmp = smkarthm(0,2*!dpi,20,'n')
        circ_xs = cos(tmp)
        circ_ys = sin(tmp)
        usersym, circ_xs, circ_ys, fill=1
        plots, txs[1],tys[1],normal=1, psym=8, symsize=0.5, color=sgcolor('yellow')


        yrange = [0,1]
        nvar = n_elements(vars)
        fig_labels = letters(nvar)+') '+['e-','H+']+' EN'
        foreach var, vars, pid do begin
            tpos = tposs[*,pid]
            xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            
            plot, xrange, yrange, $
                xstyle=1, xtitle='', xtickv=xtickv, xminor=xminor, xticks=xticks, $
                ystyle=5, ytitle=ytitle, ytickv=ytickv, yminor=yminor, yticks=yticks, $
                nodata=1, noerase=1, position=tpos, $
                xticklen=xticklen, yticklen=yticklen, xtickn=xtickn

            tx = tpos[0]-xchsz*9
            ty = tpos[3]-ychsz*0.7
            msg = fig_labels[pid]
            xyouts, tx,ty,msg, normal=1

            msg = 'RBSP-'+strupcase(probe)
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            xyouts, tx,ty,msg, normal=1, alignment=0, color=sgcolor('white');, charsize=label_size
        endforeach
    endforeach


    

;---E field panel.
    tpos = poss[*,-1]

    yrange = [-2,12]
    ystep = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 5
    ytitle = '(mV/m)'


    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    plot, xrange, yrange, $
        xstyle=5, ystyle=5, $
        nodata=1, noerase=1, position=tpos
    foreach mission_probe, saps_times.keys(), probe_id do begin
        probe_info = resolve_probe(mission_probe)
        prefix = probe_info['prefix']
        probe = probe_info['probe']
        lshell_var = prefix+'lshell'
        e_var = prefix+'edot0_fac'
        tr = time_double(saps_times[mission_probe])
        lshell = get_var_data(lshell_var, in=tr, times=times)
        ilat = lets_calc_ilat(lshell)
        mlt_var = prefix+'mlt'
        mlt = get_var_data(mlt_var, at=times)
        
        e_var_tmp = e_var+'_tmp'
        copy_data, e_var, e_var_tmp
        interp_time, e_var_tmp, times
        efac = get_var_data(e_var_tmp)

        ;plots, ilat, efac[*,2], color=colors[probe_id]
        plots, lshell, efac[*,2], color=colors[probe_id]

        ;time_ticks = time_string(interpol(times, lshell, xtickv), tformat='hh:mm')
        ;xtickn = xtickn+'!C'+time_ticks

        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*(1+probe_id)
        msg = time_string(tr[0],tformat='YYYY-MM-DD ')+strjoin(time_string(tr,tformat='hh:mm'),'-')+' UT'
        xyouts, tx,ty,msg, normal=1, color=colors[probe_id], charsize=label_size
    endforeach


    plot, xrange, yrange, $
        xstyle=1, xtitle='', xtickv=xtickv, xminor=xminor, xticks=xticks, $
        ystyle=1, ytitle=ytitle, ytickv=ytickv, yminor=yminor, yticks=yticks, $
        nodata=1, noerase=1, position=tpos, $
        xticklen=xticklen, yticklen=yticklen, xtickn=xtickn
    
    tx = tpos[0]-xchsz*9
    ty = tpos[3]-ychsz*0.7
    msg = 'c) E!D'+tex2str('perp')+',out'
    xyouts, tx,ty,msg, normal=1
    
    
    ; Add tick names manually.
    foreach tx, xtickv, xid do begin
        ty0 = tpos[1]-ychsz*1.3
        tx0 = tpos[0]

        tmp = convert_coord(tx,yrange[0], data=1, to_normal=1)
        tx = tmp[0]
        ty = ty0

        msg = string(xtickv[xid],format='(F3.1)')
        xyouts, tx,ty, msg, normal=1, alignment=0.5

        if xid eq 0 then begin
            msg = 'L-Shell (#)'
            xyouts, tx0+xchsz*0, ty, msg, normal=1, alignment=1, color=color
        endif

        foreach mission_probe, saps_times.keys(), probe_id do begin
            color = colors[probe_id]
            probe_info = resolve_probe(mission_probe)
            prefix = probe_info['prefix']
            probe = probe_info['probe']
            ty = ty0-ychsz*(probe_id+1)
            lshell_var = prefix+'lshell'
            tr = time_double(saps_times[mission_probe])
            lshell = get_var_data(lshell_var, in=tr, times=times)
            mlt_var = prefix+'mlt'
            the_time = interpol(times, lshell, xtickv[xid])
            ;msg = time_string(the_time,tformat='hh:mm')
            msg = string(get_var_data(mlt_var, at=the_time)+24, format='(F4.1)')
            xyouts, tx,ty, msg, normal=1, alignment=0.5, color=color

            if xid eq 0 then begin
                msg = strupcase('rbsp-'+probe)+' MLT (h)'
                xyouts, tx0+xchsz*0, ty, msg, normal=1, alignment=1, color=color
            endif
        endforeach
    endforeach
    
    return, poss

end

function fig_saps_v02_ssusi_panel, my_poss, event_info=event_info

    label_size = 0.8
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]
    

    dmsp_time_list = list()
    dmsp_time_list.add, ['2015-03-17/16:00','2015-03-17/16:40']
    dmsp_time_list.add, ['2015-03-17/15:15','2015-03-17/15:45']
    npanel = n_elements(dmsp_time_list)

    sc_name = 'dmsp'
    dmsp_color = sgcolor('red')
    probe = 'f18'
    prefix = 'dmsp'+probe+'_'
    mlt_image_var = prefix+'mlt_image'
    mlt_images = get_var_data(mlt_image_var, in=time_range, times=times, limits=lim)

    tmp = smkarthm(0,2*!dpi,50,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1


;---SSUSI.
    color_top = 254
    ct = 70
    ssusi_zlog = 0
    ssusi_zrange = [-1,1]*2.5
    ssusi_min_mlat = 50d
    ssusi_id = 'lhbs'
    ssusi_unit = lim.unit
    
    mlt_range = [-1d,1]*12
    mlt_range = [-12,0]
    min_mlat = 45d
    mlat_range = [min_mlat,90]

    nangle = 50
    dangle = !dpi/nangle
    angles = make_bins(mlt_range*12, dangle)

    ; xlabels.
    xtickn_pos = (90.-min_mlat)/(90.-min_mlat+2)
    xrange = mlt_range
    xstep = 6d
    xtickv = smkarthm(xrange[0],xrange[1]-1, 2, 'dx')
    xtick_minor = make_bins(xrange, 1, inner=1)
    xtickn = xtickv
    index = where(xtickv lt 0, count)
    if count ne 0 then xtickn[index] += 24
    xtickn = string(xtickn,format='(I02)')
    
    ; ylabels.
    ytick_pos = (-1)*!dpi
    ytick_pos = (0.75)*!dpi
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange, 10, inner=1)
    ytick_minor = make_bins(yrange, ystep, inner=1)
    ytickn = string(ytickv,format='(I0)')
    
    for pid=0,npanel-1 do begin
        tpos = my_poss[*,pid]
        ssusi_tpos = tpos
        ratio = abs((90-ssusi_min_mlat)/(90-min_mlat))
        ssusi_tpos[[0,2]] = mean(ssusi_tpos[[0,2]])+total(ssusi_tpos[[0,2]]*[-1,1])*ratio*0.5*[-1,1]
        ssusi_tpos[[0,2]] += (tpos[2]-ssusi_tpos[2])
        ssusi_tpos[[1,3]] = mean(ssusi_tpos[[1,3]])+total(ssusi_tpos[[1,3]]*[-1,1])*ratio*0.5*[-1,1]

        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

        dmsp_plot_time_range = time_double(dmsp_time_list[pid])
        tmp = min(times-mean(dmsp_plot_time_range), abs=1, time_id)
        npx = n_elements(mlt_images[0,0,*])
        mlt_image = reform(mlt_images[time_id,0:npx*0.5-1,*])
        ssusi_time = times[time_id]
        zzs = bytscl(mlt_image, min=ssusi_zrange[0], max=ssusi_zrange[1], top=color_top)
        ssusi_time_range = reform(lim.time_range[time_id,*])

        sgtv, fltarr(5,10)+128, ct=ct, position=tpos
        sgtv, zzs, ct=ct, position=ssusi_tpos

        ; Add labels, etc.
        ; Draw axes.
        plot, [-1,0], [-1,1], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos
        plots, [0,0], [-1,0], color=sgcolor('silver')
        ;plots, [0,0], [0,-1], color=sgcolor('silver')
        ;plots, [0,1], [0,0], color=sgcolor('silver')
        

        ; circles for ytickv.
        foreach yminor, ytick_minor, val_id do begin
            rr = (yminor-min_mlat)/(90-min_mlat)
            sc_xs = rr*cos(angles)
            sc_ys = rr*sin(angles)
            linestyle = 1
            index = where(ytickv eq yminor, count)
            if count ne 0 then linestyle = 0
            oplot, sc_xs,sc_ys, linestyle=linestyle, color=sgcolor('silver')
        endforeach


        ; lines for xickv.
        foreach xminor, xtick_minor, val_id do begin
            linestyle = 1
            index = where(xtickv eq xminor, count)
            if count ne 0 then linestyle = 0
            
            tt = (xminor*15-90)*constant('rad')
            sc_xs = [0,1]*cos(tt)
            sc_ys = [0,1]*sin(tt)

            plots, sc_xs,sc_ys, data=1, linestyle=linestyle, color=sgcolor('silver')
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
    
    
    ;---Add vector for velocity.
        vel_color = sgcolor('salmon')
        dmsp_orbit_time_range = ssusi_time_range
        mlt_var = prefix+'mlt'
        mlat_var = prefix+'mlat'
        mlon_var = prefix+'mlon'
        line_color = sgcolor('silver')
        mlts = get_var_data(mlt_var, in=dmsp_orbit_time_range, times=the_times)
        mlats = get_var_data(mlat_var, in=dmsp_orbit_time_range, times=the_times)
        mlons = get_var_data(mlon_var, at=the_times)
        ; use aacgm mlon.
        mlts = aacgm_mlon2mlt(mlons, the_times)
        index = where(abs(mlats) ge 50 and abs(mlats) le 65 and mlts le 0)
        dmsp_orbit_time_range = minmax(the_times[index])
        dmsp_orbit_time_range = dmsp_orbit_time_range-(dmsp_orbit_time_range mod 60)+[0,60]
        mlts = get_var_data(mlt_var, in=dmsp_orbit_time_range, times=the_times)
        mlats = get_var_data(mlat_var, in=dmsp_orbit_time_range, times=the_times)
        mlons = get_var_data(mlon_var, at=the_times)
        rad = constant('rad')
        tts = (mlts*15-90)*rad
        rrs = abs((90-abs(mlats))/total(mlat_range*[-1,1]))
        sc_xs = rrs*cos(tts)
        sc_ys = rrs*sin(tts)
        
        
        vel_var = prefix+'v_dmsp_xyz'
        vel_xyz = get_var_data(vel_var, in=dmsp_orbit_time_range, times=vel_times)
        vel_fac = vel_xyz

        fmlt = get_var_data(mlt_var, at=vel_times)
        fmlat = get_var_data(mlat_var, at=vel_times)
        theta = (fmlt*15-90)*constant('rad')
        n_hat = -[[cos(theta)],[sin(theta)]]
        w_hat = [[cos(theta-0.5*!dpi)],[sin(theta-0.5*!dpi)]]
        ; out, west, down.
        vel_fac[*,0] = n_hat[*,0]*vel_xyz[*,0]+n_hat[*,1]*vel_xyz[*,1]
        vel_fac[*,1] = w_hat[*,0]*vel_xyz[*,0]+w_hat[*,1]*vel_xyz[*,1]

        vel_scale = ychsz*3/1   ; convert 1 km/s to normal.
        x0s = interpol(sc_xs, the_times, vel_times)
        y0s = interpol(sc_ys, the_times, vel_times)
        x1s = x0s+(n_hat[*,0]*vel_fac[*,0]+w_hat[*,0]*vel_fac[*,1])*vel_scale
        y1s = y0s+(n_hat[*,1]*vel_fac[*,0]+w_hat[*,1]*vel_fac[*,1])*vel_scale

        foreach time, vel_times, time_id do begin
            if time_id mod 10 ne 0 then continue
            plots, [x0s[time_id],x1s[time_id]], [y0s[time_id],y1s[time_id]], color=vel_color
        endforeach
        plots, x1s, y1s, color=vel_color;, normal=1


        ; Add scale.
        tx = tpos[0]+xchsz*1
        ty = tpos[3]-ychsz*1.5
        tmp = convert_coord(tx,ty, normal=1, to_data=1)
        tx = tmp[0]
        ty = tmp[1]
        txs = tx+[0,vel_scale]
        tys = ty+[0,0]
        plots, txs, tys, color=vel_color
        msg = '1 km/s'
        tmp = convert_coord(mean(txs),ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]+ychsz*0.3
        xyouts, tx,ty, msg,normal=1, charsize=label_size, alignment=0.5, color=vel_color

        

    ;---Add SC track.
        oplot, sc_xs, sc_ys, color=line_color
        
;        tx = sc_xs[-1]
;        ty = sc_ys[-1]
;        tmp = convert_coord(tx,ty, data=1, to_normal=1)
;        tx = tmp[0]
;        ty = tmp[1]
;        msg = strupcase(sc_name+' '+probe)
;        xyouts, tx-xchsz*1,ty-ychsz*1.2,normal=1, msg, alignment=0.7, color=dmsp_color

        minor_times = make_bins(dmsp_orbit_time_range, 60, inner=1)
        minor_tts = interpol(tts, the_times, minor_times)
        minor_rrs = interpol(rrs, the_times, minor_times)
        plots, minor_rrs*cos(minor_tts), minor_rrs*sin(minor_tts), psym=8, symsize=0.5, color=line_color

        ; The times.
        major_times = smkarthm(dmsp_orbit_time_range[0],dmsp_orbit_time_range[1], 180, 'dx')
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
    

        
        
    ;---Add panel label.
        mlat_var = prefix+'mlat'
        mlat = get_var_data(mlat_var, at=ssusi_time)
        hem = (mlat ge 0)? 'North': 'South'
        msgs = [$
            'd-'+string(pid+1,format='(I0)')+') '+hem,$
            'DMSP '+strupcase(probe)+' | '+strjoin(time_string(ssusi_time_range,tformat='hh:mm'),'-')+' UT']
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.2
        xyouts, tx,ty+ychsz*1,normal=1, msgs[0]
        xyouts, tx,ty,normal=1, msgs[1], charsize=label_size



    ;---Add RBSP footpoint.
        internal_model = 'dipole'
        psyms = [1,4,6,7]
;        foreach external_model, ['t89','t96','t01','t04s'], id do begin
;            rbsp_prefix = 'rbspa_'
;            suffix = '_'+internal_model+'_'+external_model+'_'+strlowcase(hem)
;            config_time = time_double('2015-03-17/15:50')
;            fmlats = get_var_data(rbsp_prefix+'fmlat'+suffix, at=config_time)
;            fmlts = get_var_data(rbsp_prefix+'fmlt'+suffix, at=config_time)
;
;            tts = (fmlts*15-90)*rad
;            rrs = abs((90-abs(fmlats))/total(mlat_range*[-1,1]))
;            sc_xs = rrs*cos(tts)
;            sc_ys = rrs*sin(tts)
;
;            plots, sc_xs, sc_ys, data=1, color=sgcolor('purple'), psym=psyms[id]
;        endforeach
;        stop
        
    endfor


    ; Colorbar.
    if keyword_set(horizontal_cb) then begin
        ssusi_cbpos = my_poss[*,0]
        ssusi_cbpos[1] = ssusi_cbpos[3]+ychsz*0.4
        ssusi_cbpos[3] = ssusi_cbpos[1]+ychsz*0.4
        ssusi_cbpos[0] += xchsz*1
        ssusi_cbpos[2] -= xchsz*1
        ztitle = 'SSUSI '+strupcase(ssusi_id)+' ('+ssusi_unit+')'
        zrange = ssusi_zrange
        zticklen = uniform_ticklen/(ssusi_cbpos[3]-ssusi_cbpos[1])/fig_size[1]
        sgcolorbar, findgen(color_top), horizontal=1, $
            ztitle=ztitle, zrange=zrange, ct=ct, position=ssusi_cbpos, zticklen=zticklen, zminor=5
    endif else begin
        ssusi_cbpos = my_poss[*,1]
        ssusi_cbpos[0] = ssusi_cbpos[2]+xchsz*1
        ssusi_cbpos[2] = ssusi_cbpos[0]+xchsz*0.7
        ztitle = 'SSUSI '+strupcase(ssusi_id)+' ('+ssusi_unit+')'
        zrange = ssusi_zrange
        zticklen = uniform_ticklen/(ssusi_cbpos[3]-ssusi_cbpos[1])/fig_size[1]
        sgcolorbar, findgen(color_top), horizontal=0, $
            ztitle=ztitle, zrange=zrange, ct=ct, position=ssusi_cbpos, zticklen=zticklen, zminor=5
    endelse
    
    return, my_poss

end

function fig_saps_v02_config_panel, my_pos, event_info=event_info, $
    pos_xrange=pos_xrange, pos_yrange=pos_yrange, config_time=config_time, _extra=ex

;---SC position.
    probe = 'a'
    prefix = 'rbsp'+probe+'_'
    tpos = my_pos
    fig_label = ' '
    label_size = 0.8
    rbsp_info = dictionary($
        'sc_color', sgcolor('red'), $
        'sc_name', strupcase('rbsp'), $
        'probe', probe )


    the_time = config_time
    fline_thick = (!d.name eq 'X')? 4: 8
    ;fline_colors = sgcolor(['red','blue'])
    fline_colors = sgcolor(['tan','gray'])

    external_model = 't96'
    internal_model = 'dipole'
    dir = 1
    xrange = pos_xrange
    yrange = pos_yrange


    line_mlats = [45,50,55,60,65,70]
    line_mlats = [48,50,52,54,56,58,60]
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

    par_var = geopack_read_par(the_time+[-1,1]*60*20, model=external_model, t89_par=t89_par)
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

    xminor = 5
    xstep = 0.5
    xtickv = make_bins(xr, xstep, inner=1)
    xticks = n_elements(xtickv)-1
    xtickn = string(xtickv,format='(F4.1)')
    ;    xtitle = 'SM R!DXY!N (Re)'
    xtickn[-1] = 'SM R!DXY!N (Re)        '
    xtitle = ' '
    yminor = 5
    ystep = 0.5
    ytickv = make_bins(yr, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    ytitle = 'SM Z (Re)'
    
    xgrids = make_bins(xr, xstep*2, inner=1)
    ygrids = make_bins(yr, ystep*2, inner=1)

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


    tmp = lets_add_earth(xrange=xr, yrange=yr, npoint=80)
    
    ; Add MLT.
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = 'MLT = '+string(the_mlt+24,format='(F4.1)')+' h'
    xyouts, tx,ty,msg, normal=1;, charsize=label_size
    msg = 'Model = '+strupcase(external_model)
    ty = tpos[3]-ychsz*2
    xyouts, tx,ty,msg, normal=1;, charsize=label_size
    
end




function fig_saps_v02, plot_dir, event_info=event_info, test=test

    version = 'v02'
    id = '2015_0317'
    if n_elements(event_info) eq 0 then event_info = low_lshell_outflow_load_data(id)
    
    time_range = event_info.time_range
    foreach probe, ['a','b'] do begin
        foreach species, ['p','o','e'] do begin
            var = rbsp_read_en_spec(time_range, probe=probe, species=species)
        endforeach
    endforeach
    
    
    dmsp_times = dictionary()
    dmsp_times['f18'] = ['2015-03-17/15:15','2015-03-17/16:40']

    foreach probe, dmsp_times.keys() do begin
        time_range = time_double(dmsp_times[probe])
        coord = 'dmsp_xyz'
        vel_var = dmsp_read_ion_vel(time_range, probe=probe, coord=coord)
        b_var = dmsp_read_bfield(time_range, probe=probe, coord=coord, read_b0=1)
        mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, id='lbhs', update=0)
        mlat_vars = dmsp_read_mlat_vars(time_range, probe=probe)
        mlt_images = get_var_data(mlt_image_var, in=time_range, times=times, limits=lim)
    endforeach
    
    

;---Figure out figure size.
    ssusi_xsize = 1.5
    ssusi_ysize = 3
    ssusi_margins = [1,1,6,1]
    xpads = [7,1]
    rbsp_xsize = 3.5
    xpans = [rbsp_xsize,ssusi_xsize*[1,1]]


    if n_elements(plot_dir) eq 0 then plot_dir = event_info.plot_dir
    plot_file = join_path([plot_dir,'fig_saps_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    poss = panel_pos(plot_file, nxpan=3, nypan=1, $
        pansize=[rbsp_xsize,ssusi_ysize], xpans=xpans, xpad=xpads, margins=ssusi_margins, fig_size=fig_size)



    sgopen, plot_file, size=fig_size;, magnify=1.5
    ssusi_poss = poss[*,1:2]
    ssusi_poss = fig_saps_v02_ssusi_panel(ssusi_poss, event_info=event_info)
    
    rbsp_pos = poss[*,0]
    rbsp_poss = fig_saps_v02_rbsp_panel(rbsp_pos, event_info=event_info)
    

    

;    fig_size = [6,6]
;    sgopen, plot_file, size=fig_size
;    my_pos = sgcalcpos(1)
;
;
;    ; Plot config.
;    pos_xrange = [-0.5,-3.5]
;    pos_yrange = [-1,1]*1.5
;    config_time = time_double('2015-03-17/15:50')
;    tpos = fig_saps_v02_config_panel(my_pos, event_info=event_info, config_time=config_time, pos_xrange=pos_xrange, pos_yrange=pos_yrange)
;    
;    ; Plot DMSP data.
;    tr_list = list()
;    tr_list.add, ['2015-03-17/16:00','2015-03-17/16:40']
;    tr_list.add, ['2015-03-17/15:15','2015-03-17/15:45']
;    colors = sgcolor(['purple','red'])
;    prefix = 'dmspf18_'
;    
;    sgopen, 2
;    tpos = sgcalcpos(1)
;    xrange = [45,55]
;    yrange = [-1,3]
;    plot, xrange, yrange, $
;        xstyle=1, ystyle=1, $
;        nodata=1, noerase=1, position=tpos
;        
;    foreach tr, tr_list, id do begin
;        tr = time_double(tr)
;        vel_var = prefix+'v_dmsp_xyz'
;        vel = get_var_data(vel_var, in=tr, times=times)
;        mlat_var = prefix+'mlat'
;        mlat = get_var_data(mlat_var, at=times)
;        oplot, abs(mlat), vel[*,1], color=colors[id]
;    endforeach
;    stop
;    
;
;    time_range = event_info.time_range
;
;    sgopen, plot_file, size=fig_size
;
;


    

    if keyword_set(test) then stop
    sgclose
    return, plot_file

end



print, fig_saps_v02(event_info=event_info, test=1)
end