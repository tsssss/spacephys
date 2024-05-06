;+
;-

function fig_dmsp_overview_v02, test=test

    id = '2013_0317'
    version = 'v02'
    plot_dir = join_path([googledir(),'works','2024_daps',id,'plot'])
    ;time_range = time_double(['2013-03-17/09:15','2013-03-17/09:45'])
    time_range = time_double(['2013-03-17/09:20','2013-03-17/09:45'])
    ssusi_time = time_double('2013-03-17/09:32')
    ssusi_id = 'energy'
    probe = 'f18'
    prefix = 'dmsp'+probe+'_'

    r1_times = list()
    r2_times = list()

    r1r2_boundary = ['2013-03-17/09:25:50','2013-03-17/09:39:40']
    r2_times.add, ['2013-03-17/09:25:06', r1r2_boundary[0]]
    r1_times.add, [r1r2_boundary[0], '2013-03-17/09:27:32']
    
    r1_times.add, ['2013-03-17/09:38:40', r1r2_boundary[1]]
    r2_times.add, [r1r2_boundary[1], '2013-03-17/09:42:52']
    
    r1_color = sgcolor('wheat')    
    r2_color = sgcolor('lavender')
    r1_color = sgcolor('orange')
    r2_color = sgcolor('purple')

    db_xyz_var = dmsp_read_bfield_madrigal(time_range, probe=probe)
    b0_xyz_var = dmsp_read_bfield_madrigal(time_range, probe=probe, read_b0=1)
    v_var = dmsp_read_ion_vel_madrigal(time_range, probe=probe)
    dmsp_mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, id=ssusi_id)
    mlat_vars = dmsp_read_mlat_vars(time_range, probe=probe, errmsg=errmsg)
    mlat_var = mlat_vars[0]
    mlt_var = mlat_vars[1]
    mlon_var = mlat_vars[2]

    ele_spec_var = dmsp_read_en_spec(time_range, probe=probe, species='e', errmsg=errmsg)
    ion_spec_var = dmsp_read_en_spec(time_range, probe=probe, species='p', errmsg=errmsg)
    ele_eflux_var = dmsp_read_eflux(time_range, probe=probe, species='e', errmsg=errmsg)

    b_vec = get_var_data(db_xyz_var, times=times)
    index = where(finite(snorm(b_vec)), count)
    b_vec = sinterpol(b_vec[index,*], times[index], times)
    store_data, db_xyz_var, times, b_vec

    b_vec = get_var_data(b0_xyz_var, times=times)
    bmag_var = prefix+'bmag'
    store_data, bmag_var, times, snorm(b_vec)
    add_setting, bmag_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', '|B|', $
        'unit', 'nT')
    options, bmag_var, 'ytitle', '|B| (nT)'

    var = prefix+'v_dmsp_xyz'
    options, var, 'labels', 'V'+constant('xyz')
    yrange = [-1.5,5.5]
    ytickv = [-1,2,5]
    yminor = 3
    set_ytick, var, yrange=yrange, ytickv=ytickv, yminor=yminor
    options, var, 'constant', ytickv

    var = prefix+'db_dmsp_xyz'
    options, var, 'labels', 'dB'+constant('xyz')
    yrange = [-1,1]*700
    ytickv = [-1,0,1]*500
    yminor = 5
    set_ytick, var, yrange=yrange, ytickv=ytickv, yminor=yminor
    options, var, 'constant', ytickv




;---Prepare plot.
    base = 'fig_dmsp_overview_'+id+'_'+version+'.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(test) then plot_file = 0
    ypans = [1,1.5]
    pansize = [1,0.5]*5.5
    margins = [1,0,6.5,1.5]
    all_poss = panel_pos(2, ypans=ypans, margins=margins, pansize=pansize, fig_size=fig_size)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    uniform_ticklen = -ychsz*0.3*fig_size[1]


;---SSUSI
    tpos = all_poss[*,0]

    tmp = smkarthm(0,2*!dpi,20,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs, circ_ys, fill=1

    color_top = 254
    ct = 70
    ssusi_zrange = [-1,1]*20
    ssusi_min_mlat = 50d

    mlt_range = [-1d,1]*6
    min_mlat = 50d
    mlat_range = [min_mlat,90]

    nangle = 50
    dangle = !dpi/nangle
    angles = make_bins(mlt_range*12, dangle)

    ; xlabels.
    xtickn_pos = (90.-min_mlat)/(90.-min_mlat+2)
    xrange = mlt_range
    xstep = 6d
    xtickv = smkarthm(xrange[0]+2,xrange[1]-2, 2, 'dx')
    xtick_minor = make_bins(xrange, 1, inner=1)
    xtickn = xtickv
    index = where(xtickv lt 0, count)
    if count ne 0 then xtickn[index] += 24
    xtickn = string(xtickn,format='(I02)')
    
    ; ylabels.
    ytick_pos = (1d/12-1)*!dpi
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange, 10, inner=1)
    ytick_minor = make_bins(yrange, ystep, inner=1)
    ytickn = string(ytickv,format='(I0)')


    ssusi_tpos = tpos
    ratio = abs((90-ssusi_min_mlat)/(90-min_mlat))
    ssusi_tpos[[0,2]] = mean(ssusi_tpos[[0,2]])+total(ssusi_tpos[[0,2]]*[-1,1])*ratio*0.5*[-1,1]
    ssusi_tpos[[0,2]] += (tpos[2]-ssusi_tpos[2])
    ssusi_tpos[[1,3]] = mean(ssusi_tpos[[1,3]])+total(ssusi_tpos[[1,3]]*[-1,1])*ratio*0.5*[-1,1]

    label_size = 0.8
    label_size = 1d
    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]


    mlt_images = get_var_data(dmsp_mlt_image_var, times=times, limits=lim)
    npx = n_elements(mlt_images[0,*,0])
    tmp = min(times-ssusi_time, time_index, abs=1)
    mlt_image = reform(mlt_images[time_index,*,0:npx*0.5-1])
    ssusi_time = times[time_index]
    zzs = bytscl(mlt_image, min=ssusi_zrange[0], max=ssusi_zrange[1], top=color_top)
    ssusi_time_range = reform(lim.time_range[time_index,*])
    ssusi_unit = lim.unit

    sgtv, zzs, position=tpos, ct=ct
    
    ; Add labels, etc.
    ; Draw axes.
    plot, [-1,1], [-1,0], nodata=1, noerase=1, $
        xstyle=5, ystyle=5, position=tpos
    plots, [0,0], [-1,0], color=sgcolor('silver')
    ;plots, [0,0], [0,-1], color=sgcolor('silver')
    plots, [-1,1], [0,0], color=sgcolor('silver')
    

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
    vel_color = sgcolor('green')
    dmsp_orbit_time_range = time_range
    mlt_var = prefix+'mlt'
    mlat_var = prefix+'mlat'
    mlon_var = prefix+'mlon'
    line_color = sgcolor('silver')
    mlts = get_var_data(mlt_var, in=dmsp_orbit_time_range, times=the_times)
    mlats = get_var_data(mlat_var, in=dmsp_orbit_time_range, times=the_times)
    mlons = get_var_data(mlon_var, at=the_times)
    ; use aacgm mlon.
    mlts = aacgm_mlon2mlt(mlons, the_times)
    ;index = where(abs(mlats) ge 50 and abs(mlats) le 65 and mlts le 0)
    index = where(abs(mlats) ge min_mlat and mlts le 6)
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

    data_scale = 1
    norm_scale = ychsz*2.5
    vel_scale = norm_scale/data_scale   ; convert 1 km/s to normal.
    x0s = interpol(sc_xs, the_times, vel_times)
    y0s = interpol(sc_ys, the_times, vel_times)


    ntt = n_elements(tts)
    sc_zs = fltarr(ntt)
    sc_rs = [[sc_xs],[sc_ys],[sc_zs]]
    x_hat = sunitvec(sinterpol((sc_rs[1:-1,*]-sc_rs[0:-2,*]),(tts[1:-1]+tts[0:-2])*0.5, tts))
    z_hat = fltarr(ntt,3) & z_hat[*,2] = -1
    y_hat = vec_cross(z_hat,x_hat)
    x1s = x0s+(x_hat[*,0]*vel_xyz[*,0]+y_hat[*,0]*vel_xyz[*,1])*vel_scale
    y1s = y0s+(x_hat[*,1]*vel_xyz[*,0]+y_hat[*,1]*vel_xyz[*,1])*vel_scale

    foreach time, vel_times, time_id do begin
        if time_id mod 10 ne 0 then continue
        plots, [x0s[time_id],x1s[time_id]], [y0s[time_id],y1s[time_id]], color=vel_color
    endforeach
    plots, x1s, y1s, color=vel_color;, normal=1


    ; Add scale.
    nn = 2
    tx = tpos[0]+xchsz*2
    ty = tpos[3]-ychsz*1.5
    tmp = convert_coord(tx,ty, normal=1, to_data=1)
    tx = tmp[0]
    ty = tmp[1]
    txs = tx+[0,norm_scale*nn]
    tys = ty+[0,0]
    plots, txs, tys, color=vel_color
    foreach tx,txs do begin
        tmp = convert_coord(tx,tys[0],data=1,to_normal=1)
        plots, tmp[0]+[0,0],tmp[1]+[-1,1]*ychsz*0.1, normal=1
    endforeach
    msg = 'v!D'+tex2str('perp')+'!N '+string(data_scale*nn,format='(I0)')+' km/s'
    tmp = convert_coord(mean(txs),ty, data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]+ychsz*0.5
    xyouts, tx,ty, msg,normal=1, charsize=label_size, alignment=0.5, color=vel_color




;---Add SC track.
    thick = (keyword_set(test))? 4: 8
    foreach tr, r1_times, ii do begin
        tr = time_double(tr)
        index = where_pro(the_times, '[]', tr)
        uts = sort_uniq([the_times[index],tr])
        txs = sinterpol(sc_xs, the_times, uts)
        tys = sinterpol(sc_ys, the_times, uts)
        plots, txs,tys,data=1, color=r1_color, thick=thick
        tx = mean(txs)
        ty = mean(tys)
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        if ii eq 0 then begin
            tx = tmp[0]-xchsz*1
            ty = tmp[1]-ychsz*0.5
        endif
        msg = 'R1'
        xyouts, tx,ty-ychsz*1,msg, normal=1, color=r1_color
    endforeach
    
    foreach tr, r2_times, ii do begin
        tr = time_double(tr)
        index = where_pro(the_times, '[]', tr)
        uts = sort_uniq([the_times[index],tr])
        txs = sinterpol(sc_xs, the_times, uts)
        tys = sinterpol(sc_ys, the_times, uts)
        plots, txs,tys,data=1, color=r2_color, thick=thick
        
        tx = mean(txs)
        ty = mean(tys)
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        if ii eq 0 then begin
            tx = tmp[0]+xchsz*0.5
            ty = tmp[1]+ychsz*0.2
        endif
        msg = 'R2'
        xyouts, tx,ty-ychsz*1,msg, normal=1, color=r2_color
    endforeach
    

    oplot, sc_xs, sc_ys, color=line_color

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
        if major_times[ii] lt time_double('2013-03-17/09:32') then begin
            ty = tmp[1]-ychsz*1.2
        endif
        plots, tmp[0], tmp[1], normal=1, psym=8, symsize=0.5, color=dmsp_color
        xyouts, tx,ty,normal=1, msg, alignment=0.5, color=dmsp_color, charsize=label_size
    endforeach


;---Add labels.
    tx = tpos[0]+xchsz*0.5
    ty = tpos[1]+ychsz*0.5
    msg = 'a) South | '+strupcase('dmsp-'+probe)
    xyouts, tx,ty,msg, normal=1

    ty = tpos[1]+ychsz*1.5
    rr = ychsz/(tpos[3]-tpos[1])*(1-0)*0.5
    tmp = convert_coord(tx+xchsz*1,ty+ychsz*0.7, normal=1, to_data=1)
    tx0 = tmp[0]
    ty0 = tmp[1]
    txs = tx0+circ_xs*rr
    tys = ty0+circ_ys*rr
    plots, txs,tys, data=1
    tmp = [-1,1]*rr*cos(45*constant('rad'))    
    plots, tx0+tmp,ty0+tmp, data=1
    plots, tx0+tmp,ty0-tmp, data=1
    msg = 'B!D0'
    xyouts, tx+xchsz*3,ty+ychsz*0.3,msg, normal=1, charsize=1.2
    
    
    tx = 0.5
    ty = 0
    tmp = convert_coord(tx,ty,data=1,to_normal=1)
    tx = tmp[0]-xchsz*1
    ty = tmp[1]+ychsz*0.2
    msg = 'DAPS'
    xyouts, tx,ty,msg, normal=1, charsize=1.2
    
    tx = -0.8
    ty = -0.4
    tmp = convert_coord(tx,ty,data=1,to_normal=1)
    tx = tmp[0]-xchsz*1
    ty = tmp[1]+ychsz*0.2
    msg = 'SAPS'
    xyouts, tx,ty,msg, normal=1, charsize=1.2
    

;---Colorbar.
    if keyword_set(horizontal_cb) then begin
        ssusi_cbpos = tpos
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
        ssusi_cbpos = tpos
        ssusi_cbpos[0] = ssusi_cbpos[2]+xchsz*1
        ssusi_cbpos[2] = ssusi_cbpos[0]+xchsz*0.8
        ztitle = 'SSUSI '+strupcase(ssusi_id)+' ('+ssusi_unit+')'
        zrange = ssusi_zrange
        zticklen = uniform_ticklen/(ssusi_cbpos[2]-ssusi_cbpos[0])/fig_size[0]
        sgcolorbar, findgen(color_top), horizontal=0, $
            ztitle=ztitle, zrange=zrange, ct=ct, position=ssusi_cbpos, zticklen=zticklen, zminor=5
    endelse
    

;---Other panels.
    r1_color = sgcolor('wheat')    
    r2_color = sgcolor('lavender')

    margins = [10,6,6.5,1]
    plot_vars = prefix+['e_en_spec','p_en_spec','v_dmsp_xyz','db_dmsp_xyz']
    nplot_var = n_elements(plot_vars)
    tpos = all_poss[*,1] & tpos[0] = 0 & tpos[2] = 1
    poss = sgcalcpos(nplot_var, margins=margins, region=tpos)
    tplot_options, 'tickinterval', 5*60
    var_labels = prefix+['mlt_plot','mlat','bmag']
    copy_data, prefix+'mlt', prefix+'mlt_plot'
    var = prefix+'mlt_plot'
    mlt = get_var_data(var, times=times)
    index = where(mlt le 0, count)
    if count ne 0 then begin
        mlt[index] += 24d
        store_data, var, times, mlt
    endif
    options, var, 'ytitle', 'MLT (h)'
    var = prefix+'mlat'
    options, var, 'ytitle', 'MLat (deg)'
    
    vars = prefix+['e','p']+'_en_spec'
    options, vars, 'ytitle', 'Energy!C(eV)'

    tpos = poss[*,nplot_var-2]
    tpos[1] = poss[1,nplot_var-1]
    yrange = [0,1]
    set_axis, position=tpos, xrange=time_range, yrange=yrange
    foreach tr, r1_times do begin
        tr = time_double(tr)
        polyfill, tr[[0,1,1,0,0]], yrange[[0,0,1,1,0]], color=r1_color
        tx = mean(tr)
        ty = yrange[1]
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*1
        msg = 'R1!C'
        mlt = get_var_data(prefix+'mlt', at=tr[0])
        if mlt le 0 then msg += 'U' else msg += 'D'
        xyouts, tx,ty,normal=1, msg, alignment=0.5
    endforeach
    foreach tr, r2_times do begin
        tr = time_double(tr)
        polyfill, tr[[0,1,1,0,0]], yrange[[0,0,1,1,0]], color=r2_color
        tx = mean(tr)
        ty = yrange[1]
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*1
        msg = 'R2!C'
        mlt = get_var_data(prefix+'mlt', at=tr[0])
        if mlt ge 0 then msg += 'U' else msg += 'D'
        xyouts, tx,ty,normal=1, msg, alignment=0.5
    endforeach
    
    for pid=0,nplot_var-1 do begin
        tpos = poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, plot_vars[pid], 'xticklen', xticklen
        options, plot_vars[pid], 'yticklen', yticklen
    endfor
    
    foreach var, prefix+['e','p']+'_en_spec' do begin
        pid = where(plot_vars eq var, count)
        if count eq 0 then continue
        tpos = poss[*,pid]
        cbpos = tpos
        cbpos[0] = cbpos[2]+xchsz*1
        cbpos[2] = cbpos[0]+xchsz*0.8
        zticklen = uniform_ticklen/(cbpos[2]-cbpos[0])/fig_size[0]
        options, var, 'zticklen', zticklen
        ;options, var, 'colorscale', cbpos
    endforeach

    tplot, plot_vars, trange=time_range, position=poss, var_label=var_labels, vlab_margin=9, noerase=1
    fig_labels = letters([0,nplot_var]+1)+') '+['e-','H+','Ion Vel','dB']
    for pid=0,nplot_var-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*8.5
        ty = tpos[3]-ychsz*0.7
        msg = fig_labels[pid]
        xyouts, tx,ty,msg, normal=1
        
        if pid eq 0 then begin
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            msg = strupcase('dmsp '+probe)
            xyouts, tx,ty,msg, normal=1, color=sgcolor('white')
        endif
    endfor
    
    tpos = poss[*,-1]
    yrange = [0,1]
    set_axis, position=tpos, xrange=time_range, yrange=yrange
    msg_times = list()
    msg_times.add, minmax(time_double([r2_times[0],r1_times[0]]))
    msg_times.add, minmax(time_double([r2_times[1],r1_times[1]]))
    msgs = ['Dawn','Dusk']
    foreach msg, msgs, tid do begin
        tx = mean(msg_times[tid])
        ty = yrange[1]
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*1
        if tid eq 1 then ty = tpos[1]+ychsz*0.5
        xyouts, tx,ty,normal=1, msg, alignment=0.5
    endforeach
    

    if keyword_set(test) then stop
    sgclose
    return, plot_file
    
end

print, fig_dmsp_overview_v02(test=0)
end