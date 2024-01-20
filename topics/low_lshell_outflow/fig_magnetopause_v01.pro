;+
; 
;-

function scale_pp_time_to_zz, pp_time, time_range

    return, bytscl(pp_time, min=time_range[0], max=time_range[1], top=254)
    
end

function fig_magnetopause_v01_config_panel, config_pos, $
    pos_xrange=pos_xrange, pos_yrange=pos_yrange, event_info=event_info

    time_range = event_info.time_range
    pp_symsize = event_info.pp.symsize
    pp_psym = event_info.pp.psym
    if pp_psym eq 8 then begin
        tmp = smkarthm(0,2*!dpi,50,'n')
        circ_xs = cos(tmp)
        circ_ys = sin(tmp)
    endif
    
;---Setup coord.
    tpos = config_pos
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.2*fig_size[1]
    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    xr = pos_xrange
    yr = pos_yrange

    xminor = 2
    xtickv = make_bins(xr, xminor, inner=1)
    xticks = n_elements(xtickv)
    xtickn = string(xtickv,format='(I0)')
    xtitle = 'SM X (Re)'
    ;xtickn[-1] = 'SM X (Re)        '
    ;xtitle = ' '
    yminor = 2
    ytickv = make_bins(yr, yminor, inner=1)
    yticks = n_elements(ytickv)
    ytitle = 'SM Y (Re)'
    
    xstep = 4
    xgrids = make_bins(xr, xstep, inner=1)
    ystep = 4
    ygrids = make_bins(yr, ystep, inner=1)
    
    
    ; Add axis.
    plot, xr, yr, $
        xstyle=5, xrange=xr, xtitle=xtitle, xminor=xminor, xticks=xticks, xtickv=xtickv, xtickn=xtickn, $
        ystyle=5, yrange=yr, ytitle=ytitle, yminor=yminor, yticks=yticks, ytickv=ytickv, $
        nodata=1, noerase=1, position=tpos, iso=1, $
        xticklen=xticklen, yticklen=yticklen
        
    
    ; Add circles.
    label_size = 0.8
    tmp = smkarthm(0,2*!dpi, 40, 'n')
    txs = cos(tmp)
    tys = sin(tmp)
    foreach rr, make_bins([2,10],2) do begin
        oplot, rr*txs, rr*tys, linestyle=1, color=sgcolor('silver')
    endforeach
    
    ; Add lines.
    foreach xx, [0] do begin
        oplot, xx+[0,0], yr, linestyle=1, color=sgcolor('silver')
    endforeach
    foreach yy, [0] do begin
        oplot, xr, yy+[0,0], linestyle=1, color=sgcolor('silver')
    endforeach
    
    ; Add earth.
    tmp = smkarthm(0,2*!dpi, 40, 'n')
    txs = cos(tmp)
    tys = sin(tmp)
    polyfill, txs>0, tys, color=sgcolor('white')
    polyfill, txs<0, tys, color=sgcolor('grey')
    plots, txs, tys
    
    

    ; Add axis.
    plot, xr, yr, $
        xstyle=1, xrange=xr, xtitle=xtitle, xminor=xminor, xticks=xticks, xtickv=xtickv, xtickn=xtickn, $
        ystyle=1, yrange=yr, ytitle=ytitle, yminor=yminor, yticks=yticks, ytickv=ytickv, $
        nodata=1, noerase=1, position=tpos, iso=1, $
        xticklen=xticklen, yticklen=yticklen

    pp_infos = get_var_data('pp_infos')
    ct = event_info.pp.ct
    

;---THEMIS.
    themis_probes = event_info.themis_probes
    themis_probes = 'a'
    foreach probe, themis_probes do begin
        prefix = 'th'+probe+'_'
        mlt = get_var_data(prefix+'mlt', times=times)
        lshell = get_var_data(prefix+'lshell', at=times)
        tt = (mlt*15+180)*constant('rad')
        xxs = lshell*cos(tt)
        yys = lshell*sin(tt)
        oplot, xxs, yys, linestyle=0, color=sgcolor('silver')
    endforeach
    
    tmp = smkarthm(0,2*!dpi,50,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    foreach probe, themis_probes do begin
        prefix = 'th'+probe+'_'
        r_var = prefix+'r_sm'
        foreach pp_info, pp_infos do begin
            if pp_info.mission ne 'themis' then continue
            if pp_info.probe ne probe then continue
            pp_time = pp_info.time
            zz = scale_pp_time_to_zz(pp_time, time_range)
            color = sgcolor(zz, ct=ct)
            color = sgcolor(pp_info.color)
            mlt = get_var_data(prefix+'mlt', at=pp_time)
            lshell = get_var_data(prefix+'lshell', at=pp_time)
            tt = (mlt*15+180)*constant('rad')
            xxs = lshell*cos(tt)
            yys = lshell*sin(tt)
            if pp_psym eq 8 then begin
                usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1
                plots, xxs, yys, psym=pp_psym, color=color, symsize=pp_symsize
                usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=0
                plots, xxs, yys, psym=pp_psym, symsize=pp_symsize
            endif else begin
                plots, xxs, yys, psym=pp_psym, color=color, symsize=pp_symsize
                plots, xxs, yys, psym=pp_psym, symsize=pp_symsize
            endelse
        endforeach
    endforeach

    tx = 2
    ty = 9
    msg = 'THEMIS'
    xyouts, tx,ty, data=1, msg;, charsize=label_size


;---GOES.
    goes_probes = event_info.goes_probes
    ;goes_probes = '13'
    foreach probe, goes_probes do begin
        prefix = 'g'+probe+'_'
        mlt = get_var_data(prefix+'mlt', times=times)
        lshell = get_var_data(prefix+'lshell', at=times)
        tt = (mlt*15+180)*constant('rad')
        xxs = lshell*cos(tt)
        yys = lshell*sin(tt)
        oplot, xxs, yys, linestyle=0, color=sgcolor('silver')
    endforeach
    
    foreach probe, goes_probes do begin
        prefix = 'g'+probe+'_'
        foreach pp_info, pp_infos do begin
            if pp_info.mission ne 'goes' then continue
            if pp_info.probe ne probe then continue
            pp_time = pp_info.time
            zz = scale_pp_time_to_zz(pp_time, time_range)
            color = sgcolor(zz, ct=ct)
            color = sgcolor(pp_info.color)
            mlt = get_var_data(prefix+'mlt', at=pp_time)
            lshell = get_var_data(prefix+'lshell', at=pp_time)
            tt = (mlt*15+180)*constant('rad')
            xxs = lshell*cos(tt)
            yys = lshell*sin(tt)
            if pp_psym eq 8 then begin
                usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1
                plots, xxs, yys, psym=pp_psym, color=color, symsize=pp_symsize
                usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=0
                plots, xxs, yys, psym=pp_psym, symsize=pp_symsize
            endif else begin
                plots, xxs, yys, psym=pp_psym, color=color, symsize=pp_symsize
                plots, xxs, yys, psym=pp_psym, symsize=pp_symsize
            endelse
        endforeach
    endforeach

    tx = 5
    ty = 4
    msg = 'GOES'
    xyouts, tx,ty, data=1, msg;, charsize=label_size


;--Add label.
    tx = tpos[0]-xchsz*4
    ty = tpos[3]-ychsz*0.8
    msg = 'a)'
    xyouts, tx,ty,normal=1, msg
    
    foreach rr, make_bins([2,10],2) do begin
        tt = -30*constant('rad')
        tx = rr*cos(tt)
        ty = rr*sin(tt)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*0.3
        msg = 'L = '+string(rr,format='(I0)')
        xyouts, tx,ty,msg, normal=1, alignment=0.5, charsize=label_size
    endforeach
    
end


function fig_magnetopause_v01_right_panels, right_pos, event_info=event_info

    time_range = event_info.time_range
    b_var = omni_read_sw_b(time_range, coord='sm')
    options, b_var, 'yrange', [-1,1]*35
    options, b_var, 'ytickv', [-1,0,1]*20
    options, b_var, 'yticks', 2
    options, b_var, 'yminor', 3
    options, b_var, 'ylog', 0
    options, b_var, 'ystyle', 1
    options, b_var, 'constant', [-1,0,1]*20
    options, b_var, 'labels', 'B'+constant('xyz')
    
    
    p_var = omni_read_sw_p(time_range)
    options, p_var, 'yrange', [1,70]
    options, p_var, 'ytickv', [1,10]*3
    options, p_var, 'yticks', 1
    options, p_var, 'yminor', 9
    options, p_var, 'ylog', 1
    options, p_var, 'ystyle', 1
    options, p_var, 'constant', [1,10]*3
    options, p_var, 'labels', 'SW P!Ddyn!N'
    
    sw_b_var = omni_read_sw_b(time_range, coord='sm')
    options, sw_b_var, 'labels', 'SW B!D'+constant('xyz')
    
    
    dst_var = omni_read_symh(time_range)
    options, dst_var, 'ystyle', 1
    options, dst_var, 'labels', 'SymH'
    options, dst_var, 'ytitle', '(nT)'
    
    themis_probes = event_info.themis_probes
    goes_probes = event_info.goes_probes
    
    
    b_vars = ['tha','the','thd','g13','g15']+'_b_sm'
    options, b_vars, 'yrange', [-1,1]*200
    options, b_vars, 'ytickv', [-1,0,1]*200
    options, b_vars, 'yticks', 2
    options, b_vars, 'yminor', 4
    options, b_vars, 'constant', [-1,0,1]*100
    foreach var, b_vars do begin
        prefix = get_prefix(var)
        mission_probe = strmid(prefix,0,strlen(prefix)-1)
        options, var, 'labels', 'B!D'+constant('xyz')+'!N '+[strupcase(mission_probe),'','']
    endforeach
    
    
    v_vars = ['tha','the','thd']+'_u_sm'
    options, v_vars, 'yrange', [-1,1]*600
    options, v_vars, 'ytickv', [-1,0,1]*400
    options, v_vars, 'yticks', 2
    options, v_vars, 'yminor', 4  
    options, v_vars, 'constant', [-2,-1,0,1,2]*200  
    options, v_vars, 'labels', 'V!D'+constant('xyz')
    foreach var, v_vars do begin
        prefix = get_prefix(var)
        mission_probe = strmid(prefix,0,strlen(prefix)-1)
        options, var, 'labels', 'V!D'+constant('xyz')+'!N '+[strupcase(mission_probe),'','']
    endforeach

    plot_vars = [dst_var,b_var,p_var, 'tha_'+['b_sm','u_sm'],'g'+goes_probes+'_b_sm']
    nvar = n_elements(plot_vars)
    fig_labels = letters(minmax(1+findgen(nvar+1)))+')'
    
    
    ypans = fltarr(nvar)+1
    right_poss = sgcalcpos(nvar, position=right_pos, ypans=ypans)
    
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.2*fig_size[1]
    for pid=0,nvar-1 do begin
        tpos = right_poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, plot_vars[pid], 'xticklen', xticklen
        options, plot_vars[pid], 'yticklen', yticklen
    endfor
    
    tplot, plot_vars, trange=time_range, position=right_poss, noerase=1, novtitle=1
    
    for pid=0,nvar-1 do begin
        tpos = right_poss[*,pid]
        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*0.8
        msg = fig_labels[pid]
        xyouts, tx,ty, normal=1, msg
    endfor

    pp_infos = get_var_data('pp_infos')
    ct = event_info.pp.ct
    pp_symsize = event_info.pp.symsize
    pp_psym = event_info.pp.psym
    if pp_psym eq 8 then begin
        tmp = smkarthm(0,2*!dpi,50,'n')
        circ_xs = cos(tmp)
        circ_ys = sin(tmp)
    endif

    pp_colors = get_color(n_elements(pp_infos))
    foreach pp_info, pp_infos, pp_id do begin
        probe = pp_info.probe
        prefix = pp_info.mission_probe+'_'
        pp_time = pp_info.time
        zz = scale_pp_time_to_zz(pp_time, time_range)
        color = sgcolor(pp_info.color)
        dens_var = prefix+'b_sm'
        
        pid = where(plot_vars eq dens_var, count)
        if count eq 0 then continue
        tpos = right_poss[*,pid]
        yr = get_setting(dens_var, 'yrange')
        xr = time_range
        ylog = 0
        plot, xr, yr, $
            xrange=xr, yrange=yr, xstyle=5, ystyle=5, ylog=ylog, $
            position=tpos, nodata=1, noerase=1
        
        tmp = convert_coord(pp_time, 0, data=1, to_normal=1)
        txs = tmp[0]+[0,0]
        if probe eq 'a' then begin
            tys = [tpos[3],right_poss[1,pid+1]]
        endif else begin
            tys = tpos[[1,3]]
        endelse
        plots, txs,tys, linestyle=0, normal=1, color=color
        
;        tmp = get_var_data(dens_var, at=pp_time)
;        ty = tmp[2]
;        if pp_psym eq 8 then begin
;            usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1
;            plots, pp_time, ty, psym=pp_psym, color=color, symsize=pp_symsize
;            usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=0
;            plots, pp_time, ty, psym=pp_psym, symsize=pp_symsize
;        endif else begin
;            plots, pp_time, ty, psym=pp_psym, color=color, symsize=pp_symsize
;            plots, pp_time, ty, psym=pp_psym, symsize=pp_symsize
;        endelse
    endforeach
    
    
    
;    pid = where(plot_vars eq dst_var, count)
;    if count ne 0 then begin
;        tpos = right_poss[*,pid]
;        yticklen = get_setting(dst_var, 'yticklen')
;        yr = get_setting(dst_var, 'yrange')
;        xr = time_range
;        plot, xr, yr, $
;            xrange=xr, yrange=yr, xstyle=5, ystyle=5, $
;            position=tpos, nodata=1, noerase=1
;        
;        yrange = [5.5,10]
;        ytickv = [6,8,10]
;        yticks = n_elements(ytickv)-1
;        yminor = 4
;        ytitle = 'L-shell'
;        axis, yaxis=1, save=1, ystyle=1, $
;            yticklen=yticklen, yrange=yrange, $
;            ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle
;        
;        foreach pp_info, pp_infos do begin
;            pp_time = pp_info.time
;            zz = scale_pp_time_to_zz(pp_time, time_range)
;            color = sgcolor(zz, ct=ct)
;            prefix = pp_info.mission_probe+'_'
;            tmp = get_var_data(prefix+'lshell', at=pp_time)
;            if pp_psym eq 8 then begin
;                usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1
;                plots, pp_time, tmp, psym=pp_psym, color=color, symsize=pp_symsize
;                usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=0
;                plots, pp_time, tmp, psym=pp_psym, symsize=pp_symsize
;            endif else begin
;                plots, pp_time, tmp, psym=pp_psym, color=color, symsize=pp_symsize
;                plots, pp_time, tmp, psym=pp_psym, symsize=pp_symsize
;            endelse
;        endforeach
;    endif


end



function fig_magnetopause_v01, plot_file, event_info=event_info, test=test

    version = 'v01'
    id = '2015_0317'
    event_info = low_lshell_outflow_load_data(id)
        
    time_range = event_info.time_range
    goes_probes = event_info.goes_probes
    themis_probes = event_info.themis_probes
    
    plot_dir = event_info['plot_dir']
    base = 'fig_magnetopause_'+time_string(time_range[0],tformat='YYYY_MMDD')+'_'+version+'.pdf'
    if n_elements(plot_file) eq 0 then plot_file = join_path([plot_dir,base])
    if keyword_set(test) then plot_file = 0

    mp_infos = list()
    coord = 'sm'

    dst_var = omni_read_symh(time_range)
    sw_p_var = omni_read_sw_p(time_range)
    sw_b_var = omni_read_sw_b(time_range, coord='sm')
    options, sw_b_var, 'constant', 0

;---THEMIS.
    foreach probe, themis_probes do begin
        prefix = 'th'+probe+'_'
        b_var = themis_read_bfield(time_range, probe=probe, coord=coord)
        v_var = themis_read_ion_vel(time_range, probe=probe, coord=coord)
        r_var = themis_read_orbit(time_range, probe=probe, coord=coord)
        mlt_var = themis_read_mlt(time_range, probe=probe)
        lshell_var = themis_read_lshell(time_range, probe=probe)
        ;en_var = themis_read_en_spec(time_range, probe=probe, species='e')
        ;en_var = themis_read_en_spec(time_range, probe=probe, species='p')
        options, b_var, 'constant', 0
        options, b_var, 'yrange', [-1,1]*200
        options, v_var, 'constant', 0
    endforeach


;---GOES.
    foreach probe, goes_probes do begin
        prefix = 'g'+probe+'_'
        b_var = prefix+'b_'+coord
        if check_if_update(b_var) then goes_read_bfield, time_range, probe=probe, coord=coord
        options, b_var, 'constant', 0
        options, b_var, 'yrange', [-1,1]*200
        r_var = goes_read_orbit(time_range, probe=probe, coord=coord)
        mlt_var = goes_read_mlt(time_range, probe=probe)
        lshell_var = goes_read_lshell(time_range, probe=probe)
    endforeach


;---Dst.
    dst_var = omni_read_symh(time_range)
    ystep = 50
    exact_yrange = minmax(get_var_data(dst_var))
    yrange = minmax(make_bins(exact_yrange, ystep))
    for ii=1,10 do begin
        ytickv = make_bins(exact_yrange, ystep*ii, inner=1)
        yticks = n_elements(ytickv)-1
        if yticks le 2 then break
    endfor
    options, dst_var, 'yrange', yrange
    options, dst_var, 'ytickv', ytickv
    options, dst_var, 'yticks', yticks
    options, dst_var, 'yminor', 5
    options, dst_var, 'constant', ytickv
    options, dst_var, 'ystyle', 9
    options, dst_var, 'labels', ' '
    options, dst_var, 'ytitle', 'SymH!C (nT)'
    

;---mp_infos.
    mp_infos.add, dictionary($
        'color', 'gold', $
        'mission_probe', 'tha', $
        'mission', 'themis', $
        'probe', 'a', $
        'time', time_double('2015-03-17/09:17:00') )
    mp_infos.add, dictionary($
        'color', 'red', $
        'mission_probe', 'tha', $
        'mission', 'themis', $
        'probe', 'a', $
        'time', time_double('2015-03-17/16:13:00') )
;    mp_infos.add, dictionary($
;        'mission_probe', 'the', $
;        'mission', 'themis', $
;        'probe', 'e', $
;        'time', time_double('2015-03-17/10:12:00') )
;    mp_infos.add, dictionary($
;        'mission_probe', 'the', $
;        'mission', 'themis', $
;        'probe', 'e', $
;        'time', time_double('2015-03-17/18:51:58') )
;    mp_infos.add, dictionary($
;        'mission_probe', 'thd', $
;        'mission', 'themis', $
;        'probe', 'd', $
;        'time', time_double('2015-03-17/11:57:00') )
;    mp_infos.add, dictionary($
;        'mission_probe', 'thd', $
;        'mission', 'themis', $
;        'probe', 'd', $
;        'time', time_double('2015-03-17/21:00:00') )
    mp_infos.add, dictionary($
        'color', 'brown', $
        'mission_probe', 'g13', $
        'mission', 'goes', $
        'probe', '13', $
        'time', time_double('2015-03-17/13:11:58') )
    mp_infos.add, dictionary($
        'color', 'teal', $
        'mission_probe', 'g13', $
        'mission', 'goes', $
        'probe', '13', $
        'time', time_double('2015-03-17/17:50:00') )
    mp_infos.add, dictionary($
        'color', 'sea_green', $
        'mission_probe', 'g15', $
        'mission', 'goes', $
        'probe', '15', $
        'time', time_double('2015-03-17/17:50:00') )
    store_data, 'pp_infos', 0, mp_infos


;---Plot settings.
    tmp = get_abs_chsz()
    abs_xchsz = tmp[0]
    abs_ychsz = tmp[1]

    pos_xrange = [10,-2]
    pos_yrange = [10,-7]
    margins = [6,4,8,1]

    pos_aspect_ratio = abs(total(pos_xrange*[-1,1])/total(pos_yrange*[-1,1]))
    pos_xpan_size = 3
    pos_ypan_size = pos_xpan_size/pos_aspect_ratio

    right_xpan_size = pos_xpan_size*1.5


    all_poss = panel_pos(plot_file, xpans=[pos_xpan_size,right_xpan_size], xpad=10, margins=margins, panid=[0,0], pansize=[pos_xpan_size,pos_ypan_size], fig_size=fig_size)

    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    

;---Pos.
    config_pos = all_poss[*,0]
    tpos = fig_magnetopause_v01_config_panel(config_pos, $
        pos_xrange=pos_xrange, pos_yrange=pos_yrange, event_info=event_info)
    right_pos = all_poss[*,1]
    poss = fig_magnetopause_v01_right_panels(right_pos, event_info=event_info)
    
    
    
    if keyword_set(test) then stop
    sgclose

    return, plot_file

    
end


print, fig_magnetopause_v01(test=0)
end