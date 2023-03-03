;+
; An overview figure, showing the conjunction of wave and arc
;-

function fig_2013_0501_fig_overview_v04, plot_file, event_info=event_info

test = 1



;---Load data.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
    plasma_param = event_info['plasma_param']



;---Settings.
    prefix = event_info['prefix']
    probe = event_info['probe']
    time_range = event_info['time_range']
    snapshot_time = event_info['snapshot_time']

    poleward_expansion_time_range = time_double(['2013-05-01/07:38','2013-05-01/07:41'])
    poleward_expansion_mlat_range = [63.5,65.2]

    sc_color = event_info['sc_color']
    pe_color = sgcolor('orange')
    pe_symsize = 0.4
    pe_psym = 8
    

    tmp = smkarthm(0,2*!dpi,20,'n')
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs,tys, fill=1
    

;---PF FAC and spec.
    pf_var = prefix+'pf_fac_plot'
    var = duplicate_var(prefix+'pfdot0_fac_map', output=pf_var)
    yrange = [-100,500]
    ystep = 200
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 2
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', [0,200,400]
    pf_labels = 'S!D'+['||',tex2str('perp')+','+['west','out']]+'!N'
    options, var, 'labels', [' ',' ',' ']

    pf_spec_var = prefix+'pf_spec_plot'
    var = duplicate_var(prefix+'pfdot0_fac_mor_spec_1', output=pf_spec_var)
    pflux_setting = event_info['pflux_setting']
    filter = pflux_setting['filter']
    f_filter = minmax(1d/filter)
    yrange = f_filter*1e3
    log_yrange = alog10(yrange)
    log_yrange = [ceil(log_yrange[0]),floor(log_yrange[1])]
    log_ytickv = make_bins(log_yrange,1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    index = lazy_where(ytickv, '()', yrange)
    constant = ytickv[index]
    foreach tx, ytickv, ii do begin
        if tx eq 1 then begin
            ytickn[ii] = '1'
        endif else if tx eq 10 then begin
            ytickn[ii] = '10'
        endif
    endforeach
    ;ytickn[0:*:2] = ' '
    yminor = 10
    var = pf_spec_var
;    options, var, 'constant', [10,100,1000]
    options, var, 'yrange', yrange
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'ytickname', ytickn
    options, var, 'constant', constant
    
    vars = [pf_var]
    foreach var, vars do begin
        labels = get_setting(var,'labels')
        foreach label, labels, ii do begin
            index = strpos(label,'FAC')
            if index[0] eq -1 then continue
            labels[ii] = strmid(label,index[0]+4)
        endforeach
        options, var, 'labels', labels
    endforeach


;---Fpt ASI count.
    model_setting = event_info['model_setting']
    external_model = model_setting['external_model']
    internal_model = model_setting['internal_model']
    fmlt_var = prefix+'fmlt_'+internal_model
    fmlat_var = prefix+'fmlat_'+internal_model
    models = model_setting['models']
    model_index = where(models eq external_model)

    mlt_image_var = 'thg_asf_mlt_image_rect'
    get_data, mlt_image_var, times, mlt_images, limits=lim
    mlt_bins = lim.mlt_bins
    mlat_bins = lim.mlat_bins
    ntime = n_elements(times)
    asi_count = fltarr(ntime)

    fmlt = (get_var_data(fmlt_var, at=times))[*,model_index]
    fmlat = (get_var_data(fmlat_var, at=times))[*,model_index]
    dmlat = 1d      ; deg.
    dmlt = dmlat/15 ; h.

    foreach time, times, time_id do begin
        mlt_index = lazy_where(mlt_bins, '[]', fmlt[time_id]+[-1,1]*dmlt*0.5, count=count)
        if count eq 0 then continue
        mlat_index = lazy_where(mlat_bins, '[]', fmlat[time_id]+[-1,1]*dmlat*0.5, count=count)
        if count eq 0 then continue

        asi_count[time_id] = mean(mlt_images[time_id,mlt_index,mlat_index])
    endforeach
    asi_count_var = 'asi_count'
    store_data, asi_count_var, times, asi_count
    add_setting, asi_count_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', ' ', $
        'unit', '#')
    yrange = [0,2.5e4]
    ytickv = [1,2]*1e4
    yticks = n_elements(ytickv)-1
    yminor = 5
    var = asi_count_var
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', ytickv


    ; ASI keogram.
    asi_keo_var = 'thg_asf_keo'
    asi_mlat_range = [60,67]
    asi_mlat_step = 5
    yrange = asi_mlat_range
    ystep = asi_mlat_step
    ytickv = make_bins(yrange,ystep,inner=1)
    yminor = ystep
    yticks = n_elements(ytickv)-1
    var = asi_keo_var
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    constant = ytickv[lazy_where(ytickv, '()', yrange)]
    options, var, 'constant', constant
    
    asi_zrange = [0,1e4]
    zrange = asi_zrange
    options, var, 'zrange', zrange
    options, var, 'no_color_scale', 1

    
    
;---Plot.
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,'fig_2013_0501_fig_overview_v04.pdf'])
    if keyword_set(test) then plot_file = 0


    left_vars = [asi_keo_var,asi_count_var,pf_var,pf_spec_var]
    nleft_var = n_elements(left_vars)
    fig_letters = letters(nleft_var+1)
    fig_labels = fig_letters+') '+['Keo','Aurora','FAC S','S!D||!N spec']


    ; Calc the figure size.
    ; use the xyrange of SC position and # of left panels.
    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz
    pos_xrange = [2,-10]
    pos_yrange = [-1,1]*2.5
    pos_xpan_size = 3d
    pos_aspect_ratio = abs(total(pos_yrange*[-1,1])/total(pos_xrange*[-1,1]))
    pos_ypan_size = pos_aspect_ratio*pos_xpan_size

    asi_aspect_ratio = 0.6
    asi_ypan_size = asi_aspect_ratio*pos_xpan_size
    
    ; Now get the panel size of left panels.
    ; ASI and POS panels + 2*ycharsize for each = same height of 3 left panels.
    ypans = fltarr(nleft_var)+1
    ypans[[1,3]] = 0.7
    left_xpan_size = pos_xpan_size
    left_ypan_size = (pos_ypan_size+asi_ypan_size+abs_ychsz*4)/total(ypans[0:2])
    
    xpans = [left_xpan_size,pos_xpan_size]
    ypads = 0.4
    


    margins = [10,4,2,1]
    label_size = 0.8

    tmp = panel_pos(plot_file, pansize=[left_xpan_size,left_ypan_size], $
        xpans=xpans, ypans=ypans, ypads=ypads, margins=margins, $
        fig_size=fig_size)

    poss_left = reform(tmp[*,0,*])
    poss_right = reform(tmp[*,1,*])

    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    uniform_ticklen = -ychsz*0.15*fig_size[0]


;---POS panel.
    fig_label = fig_letters[-1]+')'
    tpos = poss_right[*,2]
    pos_xsize = (tpos[2]-tpos[0])*fig_size[0]
    pos_ysize = pos_xsize*pos_aspect_ratio
    tpos[3] = tpos[1]+pos_ysize/fig_size[1]
    tpos[1] += ychsz*4
    tpos[3] += ychsz*4
    
    the_time = snapshot_time
    fline_thick = (keyword_set(test))? 2: 4
    fline_colors = sgcolor(['red','blue'])
    
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

    par_var = geopack_read_par(time_range, model=model, t89_par=t89_par)
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
    
    
    
    ; Get the poleward expansion info.
    line_mlats = poleward_expansion_mlat_range
    nline = n_elements(line_mlats)
    pe_lines = list()
    pe_points = list()
    ps = geopack_recalc(the_time)

    r_gsm_var = prefix+'r_gsm'
    r_sm = transpose(get_var_data(r_gsm_var, at=the_time))
    r_sm = cotran(r_gsm, the_time, 'gsm2sm')
    the_mlt = pseudo_mlt(r_sm)
    line_mlts = the_mlt+dblarr(nline)

    par_var = geopack_read_par(time_range, model=model, t89_par=t89_par)
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
        pe_lines.add, fline
        
        sc_pos = r_sm
        srotate, sc_pos, (24-the_mlt)*15*rad, 2
        sc_mlat = atan(sc_pos[2],sc_pos[0])*deg
        fline_mlat = atan(fline[*,2],fline[*,0])*deg
        index = where(fline_mlat ge 0)
        pe_points.add, reform(sinterpol(fline[index,*], fline_mlat[index], sc_mlat))
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


    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    xminor = 2
    xtickv = make_bins(xr, xminor, inner=1)
    xticks = n_elements(xtickv)
    xtickn = string(xtickv,format='(I0)')
;    xtickn[-1] = 'SM X (Re)        '
    xtitle = 'SM R!DXY!N (Re)'
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
;            tx = tmp[0]+xchsz*0.5
;            ty = tmp[1]-ychsz*1
;            msg = 'F/MLat: '+strtrim(string(fmlat,format='(F4.1)'))+' deg'
;            xyouts, tx,ty,normal=1, msg, charsize=label_size, color=sgcolor('red')
;            ty = tmp[1]-ychsz*1
;            msg = string(the_mlt,format='(F4.1)')
;            if the_mlt le 0 then msg = string(the_mlt+24,format='(F4.1)')
;            msg = 'MLT: '+strtrim(msg,2)+' h'
;            xyouts, tx,ty,normal=1, msg, charsize=label_size, color=sgcolor('red')
        endif
    endforeach


    foreach line, lines do begin
        oplot, line[*,0], line[*,2]
    endforeach
    

    ; Add points, need to rotate a little.
    p1 = (pe_points[0])[[0,2]]
    p2 = (pe_points[1])[[0,2]]
    
    
    sc_p0 = sc_pos[[0,2]]
    sc_p1 = (sc_neighbors[0])[[0,2]]
    sc_p2 = (sc_neighbors[1])[[0,2]]
    dir1 = sunitvec(sc_p1-sc_p2)
    theta = atan(dir1[1],dir1[0])
    dir2 = [cos(theta+!pi*0.5),sin(theta+!pi*0.5)]

    dp1 = p1-sc_p0
    foreach ii, [1,-1] do begin
        the_dir = dir2*ii
        if sang(the_dir, dp1, deg=1) lt 90 then break
    endforeach
    p1 = sc_p0+snorm(dp1)*the_dir
    
    dp2 = p2-sc_p0
    foreach ii, [1,-1] do begin
        the_dir = dir2*ii
        if sang(the_dir, dp2, deg=1) lt 90 then break
    endforeach
    p2 = sc_p0+snorm(dp2)*the_dir
    
    
    oplot, [p1[0],p2[0]],[p1[1],p2[1]], color=pe_color, psym=-pe_psym, symsize=pe_symsize
    pe_dis = snorm(p1[[0,1]]-p2[[0,1]])
    msg = string(pe_dis,format='(F3.1)')+' Re'
    tmp = convert_coord(sc_p0[0], sc_p0[1], data=1, to_normal=1)
    tx = tmp[0]-xchsz*4
    ty = tmp[1]+ychsz*0.
    xyouts, tx,ty,normal=1, msg, color=pe_color, charsize=label_size
    
    
    ; Add SC.
    plots, sc_p0[0], sc_p0[1], psym=6, color=sc_color, symsize=label_size
    tx = tmp[0]+xchsz*1
    ty = tmp[1]-ychsz*0.3
    msg = 'RBSP-'+strupcase(probe)
    xyouts, tx,ty,normal=1, msg, charsize=label_size, color=sc_color
    

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


    

;---tplot panels.
    plot_poss = poss_left

    for ii=0,nleft_var-1 do begin
        tpos = plot_poss[*,ii]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, left_vars[ii], 'xticklen', xticklen
        options, left_vars[ii], 'yticklen', yticklen
    endfor

    tplot, left_vars, trange=time_range, position=plot_poss, noerase=1
    for ii=0,nleft_var-1 do begin
        tpos = plot_poss[*,ii]
        fig_label = fig_labels[ii]
        tx = 2*xchsz
        ty = tpos[3]-ychsz*0.6
        xyouts, tx,ty,fig_label, normal=1
    endfor
    
    
    ; Add labels for pflux.
    pid = where(left_vars eq pf_var, count)
    if count ne 0 then begin
        tpos = plot_poss[*,pid]
        xrange = time_range
        yrange = get_setting(pf_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            nodata=1, noerase=1, position=tpos
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = 'Norm to 100 km altitude!CFiltered in '+$
            string(f_filter[0]*1e3,format='(I0)')+'mHz-'+$
            string(f_filter[1],format='(I0)')+'Hz'
        xyouts, tx,ty,normal=1, msg, charsize=label_size
        
        tmp = convert_coord(xrange[0],0, data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.4
        msg = 'To N-hem'
        xyouts, tx,ty,normal=1, msg, charsize=label_size
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]-ychsz*0.8
        msg = 'To S-hem'
        xyouts, tx,ty,normal=1, msg, charsize=label_size
        
        colors = get_setting(pf_var, 'colors')
        labels = pf_labels
        foreach color, colors, ii do begin
            tx = tpos[2]-xchsz*((2-ii)*5+1)
            ty = tpos[3]-ychsz
            xyouts, tx,ty,labels[ii], normal=1, color=color, alignment=1
        endforeach
    endif
    
    

    ; Add label for keogram.
    pid = where(left_vars eq asi_keo_var, count)
    if count ne 0 then begin
        tpos = plot_poss[*,pid]
        xrange = time_range
        yrange = get_setting(asi_keo_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            nodata=1, noerase=1, position=tpos
        ; Add Fpt.
        fmlat_var = prefix+'fmlat_'+internal_model
        fmlat = get_var_data(fmlat_var, times=times, in=time_range)
        models = model_setting['models']
        model_index = where(models eq external_model)
        the_fmlat = fmlat[*,model_index]
        oplot, times, the_fmlat, linestyle=2, color=sc_color
        tx = times[0]
        ty = the_fmlat[0]
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.3
        msg = 'RBSP-'+strupcase(probe)
        xyouts, tx,ty, normal=1, msg, charsize=label_size, color=sc_color
        
        txs = poleward_expansion_time_range
        tys = poleward_expansion_mlat_range
        oplot, txs,tys, color=pe_color, psym=-pe_psym, symsize=pe_symsize
    endif
    
    
    ; Add label for ASI count.
    pid = where(left_vars eq asi_count_var, count)
    if count ne 0 then begin
        tpos = plot_poss[*,pid]
        xrange = time_range
        yrange = get_setting(asi_count_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            nodata=1, noerase=1, position=tpos
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        tmp = string(dmlat,format='(I0)')
        msg = tmp+tex2str('times')+tmp+' deg around!CRBSP-'+strupcase(probe)+' footpoint'; ('+string(external_model)+')'
        xyouts, tx,ty,normal=1, msg, charsize=label_size
    endif
    
    ; Add vertical grids and times.
    bar_times = make_bins(time_range,600, inner=1)
    timebar, bar_times, linestyle=1, color=sgcolor('black')
    timebar, snapshot_time, color=sgcolor('red'), linestyle=3

    

;---ASI snapshot.
    tpos = poss_right[*,0]
    fig_label = 'a-1)'
    
    tpos[2] = 1-xchsz*8
    cbpos = tpos
    cbpos[0] = tpos[2]+xchsz*0.5
    cbpos[2] = cbpos[0]+xchsz*0.8
    horizontal = 0
    
    model_setting = event_info['model_setting']
    external_model = model_setting['external_model']
    internal_model = model_setting['internal_model']
    models = model_setting['models']
    asi_setting = event_info['asi_setting']
    site = (asi_setting['sites'])[0]
    
    asi_mlt_range = [-2,0.5]
    asi_ct = 49
    top_color = 254

    xtitle = 'MLT (h)'
    xrange = asi_mlt_range+24
    xstep = 0.5
    xtickv = make_bins(xrange,xstep,/inner)
    xticks = n_elements(xtickv)-1
    xminor = 5
    xtickn = string(xtickv,format='(F4.1)')
    foreach tx, xtickv, ii do begin
        if tx lt 24 then continue
        xtickn[ii] = strtrim(string(xtickv[ii]-24,format='(F4.1)'),2)
    endforeach
;    xtickn[0] = 'MLT (h)    '
    

    ytitle = 'MLat (deg)'
    yrange = asi_mlat_range
    ystep = 5
    ytickv = make_bins(yrange,ystep,/inner)
    yticks = n_elements(ytickv)-1
    yminor = 5


    zrange = asi_zrange
    ztitle = 'ASI Count (#)'
    zstep = zrange[1]/2
    ztickv = make_bins(zrange,zstep)
    zticks = n_elements(ztickv)-1
    zticklen = uniform_ticklen/(cbpos[2]-cbpos[0])/fig_size[0]
    zminor = 5

    mlt_image_var = 'thg_asf_mlt_image_rect'
    get_data, mlt_image_var, times, mlt_images
    tmp = min(times-snapshot_time, index, abs=1)
    time = times[index]
    mlt_image = reform(mlt_images[index,*,*])

    mlt_bins = get_setting(mlt_image_var, 'mlt_bins')
    mlt_index = lazy_where(mlt_bins, '[]', asi_mlt_range)
    mlt_bins = mlt_bins[mlt_index]
    mlt_image = mlt_image[mlt_index,*]
    mlat_bins = get_setting(mlt_image_var, 'mlat_bins')
    mlat_index = lazy_where(mlat_bins, '[]', asi_mlat_range)
    mlat_bins = mlat_bins[mlat_index]
    mlt_image = mlt_image[*,mlat_index]
    

    zzs = bytscl(mlt_image, min=zrange[0],max=zrange[1], top=top_color)
    sgtv, zzs, ct=asi_ct, position=tpos, resize=1
    sgcolorbar, findgen(top_color), horizontal=horizontal, $
        ztitle=ztitle, zrange=zrange, ct=asi_ct, position=cbpos, $
        zticks=zticks, ztickv=ztickv, zminor=zminor, zticklen=zticklen

    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    ; Add grid.
    plot, xrange, yrange, $
        xstyle=1, xlog=0, xrange=xrange, xtitle='', xtickformat='(A1)', xtickv=xtickv, xticks=xticks, xminor=1, xticklen=1, xgridstyle=1, xtickname=xtickn, $
        ystyle=1, ylog=0, yrange=yrange, ytitle='', ytickformat='(A1)', ytickv=ytickv, yticks=yticks, yminor=1, yticklen=1, ygridstyle=1, $
        position=tpos, nodata=1, noerase=1, ynozero=1, color=sgcolor('gray')

    ; Draw axes.
    plot, xrange, yrange, $
        xstyle=1, xlog=0, xrange=xrange, xtitle=xtitle, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
        ystyle=1, ylog=0, yrange=yrange, ytitle=ytitle, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
        position=tpos, nodata=1, noerase=1, ynozero=1
        
        
    ; Add footpoint.
    fmlt = get_var_data(prefix+'fmlt_'+internal_model, at=time)+24
    fmlat = get_var_data(prefix+'fmlat_'+internal_model, at=time)
    model_index = where(models eq external_model)
    color = sc_color
    fmlt = fmlt[model_index]
    fmlat = fmlat[model_index]

    tmp = convert_coord(fmlt, fmlat, /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    plots, tx,ty,normal=1, psym=6, symsize=label_size, color=color
    msg = 'RBSP-B'
    xyouts, tx,ty+ychsz*0.5,normal=1, msg, color=color, alignment=0.5, charsize=label_size

    ; Add label.
    tx = tpos[0]+xchsz*1
    ty = tpos[1]+ychsz*0.3
    msg = strupcase(site)+' '+time_string(time,tformat='YYYY-MM-DD/hh:mm:ss')+' UT'
    xyouts, tx,ty,normal=1, msg, color=sgcolor('black'), charsize=label_size
    
    tx = tpos[0]-xchsz*4
    ty = tpos[3]-ychsz*0.8
    msg = fig_label
    xyouts, tx,ty,normal=1, msg


;    j_ver = get_var_data('thg_j_ver', at=time, limits=lim)
;    pixel_mlons = lim.mlon_grids
;    pixel_mlats = lim.mlat_grids
;    pixel_mlts = mlon2mlt(pixel_mlons, time)+24
;    index = where(pixel_mlts gt xrange[0] and pixel_mlts lt xrange[1] and $
;        pixel_mlats gt yrange[0] and pixel_mlats lt yrange[1], npixel)
;    xxs = pixel_mlts[index]
;    yys = pixel_mlats[index]
;    zzs = j_ver[index]
;    zrange = [-1,1]*1e5
;    ccs = bytscl(zzs, min=zrange[0], max=zrange[1])
;    dc = 80
;    ccs[where(zzs ge 0)] = 128+dc
;    ccs[where(zzs lt 0)] = 128-dc
;    ct = 70
;    tmp = smkarthm(0,2*!dpi,30,'n')
;    txs = cos(tmp)
;    tys = sin(tmp)
;    usersym, txs, tys
;    
;    symszs = (abs(zzs/zrange[1]))^0.25*0.5
;    for ii=0,npixel-1 do begin
;        cc = sgcolor(ccs[ii], ct=ct)
;        symsz = symszs[ii]
;        ;symsz = 0.5
;        ;cc = (zzs[ii] ge 0)? sgcolor('blue'): sgcolor('red')
;        plots, xxs[ii], yys[ii], color=cc, psym=8, symsize=symsz
;    endfor
;    stop

    
    
    sub_poss = sgcalcpos(1,2, xpad=3, region=poss_right[*,3], margins=[3,0,6,0])

;---S power.
    tpos = sub_poss[*,0]
    fig_label = 'd-1)'
    

    pf_insitu_var = prefix+'pfdot0_fac'
    pf_mor_var = pf_insitu_var+'_mor'
    scale_info = pflux_setting['scale_info']
    stplot_mor_new, pf_insitu_var, newname=pf_mor_var, frequency=1, scaleinfo=scale_info
    get_data, pf_mor_var, uts, dat, val
    pf_mor_info_var = pf_mor_var+'_fft_info'
    get_data, pf_mor_info_var, tmp, info
    nrec = n_elements(uts)
    gws = total(dat,1)/nrec
    psd = 2*info.c_tau*info.dt/info.cdelta*gws*1e3
    fs = info.fs*1e3*0.5
    
    log_xrange = alog10(minmax(psd))
    log_xrange = [floor(log_xrange[0]),ceil(log_xrange[1])]
    xrange = 10.^log_xrange
    log_xtickv = make_bins(log_xrange, 1)
    xtickv = 10d^log_xtickv
    xticks = n_elements(xtickv)-1
    xminor = 10
    xtitle = 'PS ('+tex2str('mu')+'W/m!U2!N)/Hz'
    xtickn = '10!U'+string(log_xtickv,format='(I0)')
    xtickn[0:*:2] = ' '
    xlog = 1
    yrange = get_setting(pf_spec_var, 'yrange')
    ytickv = get_setting(pf_spec_var, 'ytickv')
    yticks = get_setting(pf_spec_var, 'yticks')
    ytickn = get_setting(pf_spec_var, 'ytickname')
    yminor = get_setting(pf_spec_var, 'yminor')
    constant = get_setting(pf_spec_var, 'constant')
    ylog = 1
    ytitle = get_setting(pf_spec_var, 'ytitle')

    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    
    plot, psd, fs, $
        ystyle=1, ylog=ylog, yrange=yrange, $
        ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, $
        xstyle=1, xlog=xlog, xrange=xrange, $
        xtickv=xtickv, xticks=xticks, xminor=xminor, xtitle=xtitle, xtickname=xtickn, $
        xticklen=xticklen, yticklen=yticklen, $
        ytickformat='(A1)', $
        position=tpos, noerase=1
        
        
    ; Add f_max.
    f_target = 10
    f_color = sgcolor('salmon')
    index = lazy_where(fs, '()', f_target*[0.5,2])
    txs = psd[index]
    tys = fs[index]
    tmp = max(txs,index)
    f_max = tys[index]
    plots, xrange, f_max+[0,0], linestyle=2, color=f_color
    tmp = convert_coord(xrange[0],f_max, data=1, to_normal=1)
    tx = tmp[0]+xchsz*0.5
    ty = tmp[1]+ychsz*0.2
    xyouts, tx,ty,normal=1, string(f_max,format='(F4.1)')+' mHz', charsize=label_size, color=f_color

    pid = where(left_vars eq pf_spec_var, count)
    tmp_tpos = tpos
    if count ne 0 then begin
        tpos = plot_poss[*,pid]
        xrange = time_range
        yrange = get_setting(pf_spec_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        oplot, xrange, f_max+[0,0], color=f_color, linestyle=2
    endif
    tpos = tmp_tpos

    
    ; Add labels.
    tx = tpos[2]-xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = fig_label
    xyouts, tx,ty, msg, normal=1, alignment=1
    
    
;---E/B ratio.
    tpos = sub_poss[*,1]
    fig_label = 'd-2)'

    ebr_var = prefix+'ebr_component'
    get_data, ebr_var, ps, ebr
    fs = 1d3/ps
    xrange = [100,1e5]
    log_xrange = alog10(xrange)
    log_xrange = [ceil(log_xrange),floor(log_xrange)]
    log_xtickv = make_bins(log_xrange, 1)
    xtickv = 10d^log_xtickv
    xticks = n_elements(xtickv)-1
    xminor = 10
    xtitle = 'E/B ratio (km/s)'
    xtickn = '10!U'+string(log_xtickv,format='(I0)')
;    xtickn[0] = xtitle
;    xtitle = ' '
    xlog = 1
    yrange = get_setting(pf_spec_var, 'yrange')
    ytickv = get_setting(pf_spec_var, 'ytickv')
    yticks = get_setting(pf_spec_var, 'yticks')
    ytickn = get_setting(pf_spec_var, 'ytickname')
    yminor = get_setting(pf_spec_var, 'yminor')
    constant = get_setting(pf_spec_var, 'constant')
    ylog = 1
    ytitle = get_setting(pf_spec_var, 'ytitle')

    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    plot, ebr, fs, $
        ystyle=9, ylog=ylog, yrange=yrange, $
        ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, $
        xstyle=1, xlog=xlog, xrange=xrange, $
        xtickv=xtickv, xticks=xticks, xminor=xminor, xtitle=xtitle, xtickname=xtickn, $
        xticklen=xticklen, yticklen=yticklen, $
        ytickformat='(A1)', $
        position=tpos, noerase=1
    axis, yaxis=1, save=1, $
        ylog=ylog, yrange=yrange, yticklen=yticklen, ytitle=ytitle, $
        ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn
    foreach ty, constant do begin
        oplot, xrange, ty+[0,0], linestyle=1
    endforeach
    

    ; Add E/B ratio.
    suf = ''
    va = plasma_param['va'+suf]
    f_gi = plasma_param['fg'+suf]
    vi = plasma_param['vi'+suf]
    vf_fac = plasma_param['vf_exb_fac']
    vf = vf_fac[1]  ; about 100 km/s.

    va = plasma_param['va'+suf]/sqrt(2.3)
    vf = vf_fac[2]

    f_sc = 1/ps
    ebr_theory = va*sqrt(1+(f_sc/f_gi*(vi/vf))^2)
    oplot, ebr_theory, fs, color=sgcolor('salmon'), linestyle=0


    ; Add Va.
    plots, va+[0,0], yrange, color=sgcolor('red'), linestyle=1
    tx = va
    ty = yrange[0]
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]
    ty = tpos[3]-ychsz*1
    msg = 'v!DA!N'
    xyouts, tx,ty,normal=1, msg, alignment=0.5, color=sgcolor('red')



    ; Add labeling.
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = fig_label
    xyouts, tx,ty,normal=1, msg

    
    if keyword_set(test) then stop
    sgclose
    
    return, plot_file
    

end

print, fig_2013_0501_fig_overview_v04(event_info=event_info)
end