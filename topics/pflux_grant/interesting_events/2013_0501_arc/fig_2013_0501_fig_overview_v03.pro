;+
; An overview figure, showing the conjunction of wave and arc, and ion beams.
;-

function fig_2013_0501_overview_v03, plot_file, event_info=event_info

test = 1



;---Load data.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
    plasma_param = event_info['plasma_param']


;---Settings.
    prefix = event_info['prefix']
    probe = event_info['probe']
    time_range = event_info['time_range']
    snapshot_time = event_info['snapshot_time']
    
    xticklen_chsz = -0.2
    yticklen_chsz = -0.35

    e_var = prefix+'edot0_fac'
    var = e_var
    yrange = [-1,1]*190
    ystep = 100
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 5
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', [0,-100,100]
    get_data, var, limits=lim
    labels = lim.labels
    foreach label, labels, ii do begin
        index = strpos(label,'dot0')
        if index[0] eq -1 then continue
        labels[ii] = strmid(label,0,index)+strmid(label,index+4)
    endforeach
    options, var, 'labels', labels

    b_var = prefix+'b1_fac'
    var = b_var
    yrange = [-1,1]*60
    ystep = 50
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 5
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'constant', 0

    pf_var = prefix+'pfdot0_fac_map'
    var = pf_var
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

    pf_spec_var = prefix+'pfdot0_fac_mor_spec_1'
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
    
    vars = [e_var,b_var,pf_var]
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
;    ztitle = 'ASI Count (#)'
;    zstep = zrange[1]/2
;    ztickv = make_bins(zrange,zstep)
;    zticks = n_elements(ztickv)-1
;    zticklen = yticklen_chsz*xchsz/(cbpos[2]-cbpos[0])
;    zminor = 5
    options, var, 'zrange', zrange
    options, var, 'no_color_scale', 1


    ; O+ outflow.
    vars = rbsp_read_en_spec_combo(time_range, probe=probe, species='o')
    o_var = vars['anti']
    get_data, o_var, times, data, vals
    index = where(data eq 0 or finite(data,nan=1), count)
    if count ne 0 then begin
        data[index] = 1e-2
        store_data, o_var, times, data, vals
    endif
    var = o_var
    yrange = [5,5e4]
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
    ;ytickn[1:*:2] = ' '
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor
    options, var, 'ytickname', ytickn
    options, var, 'constant', ytickv
;    options, var, 'color_table', 40
;    options, var, 'color_table', 64
    options, var, 'color_table', 64
    options, var, 'ytitle', 'Energy!C(eV)'
    options, var, 'zrange', [1e4,1e6]/2
    options, var, 'zticklen', -0.3
    options, var, 'zminor', 10
    
    
;---Plot.
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,'fig_2013_0501_fig_overview_v03.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz

    plot_vars = [asi_keo_var,asi_count_var,pf_var,pf_spec_var,o_var]
    nplot_var = n_elements(plot_vars)
    nother_var = 4
    nvar = nplot_var+nother_var
    fig_letters = letters(nvar)
    fig_labels = fig_letters+') '+['Keo','Aurora','FAC S','S!D||!N spec','O+!Coutflow',strarr(nother_var)+'']

    pos_xrange = [2,-10]
    pos_yrange = [-1,1]*2.5
    pos_xpan_size = 3d
    pos_aspect_ratio = abs(total(pos_yrange*[-1,1])/total(pos_xrange*[-1,1]))
    pos_ypan_size = pos_aspect_ratio*pos_xpan_size
    
    asi_aspect_ratio = 0.6
    asi_ypan_size = asi_aspect_ratio*pos_xpan_size
    
    left_ypan_size = (pos_ypan_size+asi_ypan_size+abs_ychsz*4)/3
    
    
    xpans = [pos_xpan_size,pos_xpan_size]
    ypans = fltarr(nplot_var)+1
    ypads = 0.4

    margins = [10,4,2,1]
    label_size = 0.8

    tmp = panel_pos(plot_file, pansize=[pos_xpan_size,left_ypan_size], $
        xpans=xpans, ypans=ypans, ypads=ypads, margins=margins, $
        fig_size=fig_size)

    poss_left = reform(tmp[*,0,*])
    poss_right = reform(tmp[*,1,*])

    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    

;---Field line.
    fig_label = fig_labels[nplot_var+1]
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

    line_mlats = [60,65,70,75,80]
    line_mlats = [60,62,64,66,68,70]
    nline = n_elements(line_mlats)
    lines = list()
    ps = geopack_recalc(the_time)
    label_size = 0.8
    if n_elements(h0) eq 0 then h0 = 100d
    r0 = h0/constant('re')+1
    rad = constant('rad')
    deg = constant('deg')
    
    r_gsm_var = prefix+'r_gsm'
    
    r_sm = transpose(get_var_data(r_gsm_var, at=the_time))
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

;    pansize = abs([total(xr*[-1,1]),total(yr*[-1,1])])
;    pansize = pansize/pansize[0]*4


    xticklen_chsz = -0.15   ; in ychsz.
    yticklen_chsz = -0.30   ; in xchsz.
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

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
    plots, the_r_sm[0], the_r_sm[2], psym=6, color=sgcolor('red'), symsize=label_size
    tmp = convert_coord(the_r_sm[0],the_r_sm[2], data=1, to_normal=1)
    tx = tmp[0]+xchsz*0.5
    ty = tmp[1]
    msg = 'RBSP-'+strupcase(probe)
    xyouts, tx,ty,normal=1, msg, charsize=label_size, color=sgcolor('red')

    foreach dir, [-1,1], ii do begin
        geopack_trace, the_r_gsm[0],the_r_gsm[1],the_r_gsm[2], dir, par, xf,yf,zf, fline=fline, $
            t89=t89, t96=t96, t01=t01, ts04=t04s, storm=storm, igrf=igrf, r0=r0
        ntime = n_elements(fline)/3
        fline = cotran(fline, fltarr(ntime)+the_time, 'gsm2sm')
        srotate, fline, (24-the_mlt)*15*rad, 2
        oplot, fline[*,0], fline[*,2], linestyle=0, color=fline_colors[ii], thick=fline_thick
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
    ty = tpos[3]-ychsz*0.8
    xyouts, tx,ty,fig_label, normal=1


    

;---tplot panels.
    plot_poss = poss_left

    for ii=0,nplot_var-1 do begin
        tpos = plot_poss[*,ii]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        options, plot_vars[ii], 'xticklen', xticklen
        options, plot_vars[ii], 'yticklen', yticklen
    endfor

    tplot, plot_vars, trange=time_range, position=plot_poss, noerase=1
    for ii=0,nplot_var-1 do begin
        tpos = plot_poss[*,ii]
        fig_label = fig_labels[ii]
        tx = 2*xchsz
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,fig_label, normal=1
    endfor
    
    
    ; Add labels for pflux.
    pid = where(plot_vars eq pf_var, count)
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
        msg = 'Norm to 110 km altitude!CFiltered in '+$
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
    
    
    ; Add labels for o_var.
    pid = where(plot_vars eq o_var, count)
    if count ne 0 then begin
        tpos = plot_poss[*,pid]
        xrange = time_range
        yrange = get_setting(o_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            nodata=1, noerase=1, position=tpos
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = 'PA [135,180] deg'
        xyouts, tx,ty,normal=1, msg, charsize=label_size
    endif

    ; Add label for keogram.
    pid = where(plot_vars eq asi_keo_var, count)
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
        oplot, times, the_fmlat, linestyle=2, color=sgcolor('red')
        tx = times[0]
        ty = the_fmlat[0]
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.3
        msg = 'RBSP-'+strupcase(probe)
        xyouts, tx,ty, normal=1, msg, charsize=label_size, color=sgcolor('red')
    endif
    
    
    ; Add label for ASI count.
    pid = where(plot_vars eq asi_count_var, count)
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
    ;    msg = tmp+tex2str('times')+tmp+' deg around!Cfootpoint ('+string(external_model)+')'
        xyouts, tx,ty,normal=1, msg, charsize=label_size
    endif
    
    ; Add vertical grids and times.
    bar_times = make_bins(time_range,600, inner=1)
    timebar, bar_times, linestyle=1, color=sgcolor('black')
    timebar, snapshot_time, color=sgcolor('red'), linestyle=3

    

;---ASI snapshot.
    tpos = poss_right[*,0]
    fig_label = fig_labels[nplot_var]
    
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
    zticklen = yticklen_chsz*xchsz/(cbpos[2]-cbpos[0])
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

    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

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
    color = sgcolor('red')
    fmlt = fmlt[model_index]
    fmlat = fmlat[model_index]

    tmp = convert_coord(fmlt, fmlat, /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    plots, tx,ty,normal=1, psym=6, symsize=0.5, color=color
    msg = 'RBSP-B'
    xyouts, tx,ty+ychsz*0.5,normal=1, msg, color=color, alignment=0.5

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

    
    
;---E/B ratio.
    tpos = poss_right[*,3]
    tpos[0] += xchsz*8
    tpos[2] -= xchsz*12
    fig_label = fig_letters[nplot_var+2]+')'

    ebr_var = prefix+'ebr_component'
    get_data, ebr_var, ps, ebr
    fs = 1d3/ps
    xrange = [100,1e4]
    log_xrange = alog10(xrange)
    log_xrange = [ceil(log_xrange),floor(log_xrange)]
    log_xtickv = make_bins(log_xrange, 1)
    xtickv = 10d^log_xtickv
    xticks = n_elements(xtickv)-1
    xminor = 10
    xtitle = 'E/B ratio (km/s)'
    xtickn = string(xtickv,format='(I0)')
    xtickn[0] = xtitle
    xtitle = ' '
    xlog = 1
    yrange = get_setting(pf_spec_var, 'yrange')
    ytickv = get_setting(pf_spec_var, 'ytickv')
    yticks = get_setting(pf_spec_var, 'yticks')
    ytickn = get_setting(pf_spec_var, 'ytickname')
    yminor = get_setting(pf_spec_var, 'yminor')
    constant = get_setting(pf_spec_var, 'constant')
    ylog = 1
    ytitle = get_setting(pf_spec_var, 'ytitle')

    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

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
    
    ; Add Va.
    va = plasma_param['va']
    plots, va+[0,0], yrange, color=sgcolor('red')
    tx = va
    ty = yrange[0]
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]-ychsz*label_size
    msg = 'v!DA!N'
    xyouts, tx,ty,normal=1, msg, alignment=0.5, color=sgcolor('red')

;    ; Add E/B ratio.
;    suf = ''
;    va = plasma_param['va'+suf]
;    f_gi = plasma_param['fg'+suf]
;    vi = plasma_param['vi'+suf]
;    vf_fac = plasma_param['vf_fac']
;    vf = 100
;    ebr_theory = va*sqrt(1+(1/ps/f_gi*(vi/vf))^2)
;    plots, ebr_theory, fs, color=sgcolor('red')


    ; Add labeling.
    tx = tpos[0]-xchsz*3
    ty = tpos[3]-ychsz*0.8
    msg = fig_label
    xyouts, tx,ty,normal=1, msg


;---Add O+ PA.
    ; Adopted from plot_hope_l3_pitch2d_polygon.
    pa_pos = poss_right[*,-1]
    pa_pos[3] -= ychsz*2.5
    pa_pos[0] += xchsz*8
    pa_times = time_double(['2013-05-01/07:40:53','2013-05-01/07:41:16'])
;    pa_times = time_double(['2013-05-01/07:41:16','2013-05-01/07:43:09'])
;   07:40:53, 07:41:16, 07:42:01, 07:42:24, 07:42:47, 07:43:09
    unit = 'energy'
    log = 1
    mass0 = 1*(1.67e-27/1.6e-19)
    zrng = [4,6]
    
    pa_epochs = stoepoch(pa_times,'unix')
    npa_time = n_elements(pa_times)
    pa_poss = sgcalcpos(1,npa_time,position=pa_pos, xpad=0.5)
    hopel3 = sread_rbsp_hope_l3(time_range, probe=probe)
    hope_epochs = hopel3.Epoch_Ion
    hope_dtimes = hopel3.Epoch_Ion_DELTA*1e-3

    the_letters = letters(npa_time)
    for pa_id=0,npa_time-1 do begin
        tpos = pa_poss[*,pa_id]
        dy = (tpos[3]-tpos[1])*fig_size[1]
        dx = (tpos[2]-tpos[0])*fig_size[0]
        dsize = min([dx,dy])
        cx = (tpos[0]+tpos[2])*0.5
        cy = (tpos[1]+tpos[3])*0.5
        dx = dsize/fig_size[0]*0.5
        dy = dsize/fig_size[1]*0.5
        tpos = [cx-dx,cy-dy,cx+dx,cy+dy]

        tmp = min(hope_epochs-pa_epochs[pa_id], rec, abs=1)
        ut0 = sfmepoch(hope_epochs[rec],'unix')
        type_dt = 2*hope_dtimes[rec]
        type = 'oxygen'
        
        datl3 = reform((hopel3.FODU)[rec,*,*])
        enl3s = reform((hopel3.HOPE_ENERGY_Ion)[rec,*])
        pal3s = reform((hopel3.pitch_angle))
        npal3 = n_elements(pal3s)
        idx = where(datl3 eq -1e31, cnt)
        if cnt ne 0 then datl3[idx] = !values.d_nan

        idx = uniq(enl3s,sort(enl3s))
        enl3s = enl3s[idx]
        nenl3 = n_elements(enl3s)
        datl3 = datl3[idx,*]

        tdat = transpose([[datl3]])     ; in [npa,nen].
        tang = [pal3s]
        case unit of
            'energy': begin
                tdis = enl3s
                xtitl = 'E (eV)'
                end
            'velocity': begin
                tdis = sqrt(2*enl3s/mass0)*1e-3
                xtitl = 'V (km/s)'
                end
        endcase
        if keyword_set(log) then begin
            tdis = alog10(tdis)
            xtitl = 'Log!D10!N '+xtitl
        endif

        tang = tang*rad
        ndis = n_elements(tdis)
        cangs = [0,9,27,45,63,81,99,117,135,153,171,180]*!dpi/180
        cdiss = 0.5*(tdis[0:ndis-2]+tdis[1:ndis-1])
        cdiss = [cdiss[0]*2-cdiss[1],cdiss,cdiss[ndis-2]*2-cdiss[ndis-3]]
        angs = tang
        diss = tdis

        ; remove nan.
        idx = where(finite(tdat,/nan))
        tdat[idx] = 0
        
        ; Plot the PA distribution.
        nztick = 10
        ct = get_setting(o_var, 'color_table')
        
        zttl = 'Log Eflux (s!E-1!Ncm!E-2!Nsr!E-1!NeV!E-1!N)'
        if n_elements(zrng) eq 0 then zrng = sg_autolim(alog10(tdat))
        if n_elements(nztk) eq 0 then nztk = 20
        ztks = 10^smkarthm(zrng[0],zrng[1],nztk,'n')
        zcls = floor(smkarthm(10,250,nztk,'n'))

        if n_elements(xttl0) eq 0 then xttl0 = 'Log!D10!NE (eV)'
        xttl = 'Para '+xttl0
        xrng = [-1,1]*max(diss,/nan)

        yttl = 'Perp '+xttl0
        yrng = xrng
        if pa_id ne 0 then yttl = ' '
        if pa_id ne 0 then ytickformat = '(A1)' else ytickformat = ''

        ; setup coord.
        plot, xrng, yrng, position = tpos, /normal, /noerase, /nodata, $
            xstyle = 5, xrange = xrng, $
            ystyle = 5, yrange = yrng, iso=1

        sz = size(tdat,/dimensions)
        nang = sz[0]
        ndis = sz[1]

        tmp = findgen(11)*2*!dpi/10
        txs = cos(tmp)
        tys = sin(tmp)

        for i = 0, ndis-1 do begin
            for j = 0, nang-1 do begin
                if tdat[j,i] eq 0 then continue

                tmp = where(ztks le tdat[j,i], cnt)
                if cnt eq 0 then idx = 0 else idx = tmp[cnt-1]
                ;            tc = sgcolor(zcls[idx],ct=ct,file='ct2')
                tc = sgcolor(zcls[idx],ct=ct)

                tcdis = cdiss[[i,i,i+1,i+1,j]]
                tcang = cangs[[j,j+1,j+1,j,j]]

                tx = tcdis*cos(tcang)
                ty = tcdis*sin(tcang)
                polyfill, tx, ty, /data, color = tc

                tx = tcdis*cos(2*!dpi-tcang)
                ty = tcdis*sin(2*!dpi-tcang)
                polyfill, tx, ty, /data, color = tc

            endfor
        endfor
        
        step = 2
        tickv = make_bins(xrng, step, inner=1)
        ticks = n_elements(tickv)-1
        minor = step
        

        ; plot the box.
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        plot, xrng, yrng, /normal, /noerase, /nodata, $
            xstyle = 1, xrange = xrng, xtitle = xttl, $
            ystyle = 1, yrange = yrng, ytitle = yttl, color = blck, $
            xticklen=xticklen, yticklen=yticklen, iso=1, $
            ytickformat=ytickformat, $
            xtickv=tickv, xticks=ticks, xminor=minor, $
            ytickv=tickv, yticks=ticks, yminor=minor, $
            position = tpos
        
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*1
        the_letter = string(pa_id+1,format='(I0)')
        msg = 'e-'+the_letter+')'
        xyouts, tx,ty,normal=1, msg, alignment=1
            

        ; Add a label to the o_var panel.
        pid = where(plot_vars eq o_var)
        tpos = plot_poss[*,pid]
        xrange = time_range
        yrange = get_setting(o_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            nodata=1, noerase=1, position=tpos
        tx = ut0
        ty = 5e2
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]
        plots, tx+[0,0], ty+[1,-1]*ychsz*0.25, normal=1, thick=fline_thick
        xyouts, tx,ty-ychsz*1, the_letter, normal=1, alignment=0.5, charsize=label_size
    endfor
        
    
    if keyword_set(test) then stop
    sgclose
    
    

end

print, fig_2013_0501_overview_v03(event_info=event_info)
end