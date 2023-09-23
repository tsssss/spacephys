;+
; DMSP data.
; Try to fix the bad data.
;-

function fig_2017_0309_0700_dmsp_overview_v01, plot_file, event_info=event_info, test=test

;---Settings.
    version = 'v01'
    id = '2017_0309_0700'
test=0

    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
    base_name = 'fig_'+id+'_dmsp_overview_'+version+'.pdf'
    ssusi_time = time_double('2017-03-09/07:32')
    asi_time = time_double('2017-03-09/07:32')
    time_range = time_double('2017-03-09/'+['07:30','07:35'])
    dmsp_orbit_time_range = time_double('2017-03-09/'+['07:30','07:35'])
    dmsp_plot_time_range = time_range     ; for tplot.
    invertedv_times = time_double(['2017-03-09/07:32'])
    invertedv_color = sgcolor('red')
    invertedv_text = 'Inverted-V'
    pos_xrange = [2,-15]
    pos_yrange = [-1,1]*3.5
    asi_setting = (event_info.ground)['asi_setting']
    asi_sites = sort_uniq(asi_setting.sites)
    tplot_options, 'version', 2
    

    ; Figure out figure size.
    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz
    asi_ypan_size = 2d
    asi_xpan_size = asi_ypan_size*1
    xpad = 5
    ypad = 1
    margins = [9,1,6,1]
    fig_ysize = asi_ypan_size*2+(total(ypad)+margins[1]+margins[3])*abs_ychsz
    pos_margins = [0,3,0,0]
    dmsp_margins = [0,2.5,3,0]
    pos_ypan_size = asi_ypan_size-(pos_margins[1]+pos_margins[3])*abs_ychsz
    pos_aspect_ratio = abs(total(pos_yrange*[-1,1])/total(pos_xrange*[-1,1]))
    pos_xpan_size = pos_ypan_size/pos_aspect_ratio
    xpans = [pos_xpan_size,asi_xpan_size]




;---Load data.
    sc_info = event_info.themis.thd
    sc_color = sgcolor('magenta')
    model_setting = sc_info['model_setting']
    internal_model = event_info['internal_model']
    external_model = event_info['external_model']
    external_model = 't89'
    models = model_setting.models
    model_index = where(models eq external_model)
    
    dmsp_info = event_info.dmsp.dmspf18
    dmsp_color = dmsp_info.sc_color

    probe = dmsp_info['probe']
    ssusi_id = 'lbhs'
;    ssusi_id = '1356'
;    ssusi_id = '1216'
    ssusi_id = 'lbhl'
    ssusi_id = 'energy'
    if ssusi_id eq '1356' then begin
        ssusi_wavelength = '135.6 nm'
    endif else if ssusi_id eq '1216' then begin
        ssusi_wavelength = '121.6 nm'
    endif else if ssusi_id eq 'energy' then begin
        ssusi_wavelength = 'energy'
    endif else ssusi_wavelength = strupcase(ssusi_id)

    dmsp_mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, id=ssusi_id)
;    mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, id='1216')
    mlat_vars = dmsp_read_mlat_vars(time_range, probe=probe, errmsg=errmsg)
    ele_spec_var = dmsp_read_en_spec(time_range, probe=probe, species='e', errmsg=errmsg)
    ion_spec_var = dmsp_read_en_spec(time_range, probe=probe, species='p', errmsg=errmsg)
    ele_eflux_var = dmsp_read_eflux(time_range, probe=probe, species='e', errmsg=errmsg)
    db_xyz_var = dmsp_read_bfield_madrigal(time_range, probe=probe)
    r_var = dmsp_read_orbit(time_range, probe=probe)

    rad = constant('rad')
    deg = constant('deg')
    mlat_var = mlat_vars[0]
    mlt_var = mlat_vars[1]
    mlon_var = mlat_vars[2]


    ; convert dB from xyz to fac.
    prefix = dmsp_info['prefix']
    fac_labels = [tex2str('perp')+','+['north','west'],tex2str('parallel')]
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


;---Plot settings.
    label_size = 0.8
    sc_label_size = 1
    tmp = smkarthm(0,2*!dpi,50,'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1

    plot_dir = event_info.plot_dir
    if n_elements(plot_file) eq 0 then plot_file = join_path([plot_dir,base_name])
    if keyword_set(test) then plot_file = 0
    poss = panel_pos(plot_file, $
        xpans=xpans,ypans=[1,1], xpad=xpad, ypad=ypad, pansize=[asi_xpan_size,asi_ypan_size], panid=[1,1], $
        fig_size=fig_size, margins=margins)

    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz



;---SC position.
    probe = sc_info['probe']
    prefix = sc_info['prefix']
    tpos = sgcalcpos(1, region=poss[*,0,0], margins=pos_margins)
    fig_label = 'a) Config'

    the_time = asi_time
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

    par_var = geopack_read_par(time_range, model=external_model, t89_par=t89_par)
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


    uniform_ticklen = -ychsz*0.15*fig_size[0]
    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    xminor = 2
    xtickv = make_bins(xr, xminor, inner=1)
    xticks = n_elements(xtickv)
    xtickn = string(xtickv,format='(I0)')
    ;    xtitle = 'SM R!DXY!N (Re)'
    xtickn[-1] = 'SM X (Re)        '
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


    ; Add SC.
    sc_p0 = the_sc_pos[[0,2]]
    plots, sc_p0[0], sc_p0[1], psym=6, color=sc_color;, symsize=label_size
    tx = tmp[0]+xchsz*1
    ty = tmp[1]-ychsz*0.3
    msg = 'TH-'+strupcase(probe)
    xyouts, tx,ty,normal=1, msg, charsize=sc_label_size, color=sc_color

    ; Add DMSP.
    dmsp_dis = 1+800d/constant('re')
    dmsp_mlat = (180+64)*constant('rad')
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
    xyouts, tx1+xchsz*0.5, ty1-ychsz*0.2, normal=1, strupcase(sc_name+' '+probe), charsize=sc_label_size, color=dmsp_color

    

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

    


;---DMSP panels.
    prefix = dmsp_info['prefix']
    dmsp_vars = prefix+['e_en_spec','e_eflux_map']
    dmsp_labels = 'b-'+['1) e-','2) Eflux']
    dmsp_vars = prefix+['e_en_spec','db_fac']
    dmsp_labels = 'b-'+['1) e-','2) dB']
    ndmsp_var = n_elements(dmsp_vars)
    
    ; map eflux and convert unit.
    get_data, prefix+'e_eflux', times, eflux, limits=lim
    var = prefix+'e_eflux_map'
    cmap = 1.4  ; this is for 800 km.
    theta_loss_cone = 45d
    sr = !dpi*sin(theta_loss_cone*constant('rad'))^2
    eflux_map = eflux*cmap*sr*1.6e-19*1e4*1e3  ; convert from eV/cm^2-s-sr to mW/m^2.
    store_data, var, times, eflux_map
    yrange = [1e-1,8e2]
    add_setting, var, smart=1, dictionary($
        'ylog', 1, $
        'yrange', yrange, $
        'display_type', 'scalar', $
        'short_name', tex2str('Gamma')+'!De!N', $
        'unit', 'mW/m!U2!N' )
    
    
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
    
    the_poss = sgcalcpos(ndmsp_var, region=poss[*,0,1], margins=dmsp_margins)
    for pid=0,ndmsp_var-1 do begin
        tpos = the_poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, dmsp_vars[pid], 'xticklen', xticklen
        options, dmsp_vars[pid], 'yticklen', yticklen
    endfor
    
    tplot, dmsp_vars, trange=dmsp_plot_time_range, position=the_poss, noerase=1

    ; Add inverted V labels.
    pid = where(dmsp_vars eq prefix+'e_en_spec')
    tpos = the_poss[*,pid]
    plot, dmsp_plot_time_range, [0,1], $
        xstyle=5, ystyle=5, position=tpos, nodata=1, noerase=1
    foreach invertedv_time, invertedv_times do begin
        tmp = convert_coord(invertedv_time,1, data=1, to_normal=1)
        tx = tmp[0]
        ty = tpos[3]+ychsz*0.4
        msg = invertedv_text
        xyouts, tx,ty,normal=1, msg, charsize=sc_label_size, alignment=0.5, color=invertedv_color
        ;tys = [ty+ychsz*1, tpos[1]+(tpos[3]-tpos[1])*0.4]
        ;plots, tx+[0,0],tys, normal=1, color=invertedv_color
        plots, tx,tpos[3],normal=1, psym=8, symsize=0.5, color=invertedv_color
    endforeach
    ;arrow, tx,tys[0], tx,tys[1], normal=1, color=invertedv_color, solid=1, hsize=arrow_hsize
    ;    timebar, invertedv_time, color=invertedv_color
    
    for pid=0,ndmsp_var-1 do begin
        tpos = the_poss[*,pid]
        tx = tpos[0]-xchsz*8
        ty = tpos[3]-ychsz*0.8
        msg = dmsp_labels[pid]
        xyouts, tx,ty,normal=1, msg
    endfor
    
    
    
    ; Add FAC.
    fac_times = list($
        ['2017-03-09/07:31:45','2017-03-09/07:32:10'])
    pid = where(dmsp_vars eq db_var, count)
    if count ne 0 then begin
        xrange = dmsp_plot_time_range
        yrange = get_setting(db_var,'yrange')
        tpos = the_poss[*,pid]
        plot, xrange, yrange, $
            xstyle=5, ystyle=5, position=tpos, nodata=1, noerase=1

        foreach tr, fac_times do begin
            tr = time_double(tr)
            the_mlat = get_var_data(mlat_var,at=tr)
            the_dis = snorm(get_var_data(r_var,at=tr))
            the_db = (get_var_data(db_var,at=tr))[*,1]
            the_color = sgcolor('green')
            txs = tr
            foreach tmp,tr, ii do begin
                tmp = convert_coord(tmp,yrange[0],data=1,to_normal=1)
                txs[ii] = tmp[0]
            endforeach
            ty = tpos[1]+(tpos[3]-tpos[1])*0.35
            plots, txs, ty+[0,0], normal=1, color=the_color
            foreach tx,txs do begin
                plots, tx+[0,0],ty+[-1,1]*ychsz*0.25, normal=1, color=the_color
            endforeach            
            
            delta_db = total(the_db*[-1,1])
            delta_mlat = total(the_mlat*[-1,1])
            delta_ll = delta_mlat*constant('rad')*mean(the_dis)*constant('re')
            mu0 = constant('mu0')
            j_para = delta_db*1e-9/(delta_ll*1e3)/mu0*1e6
            plots, tr, the_db, color=the_color, linestyle=2
            
            msg = 'Up J ~ '+string(abs(j_para),format='(F4.1)')+' '+tex2str('mu')+'A'
            xyouts, mean(txs),ty+ychsz*1, normal=1, alignment=0.5, msg, color=the_color, charsize=label_size
        endforeach
    endif
    


;---Auroral snapshots.
    ; Settings.
    mlt_range = [-1d,0]*6
    min_mlat = 50d
    mlat_range = [min_mlat,90]
    ct_ssusi = 70
    ;ct_ssusi = 49
    ssusi_zrange = [-1,1]*20
    ;ssusi_zrange = [0,1]*2
    ct_asi = 49
    asi_zlog = 1
    asi_zrange = [5e2,2e4]
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
    endif else begin
        asi_zzs = bytscl((mlt_image), min=asi_zrange[0],max=asi_zrange[1], top=color_top)
    endelse
    ;asi_zrange = [0,5e3]
    ;asi_zzs = bytscl((mlt_image), min=asi_zrange[0],max=asi_zrange[1], top=color_top)
    
    mlt_images = get_var_data(dmsp_mlt_image_var, times=times, limits=lim)
    time_id = 1
    npx = n_elements(mlt_images[0,0,*])
;    mlt_image = reform(mlt_images[time_id,*,0:npx*0.5-1])
    mlt_image = reform(mlt_images[time_id,0:npx*0.5-1,0:npx*0.5-1])
    index = where(mlt_image le -2, count)
    if count ne 0 then mlt_image[index] = 0
    npx = n_elements(mlt_image[0,*])
    ssusi_zzs = bytscl(mlt_image, min=ssusi_zrange[0], max=ssusi_zrange[1], top=color_top)
    ssusi_time_range = reform(lim.time_range[time_id,*])
    ssusi_unit = lim.unit

    aurora_info.asi = dictionary($
        'msg', 'c-1) North | white light!C'+strjoin(strupcase(asi_sites),' ')+' | '+time_string(asi_time,tformat='hh:mm')+' UT', $
        'hemisphere', 'north', $
        'position', poss[*,1,0], $
        'zzs', asi_zzs, $
        'ct', ct_asi )
    dmsp_name = dmsp_info['sc_name']
    aurora_info.ssusi = dictionary($
        'msg', 'c-2) South | '+ssusi_wavelength+'!C'+strupcase(dmsp_name+' '+dmsp_probe)+' | '+$
            strjoin(time_string(ssusi_time_range,tformat='hh:mm'),'-')+' UT', $
        'hemisphere', 'south', $
        'position', poss[*,1,1], $
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
        
        ; add label.
        label_x = -3
        label_y = 69
        if the_info.hemisphere eq 'south' then begin
            label_y = 68.5
        endif
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
        msgs = strsplit(the_info.msg,'!C',extract=1)
        xyouts, tx,ty,normal=1, msgs[1], charsize=label_size
        xyouts, tx,ty+ychsz*1,normal=1, msgs[0]
    endforeach


    ; Colorbar.
    asi_cbpos = aurora_info.asi.position
    asi_cbpos[0] = asi_cbpos[2]+xchsz*1
    asi_cbpos[2] = asi_cbpos[0]+xchsz*0.7
    asi_ztitle = 'N-hem ASI (#)'
    asi_log_ztickv = make_bins(asi_log_zrange,1, inner=1)
    asi_ztickv = 10^asi_log_ztickv
    asi_zticks = n_elements(asi_ztickv)-1
    asi_zminor = 9
    asi_ztickn = '10!U'+string(asi_log_ztickv,format='(I0)')
    asi_zlog = 1
    asi_linestyle = 1
    zticklen = uniform_ticklen/(asi_cbpos[2]-asi_cbpos[0])/fig_size[0]
    sgcolorbar, findgen(color_top), $
        ztitle=asi_ztitle, zrange=asi_zrange, ct=ct_asi, position=asi_cbpos, $
        ztickv=asi_ztickv, zticks=asi_zticks, ztickname=asi_ztickn, zminor=asi_zminor, zticklen=zticklen, log=asi_zlog


    ssusi_cbpos = aurora_info.ssusi.position
    ssusi_cbpos[0] = ssusi_cbpos[2]+xchsz*1
    ssusi_cbpos[2] = ssusi_cbpos[0]+xchsz*0.7
    ztitle = 'S-hem SSUSI '+strupcase(ssusi_id)+' ('+ssusi_unit+')'
    zrange = ssusi_zrange
    sgcolorbar, findgen(color_top), $
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
    ; trace.
    ;res = geopack_trace_to_ionosphere(r_var, models='t89', igrf=0, south=1)
    ;mlts = get_var_data(prefix+'fmlt_t89', in=time_range, times=the_times)
    ;mlats = get_var_data(prefix+'fmlat_t89', at=the_times)
    
    
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
        label_mlt = -2.5
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
    foreach the_info, aurora_info do begin
        tpos = the_info.position
        plot, [-1,0], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos

        foreach sc_info, event_info.themis do begin
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
            if probe eq 'd' then tx = tx-xchsz*2
            xyouts, tx-xchsz*0.5,ty-ychsz*1.0, alignment=0.5,normal=1, $
                strupcase(sc_name)+'-'+strupcase(probe), color=sc_color, charsize=sc_label_size
        endforeach
    endforeach

    if keyword_set(test) then stop
    sgclose
    
    return, plot_file
    

end

test = 0
print, fig_2017_0309_0700_dmsp_overview_v01(event_info=event_info, test=test)
end
