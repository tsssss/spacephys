;+
; mission_probe. [mission,probe]
; bin_type=. 'sm_xyz','mlat_vars'.
;-

function micro_injection_spatial_distribution_all_times_v01_mlat_vars, mission_probe, test=test

    bin_type = 'mlat_vars'
    count_zrange = [1,1e4]
    count_zrange = [1,1e3]

    mlt_bin_range = [0,24d]
    mlt_bin_step = 0.25
    mlat_bin_range = [-1,1]*50d
    mlat_bin_step = 2
    dis_bin_range = [0,15]
    dis_bin_step = 0.5

    mlt_bins = make_bins(mlt_bin_range, mlt_bin_step)
    mlat_bins = make_bins(mlat_bin_range, mlat_bin_step)
    dis_bins = make_bins(dis_bin_range, dis_bin_step)
    nmlt_bin = n_elements(mlt_bins)
    nmlat_bin = n_elements(mlat_bins)
    ndis_bin = n_elements(dis_bins)

    project_info = micro_injection_stat_load_project()
    common_time_step = project_info['common_time_step']
    plot_dir = project_info['plot_dir']
    data_dir = join_path([project_info['data_dir'],'misc'])
    mission = mission_probe[0]
    probe = mission_probe[1]
    mission_probe_str = mission+'_'+probe
    file = get_filename()
    version = get_file_version(file)

    plot_file = join_path([plot_dir,'micro_injection_spatial_distribution_'+mission_probe_str+'_all_times_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    if keyword_set(test) then magn = 2 else magn = 1
    margins = [8,4,2,4]
    margins = [3,1,8,1]
    
    ct = 49
    top_color = 254
    rad = constant('rad')
    if keyword_set(test) then begin
        hsize = 15
        thick = 2
    endif else begin
        hsize = 60
        thick = 4
    endelse


    plot_xrange = [-6,13]
    plot_yrange = [-1,1]*13
    plot_xrange = [15,-6]
    plot_yrange = [1,-1]*13
    plot_zrange = [-9,2]
    pansize = [total(plot_xrange*[-1,1]),total(plot_yrange*[-1,1])]
    pansize = abs(pansize/pansize[1]*3)
    xpans = abs([total(plot_xrange*[-1,1]),total(plot_zrange*[-1,1])])
    poss = panel_pos(plot_file, fig_size=fig_size, pansize=pansize, xpans=xpans, nypan=1, xpad=1, margins=margins)

    pansize = [total(plot_xrange*[-1,1]),total(plot_yrange*[-1,1])]
    pansize = abs(pansize/pansize[1]*2)
    ypans = abs([total(plot_yrange*[-1,1]),total(plot_zrange*[-1,1])])
    poss = panel_pos(plot_file, fig_size=fig_size, pansize=pansize, ypans=ypans, nxpan=1, ypad=0.4, margins=margins)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz, magn=magn
    abs_ticklen = -0.3*ychsz*fig_size[1]

;---Collect counts.
    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = 'mi_'+mission_probe_str+'_'

    search_trs = micro_injection_load_search_time_range()
    nsearch_tr = n_elements(search_trs[*,0])
    all_time_var = prefix+'spatial_bin_'+bin_type+'_all_times'
;    all_time_file = join_path([data_dir,all_time_var+'_v01.cdf'])
;    if file_test(all_time_file) eq 0 then begin
        if check_if_update(all_time_var) then begin
            all_counts = fltarr(nmlt_bin, nmlat_bin, ndis_bin)

            for ii=0,nsearch_tr-1 do begin
                tr = reform(search_trs[ii,*])
                common_times = make_bins(tr+[0,-common_time_step], common_time_step)+common_time_step*0.5

                r_var = lets_read('orbit', tr+[-1,1]*common_time_step, source=mission_probe, coord='sm')
                mlat_vars = lets_read_mlat_vars(orbit_var=r_var)

                mlt_var = mlat_vars.mlt
                mlts = get_var_data(mlt_var, at=common_times)
                index = where(mlts lt 0, count)
                if count ne 0 then mlts[index] += 24
                mlt_index = floor((mlts-min(mlt_bin_range))/mlt_bin_step)
                mlat_var = mlat_vars.mlat
                mlats = get_var_data(mlat_var, at=common_times)
                mlat_index = floor((mlats-min(mlat_bin_range))/mlat_bin_step)
                dis_var = mlat_vars.dis
                diss = get_var_data(dis_var, at=common_times)
                dis_index = floor((diss-min(dis_bin_range))/dis_bin_step)

                foreach time, common_times, tid do begin
                    all_counts[mlt_index[tid],mlat_index[tid],dis_index[tid]] += 1
                endforeach
            endfor

            all_times = all_counts*common_time_step/60d
            all_time_var = var_store(all_time_var, all_times, 0)
        endif
        all_times = get_var_data(all_time_var)
;        coord_vars = ['mlt','mlat','dis']
;        cdf_save_var, all_time_var, value=all_times, filename=all_time_file, save_as_one=1, settings=dictionary('unit','min','text','duration in minutes for within each bin', 'dimensions', coord_vars)
;        cdf_save_var, all_time_var+'_mlt_bins', value=mlt_bins, filename=all_time_file, save_as_one=1, settings=dictionary('unit','h','text','mlt_bins','range',mlt_bin_range)
;        cdf_save_var, all_time_var+'_mlat_bins', value=mlat_bins, filename=all_time_file, save_as_one=1, settings=dictionary('unit','deg','text','mlat_bins','range',mlat_bin_range)
;        cdf_save_var, all_time_var+'_dis_bins', value=dis_bins, filename=all_time_file, save_as_one=1, settings=dictionary('unit','Re','text','dis_bins','range',dis_bin_range)
;        cdf_save_setting, dictionary('bin_type',bin_type), filename=all_time_file
;        lprmsg, 'Saved '+all_time_var+' to '+all_time_file+' ...'
;    endif
;    
;    all_times = cdf_read_var(all_time_var, filename=all_time_file)
;    mlt_bins = cdf_read_var(all_time_var+'_mlt_bins', filename=all_time_file)
;    mlat_bins = cdf_read_var(all_time_var+'_mlat_bins', filename=all_time_file)
;    dis_bins = cdf_read_var(all_time_var+'_dis_bins', filename=all_time_file)

;---XY plane.
    pid = 0
    tpos = poss[*,pid]
    angle_bins = mlt_bins*15*rad+!dpi
    zzs = total(all_times,2)
    if keyword_set(zlog) then begin
        zrange = alog10(count_zrange)
        zzs = bytscl(alog10(zzs), min=zrange[0], max=zrange[1], top=top_color)
    endif else begin
        zrange = count_zrange
        zzs = bytscl(zzs, min=zrange[0], max=zrange[1], top=top_color)
    endelse


    xtitle = 'X (Re)'
    ytitle = 'Y (Re)'
    xrange = plot_xrange
    yrange = plot_yrange
    xstep = 5
    xtickv = make_bins(xrange, xstep, inner=1)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    ystep = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = ystep
    xticklen = abs_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = abs_ticklen/(tpos[2]-tpos[0])/fig_size[0]
    if pid ne n_elements(ypans)-1 then begin
        xtitle = ' '
        xtickformat = '(A1)'
    endif else begin
        xtickformat = ''
    endelse
    ; Set up coord.
    plot, xrange, yrange, nodata=1, noerase=1, $
        xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
        position=tpos, iso=1
    
    ; Draw bins.
    circle_index = [0,1,3,2,0]
    for ii=0,nmlt_bin-2 do begin
        angle_range = angle_bins[ii:ii+1]
        for jj=0,ndis_bin-2 do begin
            dis_range = dis_bins[jj:jj+1]
            zz = zzs[ii,jj]
            xxs = dis_range # transpose(cos(angle_range))
            yys = dis_range # transpose(sin(angle_range))
            if min(xxs) lt min(xrange) then continue
            if max(xxs) gt max(xrange) then continue
            if min(yys) lt min(yrange) then continue
            if max(yys) gt max(yrange) then continue
            polyfill, xxs[circle_index], yys[circle_index], data=1, color=sgcolor(zz, ct=ct)
        endfor
    endfor


;    ; Draw axis.
;    plot, xrange, yrange, nodata=1, noerase=1, $
;        xrange=xrange, xstyle=1, xtitle=xtitle, $
;        xtickv=xtickv, xticks=xticks, xminor=xminor, $
;        yrange=yrange, ystyle=1, ytitle=ytitle, $
;        ytickv=ytickv, yticks=yticks, yminor=yminor, $
;        xtickformat=xtickformat, $
;        position=tpos, iso=1, $
;        xticklen=xticklen, yticklen=yticklen
    

    ; Add circles and lines.
    line_color = sgcolor('silver')
    linestyle = 2
    tmp = smkarthm(0,2*!dpi,40,'n')
    foreach rr, [5,10,15] do begin
        txs = rr*cos(tmp)
        tys = rr*sin(tmp)
        oplot, txs,tys, linestyle=linestyle, color=line_color
    endforeach
    rr = [1,max(dis_bin_range)]
    foreach tt, smkarthm(0,30,360/30,'x0')*rad do begin
        txs = rr*cos(tt)
        tys = rr*sin(tt)
        oplot, txs,tys, linestyle=linestyle, color=line_color
    endforeach
    
    
    ; Add MLT labels.
    rr = max(dis_bin_range)+1.5
    foreach tt, make_bins([-30,30],30) do begin
        tx = rr*cos(tt*rad)
        ty = rr*sin(tt*rad)
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*0.3
        the_mlt = tt/15+12
        msg = string(the_mlt,format='(I02)')
        xyouts, tx,ty,msg, normal=1, alignment=0.5
    endforeach
    mlt_tts = [-60,60]
    tts = mlt_tts
    rr = max(dis_bin_range)
    tmp = smkarthm(tts[0],tts[1],0.5,'dx')*rad
    txs = rr*cos(tmp)
    tys = rr*sin(tmp)
    plots, txs,tys, data=1
    txs = txs[[1,0]]
    tys = tys[[1,0]]
    arrow, txs[0],tys[0],txs[1],tys[1], data=1, hsize=hsize, solid=1, thick=thick
    tt = min(mlt_tts)*rad
    tx = rr*cos(tt)
    ty = rr*sin(tt)
    tmp = convert_coord(tx,ty,data=1,to_normal=1)
    tx = tmp[0]+xchsz*1
    ty = tmp[1]
    msg = 'MLT (h)'
    xyouts, tx,ty,msg, normal=1
    
    
    ; Add dis labels.
    tt = 0
    dis_rrs = [1,13d]
    rrs = dis_rrs[1]+[0,1]
    txs = rrs*cos(tt)
    tys = rrs*sin(tt)
    arrow, txs[0],tys[0],txs[1],tys[1], data=1, hsize=hsize, solid=1, thick=thick
    foreach rr, [5,10] do begin
        tx = rr*cos(tt)
        ty = rr*sin(tt)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        msg = string(rr,format='(I0)')
        tx = tmp[0]
        ty = tmp[1]+ychsz*0.25
        xyouts, tx,ty, msg, normal=1, alignment=0.5
    endforeach
    rr = mean(dis_rrs)
    tx = rr*cos(tt)
    ty = rr*sin(tt)
    tmp = convert_coord(tx,ty, data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]-ychsz*1
    msg = '|R| (Re)'
    xyouts, tx,ty, msg, normal=1, alignment=0

    dis_xrange = reverse(dis_rrs)
    dis_yrange = [0,yrange[1]]
    xminor = xstep
    xtickv = make_bins(dis_xrange, xstep, inner=1)
    if dis_xrange[0] le dis_xrange[1] then xtickv = reverse(xtickv)
    xticks = n_elements(xtickv)-1
    dis_pos = tpos
    foreach tx, dis_xrange, ii do begin
        ty = dis_yrange[ii]
        tmp = convert_coord(tx,ty,data=1, to_normal=1)
        if ii eq 0 then begin
            dis_pos[0] = tmp[0]
            dis_pos[1] = tmp[1]
        endif else begin
            dis_pos[2] = tmp[0]
            dis_pos[3] = tmp[1]
        endelse
    endforeach
    set_axis, xrange=dis_xrange, yrange=yrange, $
        position=dis_pos
    xticklen = abs_ticklen/(dis_pos[3]-dis_pos[1])/fig_size[1]
    yticklen = abs_ticklen/(dis_pos[2]-dis_pos[0])/fig_size[0]
    axis, xaxis=0, $
        xrange=dis_xrange, xstyle=1, xtitle=xtitle, $
        xtickv=xtickv, xticks=xticks, xminor=xminor, $
        xtickformat=xtickformat, $
        xticklen=xticklen, yticklen=yticklen

    set_axis, xrange=xrange, yrange=yrange, $
        position=tpos
    
    tmp = lets_add_earth()


;---XZ plane.
    pid = 1
    tpos = poss[*,pid]
    angle_bins = mlat_bins*rad
    zzs = total(all_times,1)
    if keyword_set(zlog) then begin
        zrange = alog10(count_zrange)
        zzs = bytscl(alog10(zzs), min=zrange[0], max=zrange[1], top=top_color)
    endif else begin
        zrange = count_zrange
        zzs = bytscl(zzs, min=zrange[0], max=zrange[1], top=top_color)
    endelse


    xtitle = 'X (Re)'
    ytitle = 'Z (Re)'
    xrange = plot_xrange
    yrange = plot_zrange
    xstep = 5
    xtickv = make_bins(xrange, xstep, inner=1)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    ystep = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = ystep
    xticklen = abs_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = abs_ticklen/(tpos[2]-tpos[0])/fig_size[0]
    if pid ne n_elements(ypans)-1 then begin
        xtitle = ' '
        xtickformat = '(A1)'
    endif else begin
        xtickformat = ''
    endelse
    ; Set up coord.
    plot, xrange, yrange, nodata=1, noerase=1, $
        xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
        position=tpos, iso=1

    ; Draw bins.
    circle_index = [0,1,3,2,0]
    for ii=0,nmlat_bin-2 do begin
        angle_range = angle_bins[ii:ii+1]
        for jj=0,ndis_bin-2 do begin
            dis_range = dis_bins[jj:jj+1]
            zz = zzs[ii,jj]
            xxs = dis_range # transpose(cos(angle_range))
            yys = dis_range # transpose(sin(angle_range))
            if min(xxs) lt min(xrange) then continue
            if max(xxs) gt max(xrange) then continue
            if min(yys) lt min(yrange) then continue
            if max(yys) gt max(yrange) then continue
            polyfill, xxs[circle_index], yys[circle_index], data=1, color=sgcolor(zz, ct=ct)
        endfor
    endfor


;    ; Draw axis.
;    plot, xrange, yrange, nodata=1, noerase=1, $
;        xrange=xrange, xstyle=1, xtitle=xtitle, $
;        xtickv=xtickv, xticks=xticks, xminor=xminor, $
;        yrange=yrange, ystyle=1, ytitle=ytitle, $
;        ytickv=ytickv, yticks=yticks, yminor=yminor, $
;        xtickformat=xtickformat, $
;        position=tpos, iso=1, $
;        xticklen=xticklen, yticklen=yticklen
    


    ; Add circles and lines.
    tmp = smkarthm(0,2*!dpi,40,'n')
    foreach rr, [5,10,15] do begin
        txs = rr*cos(tmp)
        tys = rr*sin(tmp)
        oplot, txs,tys, linestyle=linestyle, color=line_color
    endforeach
    rrs = [1,max(dis_bin_range)]
    foreach tt, smkarthm(0,30,360/30,'x0')*rad do begin
        txs = rrs*cos(tt)
        tys = rrs*sin(tt)
        oplot, txs,tys, linestyle=linestyle, color=line_color
    endforeach
    
    ; Add MLat labels.
    rr = max(dis_bin_range)+1.5
    foreach tt, [0,-30] do begin
        tx = rr*cos(tt*rad)
        ty = rr*sin(tt*rad)
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*0.3
        msg = string(tt,format='(I0)')
        xyouts, tx,ty,msg, normal=1, alignment=0.5
    endforeach
    mlat_tts = [-37,10]
    tts = mlat_tts
    rr = max(dis_bin_range)
    tmp = smkarthm(tts[0],tts[1],0.5,'dx')*rad
    txs = rr*cos(tmp)
    tys = rr*sin(tmp)
    plots, txs,tys, data=1
    txs = txs[-2:-1]
    tys = tys[-2:-1]
    arrow, txs[0],tys[0],txs[1],tys[1], data=1, hsize=hsize, solid=1, thick=thick
    tt = min(mlat_tts)*rad
    tx = rr*cos(tt)
    ty = rr*sin(tt)
    tmp = convert_coord(tx,ty,data=1,to_normal=1)
    tx = tmp[0]+xchsz*1
    ty = tmp[1]
    msg = 'MLat (deg)'
    xyouts, tx,ty,msg, normal=1

        
    ; Add dis labels.
    tt = 0
    dis_rrs = [1,13d]
    rrs = dis_rrs[1]+[0,1]
    txs = rrs*cos(tt)
    tys = rrs*sin(tt)
    arrow, txs[0],tys[0],txs[1],tys[1], data=1, hsize=hsize, solid=1, thick=thick
    foreach rr, [5,10] do begin
        tx = rr*cos(tt)
        ty = rr*sin(tt)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        msg = string(rr,format='(I0)')
        tx = tmp[0]
        ty = tmp[1]+ychsz*0.25
        xyouts, tx,ty, msg, normal=1, alignment=0.5
    endforeach
;    rr = mean(rrs)
;    tx = rr*cos(tt)
;    ty = rr*sin(tt)
;    tmp = convert_coord(tx,ty, data=1, to_normal=1)
;    tx = tmp[0]
;    ty = tmp[1]-ychsz*1
;    msg = '|R| (Re)'
;    xyouts, tx,ty, msg, normal=1, alignment=0
    
    dis_xrange = reverse(dis_rrs)
    dis_yrange = [0,yrange[1]]
    xminor = xstep
    xtickv = make_bins(dis_xrange, xstep, inner=1)
    if dis_xrange[0] le dis_xrange[1] then xtickv = reverse(xtickv)
    xticks = n_elements(xtickv)-1
    dis_pos = tpos
    foreach tx, dis_xrange, ii do begin
        ty = dis_yrange[ii]
        tmp = convert_coord(tx,ty,data=1, to_normal=1)
        if ii eq 0 then begin
            dis_pos[0] = tmp[0]
            dis_pos[1] = tmp[1]
        endif else begin
            dis_pos[2] = tmp[0]
            dis_pos[3] = tmp[1]
        endelse
    endforeach
    set_axis, xrange=dis_xrange, yrange=yrange, $
        position=dis_pos
    xticklen = abs_ticklen/(dis_pos[3]-dis_pos[1])/fig_size[1]
    yticklen = abs_ticklen/(dis_pos[2]-dis_pos[0])/fig_size[0]
    axis, xaxis=0, $
        xrange=dis_xrange, xstyle=1, xtitle='', $
        xtickv=xtickv, xticks=xticks, xminor=xminor, $
        xtickformat='(A1)', $
        xticklen=xticklen, yticklen=yticklen

    set_axis, xrange=xrange, yrange=yrange, $
        position=tpos

    tmp = lets_add_earth()


;---color bar.
    cbpos = poss[*,0]
    cbpos[1] = poss[1,-1]
    cbpos[0] = cbpos[2]+xchsz*0.8
    cbpos[2] = cbpos[0]+xchsz*0.8
    cb_hor = 0
    ztitle = 'Time (min)'
    zticklen = abs_ticklen/(cbpos[2]-cbpos[0])/fig_size[0]
    sgcolorbar, findgen(top_color), ct=ct, position=cbpos, horizontal=cb_hor, ztitle=ztitle, $
        zrange=count_zrange, zticklen=zticklen, zcharsize=0.9
    
    
;---Add labels.
    nvar = n_elements(xpans)
    fig_labels = letters([0,2]+2)+'.'
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        tx = 1*xchsz
        ty = tpos[3]+ychsz*0.1
        msg = fig_labels[ii]
        xyouts, tx,ty,msg, normal=1, alignment=0
    endfor


    if keyword_set(test) then stop
    sgclose


end

function micro_injection_spatial_distribution_all_times_v01_sm_xyz, mission_probe, test=test

    bin_type = 'sm_xyz'
    count_zrange = [0,1e3]
    
    x_bin_range = [1,-1]*15
    y_bin_range = [1,-1]*15
    z_bin_range = [-1,1]*15
    xstep = 0.5
    ystep = 0.5
    zstep = 0.5
    x_bins = make_bins(minmax(x_bin_range), xstep)
    y_bins = make_bins(minmax(y_bin_range), ystep)
    z_bins = make_bins(minmax(z_bin_range), zstep)
    nxbins = n_elements(x_bins)
    nybins = n_elements(y_bins)
    nzbins = n_elements(z_bins)
    
    project_info = micro_injection_stat_load_project()
    common_time_step = project_info['common_time_step']
    plot_dir = project_info['plot_dir']
    data_dir = join_path([project_info['data_dir'],'misc'])
    mission = mission_probe[0]
    probe = mission_probe[1]
    mission_probe_str = mission+'_'+probe
    file = get_filename()
    version = get_file_version(file)
    
    plot_file = join_path([plot_dir,'micro_injection_spatial_distribution_'+mission_probe_str+'_all_times_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    if keyword_set(test) then magn = 2 else magn = 1
    margins = [3,1,6,1]
    plot_xrange = [13,-6]
    plot_yrange = [-1,1]*13
    plot_zrange = [-9,2]
    pansize = [total(plot_xrange*[-1,1]),total(plot_yrange*[-1,1])]
    pansize = abs(pansize/pansize[1]*3)
    xpans = abs([total(plot_xrange*[-1,1]),total(plot_zrange*[-1,1])])
    poss = panel_pos(plot_file, fig_size=fig_size, pansize=pansize, xpans=xpans, nypan=1, xpad=1, margins=margins)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz, magn=magn
    abs_ticklen = -0.3*ychsz*fig_size[1]

;---Collect counts.
    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = 'mi_'+mission_probe_str+'_'

    search_trs = micro_injection_load_search_time_range()
    nsearch_tr = n_elements(search_trs[*,0])
    dates = list()
    secofday = constant('secofday')
    for ii=0,nsearch_tr-1 do begin
        tr = reform(search_trs[ii,*])
        nday = total(tr*[-1,1])/secofday
        dates.add, findgen(nday)*secofday+tr[0], extract=1
    endfor

    all_time_var = prefix+'spatial_bin_'+bin_type+'_all_times'
    all_time_file = join_path([data_dir,all_time_var+'_v01.cdf'])
    if file_test(all_time_file) eq 0 then begin
        if check_if_update(all_time_var) then begin
            all_counts = fltarr(nxbins,nybins,nzbins)

            foreach date, dates do begin
                tr = date+[0,secofday] ; to avoid overlap.
                common_times = make_bins(tr+[0,-common_time_step], common_time_step)

                r_var = lets_read('orbit', tr+[-1,1]*common_time_step, source=mission_probe, coord='sm')
                r_sm = get_var_data(r_var, at=common_times)
                xs = r_sm[*,0]
                x_index = floor((xs-min(x_bin_range))/xstep)
                ys = r_sm[*,1]
                y_index = floor((ys-min(y_bin_range))/ystep)
                zs = r_sm[*,2]
                z_index = floor((zs-min(z_bin_range))/zstep)

                foreach time, common_times, tid do begin
                    all_counts[x_index[tid],y_index[tid],z_index[tid]] += 1
                endforeach
            endforeach

            all_times = all_counts*common_time_step/60d
            all_time_var = var_store(all_time_var, all_counts, 0)
        endif
        all_times = get_var_data(all_time_var)
        coord_vars = ['x','y','z']
        cdf_save_var, all_time_var, value=all_counts, filename=all_time_file, save_as_one=1, settings=dictionary('unit','min','text','duration in minutes for within each bin', 'dimensions', coord_vars)
        cdf_save_var, all_time_var+'_x_bins', value=x_bins, filename=all_time_file, save_as_one=1, settings=dictionary('unit','Re','text','x_bins','range',x_bin_range)
        cdf_save_var, all_time_var+'_y_bins', value=y_bins, filename=all_time_file, save_as_one=1, settings=dictionary('unit','Re','text','y_bins','range',y_bin_range)
        cdf_save_var, all_time_var+'_z_bins', value=z_bins, filename=all_time_file, save_as_one=1, settings=dictionary('unit','Re','text','z_bins','range',z_bin_range)
        cdf_save_setting, dictionary('bin_type',bin_type), filename=all_time_file
        lprmsg, 'Saved '+all_time_var+' to '+all_time_file+' ...'
    endif
    
    all_times = cdf_read_var(all_time_var, filename=all_time_var)
    x_bins = cdf_read_var(all_time_var+'_x_bins', filename=all_time_file)
    y_bins = cdf_read_var(all_time_var+'_y_bins', filename=all_time_file)
    z_bins = cdf_read_var(all_time_var+'_z_bins', filename=all_time_file)

    
;---Gen plot.
    ct = 49
    top_color = 254
    x_index = where_pro(x_bins, '[]', minmax(plot_xrange), count=nxbin)
    y_index = where_pro(y_bins, '[]', minmax(plot_yrange), count=nybin)
    z_index = where_pro(z_bins, '[]', minmax(plot_zrange), count=nzbin)
    all_counts = all_counts[x_index,*,*]
    all_counts = all_counts[*,y_index,*]
    all_counts = all_counts[*,*,z_index]
    x_bins = x_bins[x_index]
    y_bins = y_bins[y_index]
    z_bins = z_bins[z_index]
    
    
    
;---XY plane.
    tpos = poss[*,0]
    the_counts = total(all_counts, 3)
    if x_bin_range[0] gt x_bin_range[1] then the_counts = reverse(the_counts, 1)
    if y_bin_range[0] gt y_bin_range[1] then the_counts = reverse(the_counts, 2)
    zzs = bytscl(the_counts, min=count_zrange[0], max=count_zrange[1], top=top_color)
    sgtv, zzs, position=tpos, ct=ct
        
    xtitle = 'SM X (Re)'
    ytitle = 'SM Y (Re)'
    xrange = plot_xrange
    yrange = plot_yrange
    xstep = 5
    xtickv = make_bins(xrange, xstep, inner=1)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    ystep = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = ystep
    xticklen = abs_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = abs_ticklen/(tpos[2]-tpos[0])/fig_size[0]
    plot, xrange, yrange, nodata=1, noerase=1, $
        xrange=xrange, xstyle=1, xtitle=xtitle, $
        xtickv=xtickv, xticks=xticks, xminor=xminor, $
        yrange=yrange, ystyle=1, ytitle=ytitle, $
        ytickv=ytickv, yticks=yticks, yminor=yminor, $
        position=tpos, iso=1, $
        xticklen=xticklen, yticklen=yticklen
    
    line_color = sgcolor('silver')
    linestyle = 0
    tmp = smkarthm(0,2*!dpi,40,'n')
    foreach rr, [5,10,15] do begin
        txs = rr*cos(tmp)
        tys = rr*sin(tmp)
        oplot, txs,tys, linestyle=linestyle, color=line_color
    endforeach
    plots, [0,0], yrange, linestyle=linestyle, color=line_color
    plots, xrange, [0,0], linestyle=linestyle, color=line_color
    
    tmp = lets_add_earth()


    
;---YZ plane.
    tpos = poss[*,1]
    the_counts = total(all_counts, 1)
    the_counts = transpose(the_counts)  ; switch to x for Y and y for Z.
    if yrange[0] gt yrange[1] then the_counts = reverse(the_counts, 2)
    zzs = bytscl(the_counts, min=count_zrange[0], max=count_zrange[1], top=top_color)
    sgtv, zzs, position=tpos, ct=ct

    xtitle = 'SM Z (Re)'
    ytitle = ' '
    ytickformat = '(A1)'
    xrange = plot_zrange
    yrange = plot_yrange
    xstep = 5
    xtickv = make_bins(xrange, xstep, inner=1)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    ystep = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = ystep
    xticklen = abs_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = abs_ticklen/(tpos[2]-tpos[0])/fig_size[0]
    plot, xrange, yrange, nodata=1, noerase=1, $
        xrange=xrange, xstyle=1, xtitle=xtitle, $
        xtickv=xtickv, xticks=xticks, xminor=xminor, $
        yrange=yrange, ystyle=1, ytitle=ytitle, $
        ytickv=ytickv, yticks=yticks, yminor=yminor, $
        position=tpos, iso=1, $
        xticklen=xticklen, yticklen=yticklen, $
        ytickformat=ytickformat
    
    foreach tx, xtickv do begin
        if tx eq 0 then begin
            plots, tx+[0,0], [yrange[0],-1], linestyle=linestyle, color=line_color
            plots, tx+[0,0], [1,yrange[1]], linestyle=linestyle, color=line_color
        endif else begin
            plots, tx+[0,0], yrange, linestyle=linestyle, color=line_color
        endelse
    endforeach
    foreach ty, [0] do begin
        if ty eq 0 then begin
            plots, [xrange[0],-1], ty+[0,0], linestyle=linestyle, color=line_color
            plots, [1,xrange[1]], ty+[0,0], linestyle=linestyle, color=line_color
        endif else begin
            plots, xrange, ty+[0,0], linestyle=linestyle, color=line_color
        endelse
    endforeach

    tmp = lets_add_earth(only_outline=1)


;---color bar.
    cbpos = poss[*,0]
    cbpos[2] = poss[2,-1]
    cbpos[1] = cbpos[3]+ychsz*0.5
    cbpos[3] = cbpos[1]+ychsz*0.5
    ztitle = 'Time (min)'
    zticklen = abs_ticklen/(cbpos[3]-cbpos[1])/fig_size[1]
    sgcolorbar, findgen(top_color), ct=ct, position=cbpos, horizontal=1, ztitle=ztitle, $
        zrange=count_zrange, zticklen=zticklen, zcharsize=0.9

    ; Add circles and lines.
    line_color = sgcolor('silver')
    linestyle = 0
    tmp = smkarthm(0,2*!dpi,40,'n')
    foreach rr, [5,10,15] do begin
        txs = rr*cos(tmp)
        tys = rr*sin(tmp)
        oplot, txs,tys, linestyle=linestyle, color=line_color
    endforeach
    plots, [0,0], yrange, linestyle=linestyle, color=line_color
    plots, xrange, [0,0], linestyle=linestyle, color=line_color
    

    
;---Add labels.
    nvar = n_elements(xpans)
    fig_labels = letters([0,2]+2)+'.'
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        tx = tpos[0]+xchsz*1
        ty = tpos[3]-ychsz*1
        msg = fig_labels[ii]
        xyouts, tx,ty,msg, normal=1
    endfor


    if keyword_set(test) then stop
    sgclose


end


function micro_injection_spatial_distribution_all_times_v01, mission_probe, id=ids, test=test

    if n_elements(mission_probe) eq 0 then mission_probe = ['mms','1']
    ;if n_elements(ids) eq 0 then ids = 'sm_xyz'
    if n_elements(ids) eq 0 then ids = 'mlat_vars'
    my_name = get_filename()
    subroutine = ''
    foreach id, ids do begin
        subroutine += '_'+id
    endforeach
    routine = get_file_stem(get_file_base(my_name))+subroutine
    return, call_function(routine, mission_probe, test=test)

end

test = 1
print, micro_injection_spatial_distribution_all_times_v01(test=test)
end