;+
; mission_probe. [mission,probe]
;-


function micro_injection_spatial_distribution_flux_simplest_v01, mission_probe, test=test, var_type=var_type


;---Settings.
    if n_elements(mission_probe) eq 0 then mission_probe = ['mms','1']
    bin_type = 'mlat_vars'
    var_zrange = [0,3]
    var_zlog = 0
    ct = 49 ; 70: blue-red


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
    plot_dir = project_info['plot_dir']
    data_dir = join_path([project_info['data_dir'],'misc'])
    mission = mission_probe[0]
    probe = mission_probe[1]
    mission_probe_str = mission+'_'+probe
    file = get_filename()
    version = get_file_version(file)

    plot_file = join_path([plot_dir,'micro_injection_spatial_distribution_'+mission_probe_str+'_'+var_type+'_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    if keyword_set(test) then magn = 2 else magn = 1
    margins = [8,4,2,4]
    margins = [3,1,8,1]
    
    top_color = 254
    rad = constant('rad')
    if keyword_set(test) then begin
        hsize = 15
        thick = 2
    endif else begin
        hsize = 120
        thick = 2
    endelse


    plot_xrange = [-6,13]
    plot_yrange = [-1,1]*13
    plot_xrange = [15,-6]
    plot_yrange = [1,-1]*13
    plot_zrange = [-9,2]
    pansize = [total(plot_xrange*[-1,1]),total(plot_yrange*[-1,1])]
    pansize = abs(pansize/pansize[1]*2)
    ypans = abs([total(plot_yrange*[-1,1]),total(plot_zrange*[-1,1])])
    nxpan = 1
    xpad = 10
    all_poss = panel_pos(plot_file, fig_size=fig_size, pansize=pansize, ypans=ypans, nxpan=nxpan, xpad=xpad, ypad=0.4, margins=margins)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz, magn=magn
    abs_ticklen = -0.3*ychsz*fig_size[1]


;---Load data and collect wanted var.
    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = mission+probe+'_'

    project_info = micro_injection_stat_load_project()
    common_time_step = project_info['common_time_step']
    sample_energys = project_info['sample_energys']
    high_energy = max(sample_energys)
    high_energy_str = string(high_energy,format='(I0)')
    
    search_trs = micro_injection_load_search_time_range()
    nsearch_tr = n_elements(search_trs[*,0])
    binned_var = prefix+'spatial_bin_'+var_type
    if check_if_update(binned_var) then begin
        binned_data = fltarr(nmlt_bin,nmlat_bin,ndis_bin)
        binned_counts = fltarr(nmlt_bin,nmlat_bin,ndis_bin)

        for ii=0,nsearch_tr-1 do begin
            tr = reform(search_trs[ii,*])

        ;---Load MI times.
            mi_trs = micro_injection_stat_read_event_times(tr, mission_probe=mission_probe)
            nmi_tr = n_elements(mi_trs[*,0])
            mi_times = list()
            for tid=0,nmi_tr-1 do begin
                mi_tr = reform(mi_trs[tid,*])
                mi_times.add, make_bins(mi_tr,common_time_step), extract=1
            endfor
            mi_times = mi_times.toarray()

        ;---Load orbit.
            r_var = lets_read('orbit', tr+[-1,1]*common_time_step, source=mission_probe, coord='sm')
            mlat_vars = lets_read_mlat_vars(orbit_var=r_var)

        ;---Load wanted data.
            flux_vars = micro_injection_read_cwt_kev_electron(tr, mission_probe=mission_probe, log=1)
            spec_vars = flux_vars+'_cwt'
            index = where(stregex(flux_vars, var_type) ne -1, count)
            wanted_var = flux_vars[index]

        ;---Collect results.
            mlt_var = mlat_vars.mlt
            mlts = get_var_data(mlt_var, at=mi_times)
            index = where(mlts lt 0, count)
            if count ne 0 then mlts[index] += 24
            mlt_index = floor((mlts-min(mlt_bin_range))/mlt_bin_step)
            mlat_var = mlat_vars.mlat
            mlats = get_var_data(mlat_var, at=mi_times)
            mlat_index = floor((mlats-min(mlat_bin_range))/mlat_bin_step)
            dis_var = mlat_vars.dis
            diss = get_var_data(dis_var, at=mi_times)
            dis_index = floor((diss-min(dis_bin_range))/dis_bin_step)

            wanted_data_in_tr = get_var_data(wanted_var, at=mi_times)

            foreach time, mi_times, tid do begin
                binned_counts[mlt_index[tid],mlat_index[tid],dis_index[tid]] += 1
                binned_data[mlt_index[tid],mlat_index[tid],dis_index[tid]] += wanted_data_in_tr[tid]
            endforeach
        endfor

        binned_data = binned_data/binned_counts
        binned_var = var_store(binned_var, binned_data, 0)
    endif
    binned_data = get_var_data(binned_var)


;---Plot.
    poss = reform(all_poss)

;---XY plane.
    pid = 0
    tpos = poss[*,pid]
    angle_bins = mlt_bins*15*rad+!dpi
    zzs = mean(binned_data,dimension=2,nan=1)
    if keyword_set(var_zlog) then begin
        zrange = alog10(var_zrange)
        zzs = bytscl(alog10(zzs), min=zrange[0], max=zrange[1], top=top_color)
    endif else begin
        zrange = var_zrange
        zzs = bytscl(zzs, min=zrange[0], max=zrange[1], top=top_color)
    endelse


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


;    pdyn_var = prefix+'spatial_bin_pdyn'
;    if check_if_update(pdyn_var) then begin
;        pdyns = []
;        times = []
;        for ii=0,nsearch_tr-1 do begin
;            tr = reform(search_trs[ii,*])
;            pdyn_var = omni_read_sw_p(tr)
;
;            mi_trs = micro_injection_stat_read_event_times(tr, mission_probe=mission_probe)
;            nmi_tr = n_elements(mi_trs[*,0])
;            mi_times = list()
;            for tid=0,nmi_tr-1 do begin
;                mi_tr = reform(mi_trs[tid,*])
;                mi_times.add, make_bins(mi_tr,common_time_step), extract=1
;            endfor
;            mi_times = mi_times.toarray()
;
;            times = [times,mi_times]
;            pdyns = [pdyns,get_var_data(pdyn_var, at=mi_times)]
;        endfor
;        store_data, pdyn_var, times, pdyns
;    endif
;    pdyns = get_var_data(pdyn_var)
;    pdyn_range = minmax(pdyns)

    ; Add Earth and magnetosphere.
    set_axis, xrange=xrange, yrange=yrange, position=tpos
    tmp = lets_add_earth()
    pdyn_range = project_info['pdyn_range']
    tts = smkarthm(0,2*!dpi,50,'n')
    rrs = 2
    xxs = rrs*cos(tts)
    yys = rrs*cos(tts)
    zzs = fltarr(n_elements(xxs))
    color = sgcolor('red')
    foreach pdyn, pdyn_range do begin
        mpause_t96, pdyn, xmgnp=xmgnp, ymgnp=ymgnp, zmgnp=zmgnp, $
            xgsm=xxs, ygsm=yys, zgsm=zzs, id=id, distan=distan
        oplot, xmgnp, ymgnp, linestyle=2, color=color
    endforeach
    

;---XZ plane.
    pid = 1
    tpos = poss[*,pid]
    angle_bins = mlat_bins*rad
    zzs = mean(binned_data,dimension=1,nan=1)
    if keyword_set(var_zlog) then begin
        zrange = alog10(var_zrange)
        zzs = bytscl(alog10(zzs), min=zrange[0], max=zrange[1], top=top_color)
    endif else begin
        zrange = var_zrange
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

    set_axis, xrange=xrange, yrange=yrange, position=tpos
    tmp = lets_add_earth()


;---color bar.
    cbpos = poss[*,0]
    cbpos[1] = poss[1,-1]
    cbpos[0] = cbpos[2]+xchsz*0.8
    cbpos[2] = cbpos[0]+xchsz*0.8
    cb_hor = 0
    ztitle = 'Phase Lag (deg)'
    zticklen = abs_ticklen/(cbpos[2]-cbpos[0])/fig_size[0]
    sgcolorbar, findgen(top_color), ct=ct, position=cbpos, horizontal=cb_hor, ztitle=ztitle, $
        zrange=var_zrange, zticklen=zticklen, zcharsize=0.9


;---Add labels.
    nvar = n_elements(ypans)
    fig_labels = letters([0,2])+'.'
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        tx = tpos[0]-xchsz*1.5
        ty = tpos[3]+ychsz*0.1
        msg = fig_labels[ii]
        xyouts, tx,ty,msg, normal=1, alignment=0
    endfor


;---Add labels.
    nvar = n_elements(ypans)
    fig_labels = letters([0,2])+'.'
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        tx = tpos[0]-xchsz*1.5
        ty = tpos[3]+ychsz*0.1
        msg = fig_labels[ii]
        xyouts, tx,ty,msg, normal=1, alignment=0
    endfor

    if keyword_set(test) then stop
    sgclose

    return, plot_file

end


test = 0
foreach var_type, 'og_kev_electron_'+['67kev','122kev'] do begin
    print, micro_injection_spatial_distribution_flux_simplest_v01(test=test, var_type=var_type)
endforeach
end