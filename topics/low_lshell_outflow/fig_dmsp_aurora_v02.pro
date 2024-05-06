;+
; Plot DMSP aurora over the storm.
;-

function fig_dmsp_aurora_v02, plot_file, event_info=event_info, test=test


    version = 'v02'
    id = '2015_0317'
    if n_elements(event_info) eq 0 then event_info = low_lshell_outflow_load_data(id)
    time_range = event_info.time_range
    dmsp_probes = 'f'+['16','17','18','19']

    mlt_type = 'night'

;---Load data.
    ssusi_id = 'energy'
    ssusi_wavelength = strupcase(ssusi_id)
    data_time_range = time_range+[-1800,0]
    foreach probe, dmsp_probes do begin
        dmsp_mlt_image_var = dmsp_read_mlt_image(data_time_range, probe=probe, id=ssusi_id, get_name=1)
        if ~check_if_update(dmsp_mlt_image_var, time_range) then continue
        dmsp_mlt_image_var = dmsp_read_mlt_image(data_time_range, probe=probe, id=ssusi_id, update=1)
        options, dmsp_mlt_image_var, 'requested_time_range', time_range
    endforeach
    
    ; Sort by time.
    common_times = []
    all_mlt_images = []
    all_probes = []
    all_time_ranges = []
    all_hems = []
    foreach probe, dmsp_probes do begin
        dmsp_mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, id=ssusi_id, get_name=1)
        mlt_images = get_var_data(dmsp_mlt_image_var, times=times, limits=lim)
        common_times = [common_times, times]
        all_mlt_images = [all_mlt_images, mlt_images]
        ntime = n_elements(times)
        all_probes = [all_probes, probe+strarr(ntime)]
        all_time_ranges = [all_time_ranges, lim.time_range]
        all_hems = [all_hems,lim.hemisphere]
    endforeach
    index = sort(common_times)
    common_times = common_times[index]
    all_mlt_images = all_mlt_images[index,*,*]
    all_probes = all_probes[index]
    all_time_ranges = all_time_ranges[index,*]
    all_hems = all_hems[index]
    ssusi_unit = lim.unit
    
    ; Filter in time.
    plot_time_range = time_double(['2015-03-17/05:00','2015-03-18/00:00'])
    time_index = where_pro(common_times, '[]', plot_time_range)
    common_times = common_times[time_index]
    all_mlt_images = all_mlt_images[time_index,*,*]
    all_probes = all_probes[time_index]
    all_time_ranges = all_time_ranges[time_index,*]
    all_hems = all_hems[time_index]
    
    ; Filter in hemisphere
    time_index = where(all_hems eq 'SOUTH')
    common_times = common_times[time_index]
    all_mlt_images = all_mlt_images[time_index,*,*]
    all_probes = all_probes[time_index]
    all_time_ranges = all_time_ranges[time_index,*]
    all_hems = all_hems[time_index]
    
    ; Manual filter.
    time_index = [3,5,7,11,36,40,42,44]-1
    ;time_index = [3,5,7,9,11,13,17,21,38,40,42,44]-1
    common_times = common_times[time_index]
    all_mlt_images = all_mlt_images[time_index,*,*]
    all_probes = all_probes[time_index]
    all_time_ranges = all_time_ranges[time_index,*]
    all_hems = all_hems[time_index]

    npx = n_elements(all_mlt_images[0,0,*])
    if mlt_type eq 'pre_midn' then begin
        all_mlt_images = all_mlt_images[*,0:npx*0.5-1,0:npx*0.5-1]
    endif else if mlt_type eq 'post_midn' then begin
        all_mlt_images = all_mlt_images[*,npx*0.5-1:npx-1,0:npx*0.5-1]
    endif else if mlt_type eq 'night' then begin
        all_mlt_images = all_mlt_images[*,*,0:npx*0.5-1]
    endif else if mlt_type eq 'all_mlt' then begin
        ; do nothing. no need to crop.
    endif
    
    npan = n_elements(common_times)
    nxpan = 4d
    pansize = 2.5
    nypan = ceil(npan/nxpan)
    dst_ysize = 0.6
    

    margins = [1,1,7,1]
    if n_elements(plot_dir) eq 0 then plot_dir = event_info.plot_dir
    plot_file = join_path([plot_dir,'fig_dmsp_aurora_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    all_poss = panel_pos(plot_file, nxpan=nxpan,nypan=nypan+1,pansize=[1,dst_ysize]*pansize, ypans=[dst_ysize,fltarr(nypan)+0.5,0.5], $
        margins=margins, fig_size=fig_size, xpad=0, ypad=0)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    
    poss = all_poss[*,*,1:*]
    tpos = all_poss[*,0,0]
    tpos[2] = all_poss[2,-1,0]
    tpos[1] += ychsz*4
    tpos[0] += xchsz*8

 
    tmp = get_abs_chsz()
    abs_xchsz = tmp[0]
    abs_ychsz = tmp[1]
    uniform_ticklen = -abs_ychsz*0.15;*fig_size[0]
    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
    dst_var = omni_read_symh(time_range)
    options, dst_var, 'xticklen', xticklen
    options, dst_var, 'yticklen', yticklen
    ystep = 50
    exact_yrange = minmax(get_var_data(dst_var))
    yrange = minmax(make_bins(exact_yrange, ystep))
    for ii=1,10 do begin
        ytickv = make_bins(exact_yrange, ystep*ii, inner=1)
        yticks = n_elements(ytickv)-1
        if yticks le 2 then break
    endfor
    yminor = 5
    options, dst_var, 'yrange', yrange
    options, dst_var, 'ytickv', ytickv
    options, dst_var, 'yticks', yticks
    options, dst_var, 'yminor', yminor
    options, dst_var, 'constant', ytickv
    options, dst_var, 'ystyle', 1
    options, dst_var, 'labels', 'SymH'
    options, dst_var, 'ytitle', '(nT)'
    options, dst_var, 'num_lab_min', 12

    ;tplot, dst_var, trange=time_range, position=tpos, single_line_uttick=1
    xrange = time_range
    xstep = 3600*4
    xtickv = make_bins(xrange, xstep)
    xticks = n_elements(xtickv)-1
    xminor = 4
    xtickn = time_string(xtickv,tformat='hhmm')
    secofday = constant('secofday')
    index = where(xtickv mod secofday eq 0, count)
    if count ne 0 then begin
        xtickn[index] += '!C'+time_string(xtickv[index],tformat='YYYY-MM-DD')
    endif
    ytitle = '(nT)'
    
    dsts = get_var_data(dst_var, times=times)
    plot, times, dsts, $
        xstyle=1, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xrange=xrange, $
        ystyle=1, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, yrange=yrange, $
        position=tpos, noerase=1, $
        xtickname=xtickn, ytitle=ytitle
    
    tx = tpos[0]-xchsz*6
    ty = tpos[3]-ychsz*0.7
    msg = 'a) SymH'
    xyouts, tx,ty,msg, normal=1

    
    ct_ssusi = 70
    ssusi_zrange = [-1,1]*20

    color_top = 254
    nangle = 50
    dangle = !dpi/nangle
    angles = make_bins([-!dpi,!dpi], dangle)
    label_size = 0.8

    mlt_range = [-1,1]*6
    min_mlat = 50d
    mlat_range = [min_mlat,90]

    ; xlabels.
    xtickn_pos = (90.-min_mlat)/(90.-min_mlat+6)
    xrange = mlt_range
    xstep = 6d
    xtickv = smkarthm(xrange[0],xrange[1], 2, 'dx')
    xtick_minor = make_bins(xrange, 1, inner=1)
    xtickn = xtickv
    index = where(xtickv lt 0, count)
    if count ne 0 then xtickn[index] += 24
    xtickn = string(xtickn,format='(I02)')
    xtickn[0] = ' '
    xtickn[-1] = ' '
    
    ; ylabels.
    ytick_pos = (-1)*!dpi
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange, 10, inner=1)
    ytick_minor = make_bins(yrange, ystep, inner=1)
    ytickn = string(ytickv,format='(I0)')


    dst_pos = tpos
    fig_letter = 'b'
    foreach time, common_times, time_id do begin
        if all_hems[time_id] eq 'NORTH' then continue
        
        tpos = dst_pos
        xrange = time_range
        yrange = get_setting(dst_var, 'yrange')
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        plots, time+[0,0], yrange
        tmp = convert_coord(time, yrange[0], data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*0.9
        xyouts, tx,ty, normal=1, alignment=0.5, $
            string(time_id+1,format='(I0)'), charsize=label_size
        
        
        tmp = array_indices([nxpan,nypan], time_id, dimension=1)
        xid = tmp[0]
        yid = tmp[1]
        tpos = poss[*,xid,yid]
        

        mlt_image = reform(all_mlt_images[time_id,*,*])
        ssusi_zzs = bytscl(mlt_image, min=ssusi_zrange[0], max=ssusi_zrange[1], top=color_top)
        ;ssusi_time_range = reform(all_time_range[time_id,*])
        ssusi_time_range = time+[0,0]
        probe = all_probes[time_id]


        ; Draw data.
        zzs = ssusi_zzs
        ct = ct_ssusi
        sgtv, zzs, ct=ct, position=tpos

        ; Add labels, etc.
        ; Draw axes.
        if mlt_type eq 'pre_midn' then begin
            xrange = [-1,0]
            yrange = [-1,0]
            plot, xrange, yrange, nodata=1, noerase=1, $
                xstyle=5, ystyle=5, position=tpos
            plots, [0,-1], [0,0], color=sgcolor('silver')
            plots, [0,0], [0,-1], color=sgcolor('silver')
        endif else if mlt_type eq 'post_midn' then begin
            xrange = [0,1]
            yrange = [-1,0]
            plot, xrange, yrange, nodata=1, noerase=1, $
                xstyle=5, ystyle=5, position=tpos
            plots, [0,1], [0,0], color=sgcolor('silver')
            plots, [0,0], [0,-1], color=sgcolor('silver')
        endif else if mlt_type eq 'night' then begin
            xrange = [-1,1]
            yrange = [-1,0]
            plot, xrange, yrange, nodata=1, noerase=1, $
                xstyle=5, ystyle=5, position=tpos
            plots, [0,-1], [0,0], color=sgcolor('silver')
            plots, [0,1], [0,0], color=sgcolor('silver')
            plots, [0,0], [0,-1], color=sgcolor('silver')
        endif else if mlt_type eq 'all_mlt' then begin
            xrange = [-1,1]
            yrange = [-1,1]
            plot, xrange, yrange, nodata=1, noerase=1, $
                xstyle=5, ystyle=5, position=tpos
            plots, [0,-1], [0,0], color=sgcolor('silver')
            plots, [0,1], [0,0], color=sgcolor('silver')
            plots, [0,0], [0,-1], color=sgcolor('silver')
            plots, [0,0], [0,1], color=sgcolor('silver')
        endif



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
            ty = tmp[1]-ychsz*label_size*0.8
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
        

        ; Add panel label.
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]
        the_time_range = reform(all_time_ranges[time_id,*])
        msg = 'DMSP '+strupcase(probe)+' | '+str_cap(strlowcase(all_hems[time_id]))
        ;msg = time_string(time, tformat='YYYY-MM-DD/hh:mm')+' UT'
        xyouts, tx,ty+ychsz*1.2,normal=1, msg, charsize=label_size
        msg = fig_letter+'-'+string(time_id+1,format='(I0)')+') '+strjoin(time_string(the_time_range,tformat='hh:mm'),'-')+' UT'
        xyouts, tx,ty+ychsz*0.2,normal=1, msg;, charsize=label_size

        ; Colorbar.
        if xid eq nxpan-1 then begin
            ssusi_cbpos = tpos
            ssusi_cbpos[0] = ssusi_cbpos[2]+xchsz*1
            ssusi_cbpos[2] = ssusi_cbpos[0]+xchsz*0.7
            ssusi_cbpos[1] += ychsz*0.2
            ssusi_cbpos[3] -= ychsz*0.2
            ztitle = 'SSUSI '+strupcase(ssusi_id)+' ('+ssusi_unit+')'
            zrange = ssusi_zrange
            sgcolorbar, findgen(color_top), $
                ztitle=ztitle, zrange=zrange, ct=ct_ssusi, position=ssusi_cbpos, zticklen=zticklen, zminor=5
        endif
    endforeach


;---ASI.
    asi_times = time_double('2015-03-17/'+[$
        '05:45','06:51','07:25','09:10'])
        ;'06:00','06:51','07:25','09:10'])
    asi_time_range = minmax(asi_times)
    asi_poss = all_poss[*,*,-1]


    add_jver = 0
    mlt_type = 'night'
    asi_zlog = 0
    if asi_zlog eq 0 then asi_zrange = [-1,1]*4e3 else asi_zrange = [5e1,1e4]
    asi_ct = 70
    
;---SC setting.
    sc_list = list()
    sc_list.add, event_info.rbsp.rbspa
    sc_list.add, event_info.rbsp.rbspb
    nprobe = n_elements(sc_list)
;    model_setting = event_info.rbsp.rbspa['model_setting']
;    internal_model = event_info['internal_model']
;    external_model = event_info['external_model']
;    models = model_setting.models
;    model_index = where(models eq external_model)


;---Auto settings.
    common_times = asi_times
    ntime = n_elements(common_times)

    asi_setting = event_info['asi_setting']
    asi_sites = sort_uniq(asi_setting.sites)


    margins = [1,0.5,1,4]
    if mlt_type eq 'pre_midn' then begin
        pansize = [3,3]
    endif else if mlt_type eq 'post_midn' then begin
        pansize = [3,3]
    endif else if mlt_type eq 'night' then begin
        pansize = [6,3]
    endif else if mlt_type eq 'all_mlt' then begin
        pansize = [6,6]
    endif
    if add_jver then begin
        poss = panel_pos(nxpan=2, $
            pansize=pansize, margins=margins, fig_size=fig_size, xpad=2)
        asi_tpos = reform(poss[*,0])
        jver_tpos = reform(poss[*,1])
    endif else begin
        poss = panel_pos(pansize=pansize, margins=margins, fig_size=fig_size)
        asi_tpos = reform(poss[*,0])
    endelse
    top_color = 254



    ; Settings.
    if mlt_type eq 'pre_midn' then begin
        mlt_range = [-1d,0]*6
    endif else if mlt_type eq 'post_midn' then begin
        mlt_range = [0d,1]*6
    endif else if mlt_type eq 'night' then begin
        mlt_range = [-1,1]*6
    endif else if mlt_type eq 'all_mlt' then begin
        mlt_range = [-1,1]*12
    endif
    min_mlat = 50d
    mlat_range = [min_mlat,90]
    color_top = 254
    nangle = 50
    dangle = !dpi/nangle
    if mlt_type eq 'pre_midn' then begin
        angles = make_bins([-!dpi*0.5,0], dangle)
    endif else if mlt_type eq 'post_midn' then begin
        angles = make_bins([-!dpi,-!dpi*0.5], dangle)
    endif else if mlt_type eq 'night' then begin
        angles = make_bins([-!dpi,0], dangle)
    endif else if mlt_type eq 'all_mlt' then begin
        angles = make_bins([-!dpi,!dpi], dangle)
    endif
    ; xlabels.
    xtickn_pos = (90.-min_mlat)/(90.-min_mlat+6)
    xrange = mlt_range
    xstep = 6d
    xtickv = smkarthm(xrange[0],xrange[1], 2, 'dx')
    xtick_minor = make_bins(xrange, 1, inner=1)
    xtickn = xtickv
    index = where(xtickv lt 0, count)
    if count ne 0 then xtickn[index] += 24
    xtickn = string(xtickn,format='(I02)')
    xtickn[0] = ' '
    xtickn[-1] = ' '
    
    ; ylabels.
    ytick_pos = (-1)*!dpi
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange, 10, inner=1)
    ytick_minor = make_bins(yrange, ystep, inner=1)
    ytickn = string(ytickv,format='(I0)')


    ;---asi.
    asi_var = 'thg_asf_mlt_image_selected'
    if check_if_update(asi_var) then begin
        all_times = []
        all_images = []
        asi_setting = event_info['asi_setting']
        sites = asi_setting['sites']
        min_elevs = asi_setting['min_elevs']
        merge_method = asi_setting['merge_method']
        calibration_method = asi_setting['calibration_method']
        foreach asi_time, asi_times, time_id do begin
            time_range = asi_time+[-1d,1d]*3
            ;the_sites = themis_asf_get_sites_have_data(asi_time, sites=sites)
            the_sites = sites
            if asi_time ge time_double('2015-03-15/09:00') then begin
                index = where(sites eq 'nrsq' or sites eq 'gbay', complement=index2)
                the_sites = sites[index2]
            endif
            the_var = themis_asf_read_mlt_image(time_range, sites=the_sites, min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method)
            get_data, the_var, times, mlt_images, limits=lim
            all_images = [all_images, mlt_images]
            all_times = [all_times, times]
        endforeach
        store_data, asi_var, all_times, all_images, limits=lim
    endif
    

    get_data, asi_var, times, mlt_images
    npx = n_elements(mlt_images[0,0,*])
    if mlt_type eq 'pre_midn' then begin
        mlt_images = mlt_images[*,0:npx*0.5-1,0:npx*0.5-1]
    endif else if mlt_type eq 'post_midn' then begin
        mlt_images = mlt_images[*,npx*0.5-1:npx-1,0:npx*0.5-1]
    endif else if mlt_type eq 'night' then begin
        mlt_images = mlt_images[*,*,0:npx*0.5-1]
    endif else if mlt_type eq 'all_mlt' then begin
        ; do nothing. no need to crop.
    endif
    
    image_size = size(reform(mlt_images[0,*,*]), dimensions=1)
    asi_mlt_images = fltarr([ntime,image_size])
    for ii=0,image_size[0]-1 do for jj=0,image_size[1]-1 do asi_mlt_images[*,ii,jj] = interpol(mlt_images[*,ii,jj],times, common_times)
    index = where(finite(asi_mlt_images,nan=1), count)
    if count ne 0 then asi_mlt_images[index] = 0

    if asi_zlog eq 1 then begin
        asi_log_zrange = alog10(asi_zrange)
        asi_log_ztickv = make_bins(asi_log_zrange,1, inner=1)
        asi_ztickv = 10^asi_log_ztickv
        asi_zticks = n_elements(asi_ztickv)-1
        asi_zminor = 9
        asi_ztickn = '10!U'+string(asi_log_ztickv,format='(I0)')

        asi_zzs = bytscl(alog10(asi_mlt_images), min=asi_log_zrange[0],max=asi_log_zrange[1], top=color_top)
    endif else begin
        asi_ztickv = sort_uniq([asi_zrange,0])
        asi_zticks = n_elements(asi_ztickv)-1
        asi_zminor = 10
        asi_ztickn = string(asi_ztickv,format='(I0)')

        asi_zzs = bytscl((asi_mlt_images), min=asi_zrange[0],max=asi_zrange[1], top=color_top)
    endelse
    asi_ztitle = 'N-hem ASI (#)'
    asi_zticks = n_elements(asi_ztickv)-1
    asi_linestyle = 1


    ; Loop through the times.
    foreach time, asi_times, time_id do begin
        asi_tpos = asi_poss[*,time_id]

    ;---ASI.
        if time_id eq nxpan-1 then begin
            asi_cbpos = asi_tpos
            asi_cbpos[0] = asi_cbpos[2]+xchsz*1
            asi_cbpos[2] = asi_cbpos[0]+xchsz*0.7
            asi_cbpos[1] += ychsz*0.2
            asi_cbpos[3] -= ychsz*0.2
            asi_zticklen = uniform_ticklen/xchsz*1/fig_size[0]
            sgcolorbar, findgen(top_color), horizontal=0, $
                ztitle=asi_ztitle, zrange=asi_zrange, ct=asi_ct, position=asi_cbpos, $
                ztickv=asi_ztickv, zticks=asi_zticks, ztickname=asi_ztickn, zminor=asi_zminor, log=asi_zlog, zticklen=asi_zticklen
        endif

        
        tpos = asi_tpos
        zzs = reform(asi_zzs[time_id,*,*])
        sgtv, zzs, ct=asi_ct, position=tpos


    ;---Add axis, labels, etc.
        info_list = list()
        str_id = string(time_id+1,format='(I0)')
        info_list.add, dictionary($
            'msg', ['c-'+str_id+') '+time_string(time,tformat='hh:mm:ss')+' UT'], $
            'hemisphere', 'north', $
            'position', asi_tpos, $
            'ct', asi_ct )

        foreach the_info, info_list do begin
            tpos = the_info.position

            ; Add labels, etc.
            ; Draw axes.
            if mlt_type eq 'pre_midn' then begin
                xrange = [-1,0]
                yrange = [-1,0]
                plot, xrange, yrange, nodata=1, noerase=1, $
                    xstyle=5, ystyle=5, position=tpos
                plots, [0,-1], [0,0], color=sgcolor('silver')
                plots, [0,0], [0,-1], color=sgcolor('silver')
            endif else if mlt_type eq 'post_midn' then begin
                xrange = [0,1]
                yrange = [-1,0]
                plot, xrange, yrange, nodata=1, noerase=1, $
                    xstyle=5, ystyle=5, position=tpos
                plots, [0,1], [0,0], color=sgcolor('silver')
                plots, [0,0], [0,-1], color=sgcolor('silver')
            endif else if mlt_type eq 'night' then begin
                xrange = [-1,1]
                yrange = [-1,0]
                plot, xrange, yrange, nodata=1, noerase=1, $
                    xstyle=5, ystyle=5, position=tpos
                plots, [0,-1], [0,0], color=sgcolor('silver')
                plots, [0,1], [0,0], color=sgcolor('silver')
                plots, [0,0], [0,-1], color=sgcolor('silver')
            endif else if mlt_type eq 'all_mlt' then begin
                xrange = [-1,1]
                yrange = [-1,1]
                plot, xrange, yrange, nodata=1, noerase=1, $
                    xstyle=5, ystyle=5, position=tpos
                plots, [0,-1], [0,0], color=sgcolor('silver')
                plots, [0,1], [0,0], color=sgcolor('silver')
                plots, [0,0], [0,-1], color=sgcolor('silver')
                plots, [0,0], [0,1], color=sgcolor('silver')
            endif

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
                ty = tmp[1]-ychsz*label_size*0.9
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

            
            
            ; Add panel label.
            tx = tpos[0]+xchsz*0.5
            ty = tpos[1]+ychsz*0.2
            msgs = the_info.msg
            xyouts, tx,ty,normal=1, msgs[0];, charsize=label_size
            ;xyouts, tx,ty+ychsz*1,normal=1, msgs[0]
            msg = 'ASI | North'
            xyouts, tx,ty+ychsz*1,normal=1, msg, charsize=label_size
        endforeach
    endforeach

    if keyword_set(test) then stop
    sgclose

    return, plot_file

end


print, fig_dmsp_aurora_v02(test=0, event_info=event_info)
end