;+
; Plot DMSP aurora over the storm.
;-

function fig_dmsp_aurora_v01, plot_file, event_info=event_info, test=test


    version = 'v01'
    id = '2015_0317'
    if n_elements(event_info) eq 0 then event_info = low_lshell_outflow_load_data(id)
    time_range = event_info.time_range
    dmsp_probes = 'f'+['16','17','18','19']


;---Load data.
    ssusi_id = 'energy'
    ssusi_wavelength = strupcase(ssusi_id)
    data_time_range = time_range+[-1800,0]
    ssusi_info = dictionary()
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
    
    npan = n_elements(common_times)
    nxpan = 4d
    pansize = 2.5
    nypan = ceil(npan/nxpan)
    
    

    margins = [1,1,7,1]
    if n_elements(plot_dir) eq 0 then plot_dir = event_info.plot_dir
    plot_file = join_path([plot_dir,'fig_dmsp_aurora_v01.pdf'])
    if keyword_set(test) then plot_file = 0
    all_poss = panel_pos(plot_file, nxpan=nxpan,nypan=nypan+1,pansize=[1,1]*pansize, $
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
    
    ct_ssusi = 70
    ssusi_zrange = [-1,1]*20

    color_top = 254
    nangle = 50
    dangle = !dpi/nangle
    angles = make_bins([-!dpi,!dpi], dangle)
    label_size = 0.8

    mlt_range = [-1,1]*12
    min_mlat = 50d
    mlat_range = [min_mlat,90]

    ; xlabels.
    xtickn_pos = (90.-min_mlat)/(90.-min_mlat+6)
    xrange = mlt_range
    xstep = 6d
    xtickv = smkarthm(xrange[0],xrange[1], 6, 'dx')
    xtick_minor = make_bins(xrange, 1, inner=1)
    xtickn = xtickv
    index = where(xtickv lt 0, count)
    if count ne 0 then xtickn[index] += 24
    xtickn = string(xtickn,format='(I02)')
    
    ; ylabels.
    ytick_pos = (-1+1d/12)*!dpi
    ytick_pos = (0.25)*!dpi
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
        plot, [-1,1], [-1,1], nodata=1, noerase=1, $
            xstyle=5, ystyle=5, position=tpos
;        plots, [0,-1], [0,0], color=sgcolor('silver')
;        plots, [0,0], [0,-1], color=sgcolor('silver')
;        plots, [0,1], [0,0], color=sgcolor('silver')
;        plots, [0,0], [0,1], color=sgcolor('silver')


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
        

        ; Add panel label.
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*2
        the_time_range = all_time_ranges[time_id,*]
        msg = time_string(the_time_range[0],tformat='YYYY-MM-DD/hh:mm')+$
            '-'+time_string(the_time_range[1],tformat='hh:mm')+' UT'
        msg = time_string(time, tformat='YYYY-MM-DD/hh:mm')+' UT'
        xyouts, tx,ty,normal=1, msg, charsize=label_size
        msg = fig_letter+'-'+string(time_id+1,format='(I0)')+') DMSP '+strupcase(probe)+' | '+all_hems[time_id]
        xyouts, tx,ty+ychsz*1,normal=1, msg, charsize=label_size

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

    if keyword_set(test) then stop
    sgclose

    return, plot_file

end


print, fig_dmsp_aurora_v01(test=0, event_info=event_info)
end