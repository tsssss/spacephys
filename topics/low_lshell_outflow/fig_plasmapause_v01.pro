;+
; Plasmapause location
;-

function scale_pp_time_to_zz, pp_time, time_range

    return, bytscl(pp_time, min=time_range[0], max=time_range[1], top=254)
    
end

function fig_plasmapause_v01_config_panel, config_pos, $
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
    foreach rr, make_bins([2,6],1) do begin
        oplot, rr*txs, rr*tys, linestyle=1, color=sgcolor('silver')
        tt = 145*constant('rad')
        tx = rr*cos(tt)
        ty = rr*sin(tt)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*0.3
        msg = 'L = '+string(rr,format='(I0)')
        xyouts, tx,ty,msg, normal=1, alignment=0.5, charsize=label_size
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
    

;---RBSP.
    rbsp_probes = event_info.rbsp_probes
    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'
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
    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'
        r_var = prefix+'r_sm'
        foreach pp_info, pp_infos do begin
            if pp_info.mission ne 'rbsp' then continue
            if pp_info.probe ne probe then continue
            pp_time = pp_info.time
            zz = scale_pp_time_to_zz(pp_time, time_range)
            color = sgcolor(zz, ct=ct)
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

    tx = -3.5
    ty = 5.5
    msg = 'RBSP-A & B'
    xyouts, tx,ty, data=1, msg;, charsize=label_size


;;---THEMIS.
;    themis_probes = event_info.themis_probes
;    foreach probe, themis_probes do begin
;        prefix = 'th'+probe+'_'
;        mlt = get_var_data(prefix+'mlt', times=times)
;        lshell = get_var_data(prefix+'lshell', at=times)
;        tt = (mlt*15+180)*constant('rad')
;        xxs = lshell*cos(tt)
;        yys = lshell*sin(tt)
;        oplot, xxs, yys, linestyle=0, color=sgcolor('silver')
;    endforeach
;    
;    foreach probe, themis_probes do begin
;        prefix = 'th'+probe+'_'
;        foreach pp_info, pp_infos do begin
;            if pp_info.mission ne 'themis' then continue
;            if pp_info.probe ne probe then continue
;            pp_time = pp_info.time
;            zz = scale_pp_time_to_zz(pp_time, time_range)
;            color = sgcolor(zz, ct=ct)
;            mlt = get_var_data(prefix+'mlt', at=pp_time)
;            lshell = get_var_data(prefix+'lshell', at=pp_time)
;            tt = (mlt*15+180)*constant('rad')
;            xxs = lshell*cos(tt)
;            yys = lshell*sin(tt)
;            plots, xxs, yys, psym=pp_psym, color=color, symsize=pp_symsize
;        endforeach
;    endforeach
;    
;    tx = 2
;    ty = 5
;    msg = 'THEMIS'
;    xyouts, tx,ty, data=1, msg


;--Add label.
    tx = tpos[0]-xchsz*4
    ty = tpos[3]-ychsz*0.8
    msg = 'a)'
    xyouts, tx,ty,normal=1, msg
    
    foreach rr, make_bins([2,6],1) do begin
        tt = 145*constant('rad')
        tx = rr*cos(tt)
        ty = rr*sin(tt)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]
        ty = tmp[1]-ychsz*0.3
        msg = 'L = '+string(rr,format='(I0)')
        xyouts, tx,ty,msg, normal=1, alignment=0.5, charsize=label_size
    endforeach
    
end



function fig_plasmapause_v01_right_panels, right_pos, event_info=event_info

    time_range = event_info.time_range
    dst_var = omni_read_symh(time_range, get_name=1)
    rbsp_probes = event_info.rbsp_probes

    plot_vars = [dst_var, 'rbsp'+rbsp_probes+'_density']
    nvar = n_elements(plot_vars)
    fig_labels = letters(minmax(1+findgen(nvar+1)))+')'
    
    
    ypans = fltarr(nvar)+1
    ypans[0] = 1.2
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

    foreach pp_info, pp_infos do begin
        probe = pp_info.probe
        prefix = pp_info.mission_probe+'_'
        pp_time = pp_info.time
        zz = scale_pp_time_to_zz(pp_time, time_range)
        color = sgcolor(zz, ct=ct)
        dens_var = prefix+'density'
        
        pid = where(plot_vars eq dens_var, count)
        if count eq 0 then continue
        tpos = right_poss[*,pid]
        yr = get_setting(dens_var, 'yrange')
        xr = time_range
        ylog = get_setting(dens_var, 'ylog')
        plot, xr, yr, $
            xrange=xr, yrange=yr, xstyle=5, ystyle=5, ylog=ylog, $
            position=tpos, nodata=1, noerase=1
        tmp = get_var_data(dens_var, at=pp_time)
        if pp_psym eq 8 then begin
            usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1
            plots, pp_time, tmp, psym=pp_psym, color=color, symsize=pp_symsize
            usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=0
            plots, pp_time, tmp, psym=pp_psym, symsize=pp_symsize
        endif else begin
            plots, pp_time, tmp, psym=pp_psym, color=color, symsize=pp_symsize
            plots, pp_time, tmp, psym=pp_psym, symsize=pp_symsize
        endelse
    endforeach
    
    pid = where(plot_vars eq dst_var, count)
    if count ne 0 then begin
        tpos = right_poss[*,pid]
        yticklen = get_setting(dst_var, 'yticklen')
        yr = get_setting(dst_var, 'yrange')
        xr = time_range
        plot, xr, yr, $
            xrange=xr, yrange=yr, xstyle=5, ystyle=5, $
            position=tpos, nodata=1, noerase=1
        
        yrange = [2,5.5]
        ytickv = [2,3,4,5]
        yticks = n_elements(ytickv)-1
        yminor = 4
        ytitle = 'L-shell'
        axis, yaxis=1, save=1, ystyle=1, $
            yticklen=yticklen, yrange=yrange, $
            ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle
        
        foreach pp_info, pp_infos do begin
            if pp_info.mission eq 'themis' then continue
            pp_time = pp_info.time
            zz = scale_pp_time_to_zz(pp_time, time_range)
            color = sgcolor(zz, ct=ct)
            prefix = pp_info.mission_probe+'_'
            tmp = get_var_data(prefix+'lshell', at=pp_time)
            if pp_psym eq 8 then begin
                usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=1
                plots, pp_time, tmp, psym=pp_psym, color=color, symsize=pp_symsize
                usersym, circ_xs[0:*:2], circ_ys[0:*:2], fill=0
                plots, pp_time, tmp, psym=pp_psym, symsize=pp_symsize
            endif else begin
                plots, pp_time, tmp, psym=pp_psym, color=color, symsize=pp_symsize
                plots, pp_time, tmp, psym=pp_psym, symsize=pp_symsize
            endelse
        endforeach
    endif


end




function fig_plasmapause_v01, plot_file, event_info=event_info, test=test

    version = 'v01'
    id = '2015_0317'
    event_info = low_lshell_outflow_load_data(id)
        
    time_range = event_info.time_range
    rbsp_probes = event_info.rbsp_probes
    themis_probes = event_info.themis_probes
    
    plot_dir = event_info['plot_dir']
    base = 'fig_plasmapause_'+time_string(time_range[0],tformat='YYYY_MMDD')+'_'+version+'.pdf'
    if n_elements(plot_file) eq 0 then plot_file = join_path([plot_dir,base])
    if keyword_set(test) then plot_file = 0

    pp_infos = list()

;---THEMIS.
    dens_range = 1000+[-1,1]*10
    pad_time = 3600

    foreach probe, themis_probes do begin
        prefix = 'th'+probe+'_'
        
        pp_trs = list()
        dens_var = themis_read_density_efi(time_range, probe=probe, suffix='')
        options, dens_var, 'constant', mean(dens_range)
        options, dens_var, 'yrange', [0.3,5e3]
        options, dens_var, 'labels', strupcase('TH-'+probe)

        pp_var = prefix+'pp_times'
        get_data, dens_var, times, dens
        ;index = where_pro(dens, '[]', dens_range, count=count)
        tmp = dens-mean(dens_range)
        index = where(tmp[1:-1]*tmp[0:-2] le 0, count)
        if count ne 0 then begin
            pp_list = list()
            pp_times = 0.5*(times[index]+times[index+1])
            for ii=0,count-1 do begin
                if n_elements(pp_list) eq 0 then begin
                    pp_list.add, pp_times[ii]
                    continue
                endif
                
                if abs(pp_list[-1]-pp_times[ii]) le pad_time then begin
                    pp_list.add, pp_times[ii]
                    continue
                endif
                
                if n_elements(pp_list) ne 0 then pp_trs.add, minmax(pp_list.toarray())
                pp_list = list(pp_times[ii])
            endfor
            if n_elements(pp_list) ne 0 then pp_trs.add, minmax(pp_list.toarray())
            pp_trs = pp_trs.toarray()
            pp_times = (pp_trs[*,0]+pp_trs[*,1])*0.5
            store_data, pp_var, pp_times, pp_trs, limits={dens_range:dens_range}
        endif

        r_var = themis_read_orbit(time_range, probe=probe, coord='sm')
        mlt_var = themis_read_mlt(time_range, probe=probe)
        lshell_var = themis_read_lshell(time_range, probe=probe)
        
        pp_times = get_var_time(pp_var)
        pp_mlts = get_var_data(mlt_var, at=pp_times)
        pp_diss = snorm(get_var_data(r_var, at=pp_times))
        
        npp_time = n_elements(pp_times)
        for ii=0,npp_time-1 do begin
            if pp_diss[ii] lt 2 then continue
            if pp_diss[ii] ge 6 then continue
            pp_infos.add, dictionary($
                'mission_probe', 'th'+probe, $
                'mission', 'themis', $
                'probe', probe, $
                'time', pp_times[ii], $
                'mlt', pp_mlts[ii], $
                'dis', pp_diss[ii] )
        endfor
    endforeach


;---RBSP.
    dens_range = 100+[-1,1]*10
    pad_time = 3600
    foreach probe, rbsp_probes do begin
        prefix = 'rbsp'+probe+'_'
        
        pp_trs = list()
        dens_var = rbsp_read_density(time_range, probe=probe, id='emfisis', suffix='')
        options, dens_var, 'constant', mean(dens_range)
        options, dens_var, 'yrange', [0.3,5e3]
        options, dens_var, 'labels', strupcase('RBSP-'+probe)
        
        pp_var = prefix+'pp_times'
        get_data, dens_var, times, dens
        ;index = where_pro(dens, '[]', dens_range, count=count)
        tmp = dens-mean(dens_range)
        index = where(tmp[1:-1]*tmp[0:-2] le 0, count)
        if count ne 0 then begin
            pp_list = list()
            pp_times = 0.5*(times[index]+times[index+1])
            for ii=0,count-1 do begin
                if n_elements(pp_list) eq 0 then begin
                    pp_list.add, pp_times[ii]
                    continue
                endif
                
                if abs(pp_list[-1]-pp_times[ii]) le pad_time then begin
                    pp_list.add, pp_times[ii]
                    continue
                endif
                
                if n_elements(pp_list) ne 0 then pp_trs.add, minmax(pp_list.toarray())
                pp_list = list(pp_times[ii])
            endfor
            if n_elements(pp_list) ne 0 then pp_trs.add, minmax(pp_list.toarray())
            pp_trs = pp_trs.toarray()
            pp_times = (pp_trs[*,0]+pp_trs[*,1])*0.5
            store_data, pp_var, pp_times, pp_trs, limits={dens_range:dens_range}
        endif
        
        
        r_var = rbsp_read_orbit(time_range, probe=probe, coord='sm')
        mlt_var = rbsp_read_mlt(time_range, probe=probe)
        lshell_var = rbsp_read_lshell(time_range, probe=probe)
        
        pp_times = get_var_time(pp_var)
        pp_mlts = get_var_data(mlt_var, at=pp_times)
        pp_diss = snorm(get_var_data(r_var, at=pp_times))
        
        npp_time = n_elements(pp_times)
        for ii=0,npp_time-1 do begin
            pp_infos.add, dictionary($
                'mission_probe', 'rbsp'+probe, $
                'mission', 'rbsp', $
                'probe', probe, $
                'time', pp_times[ii], $
                'mlt', pp_mlts[ii], $
                'dis', pp_diss[ii] )
        endfor
    endforeach
    store_data, 'pp_infos', 0, pp_infos
    
    
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
    


;---Plot settings.
    tmp = get_abs_chsz()
    abs_xchsz = tmp[0]
    abs_ychsz = tmp[1]

    pos_xrange = [2.5,-6.5]
    pos_yrange = [6,-3]
    margins = [6,4,8,1]

    pos_aspect_ratio = abs(total(pos_xrange*[-1,1])/total(pos_yrange*[-1,1]))
    pos_xpan_size = 3
    pos_ypan_size = pos_xpan_size/pos_aspect_ratio

    right_xpan_size = pos_xpan_size*1.5


    all_poss = panel_pos(plot_file, xpans=[pos_xpan_size,right_xpan_size], xpad=10, margins=margins, panid=[0,0], pansize=[pos_xpan_size,pos_ypan_size], fig_size=fig_size)

    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    

;---Pos.
    config_pos = all_poss[*,0]
    tpos = fig_plasmapause_v01_config_panel(config_pos, $
        pos_xrange=pos_xrange, pos_yrange=pos_yrange, event_info=event_info)
    right_pos = all_poss[*,1]
    poss = fig_plasmapause_v01_right_panels(right_pos, event_info=event_info)
    
    
    
    if keyword_set(test) then stop
    sgclose

    return, plot_file

end


print, fig_plasmapause_v01(test=1)
end