
pro azim_df_plot_time_lag, event_time_range, event_id=event_id, project=project

    test = 0

;---Check inputs.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(event_time_range) ne 2 then begin
        errmsg = handle_error('No input event_time_range ...')
        return
    endif
    candidates = project.candidate
    foreach candidate, candidates do begin
        full_time_range = candidate.time_range
        if max(event_time_range) lt min(full_time_range) then continue
        if min(event_time_range) gt max(full_time_range) then continue
        event_id = candidate.id
        break
    endforeach
    if n_elements(event_id) eq 0 then begin
        errmsg = handle_error('No input time_range or event_id ...')
        return
    endif
    candidate = candidates[event_id]
    full_time_range = candidate.time_range

;---Load basic data and event_info.
    azim_df_load_basic_data, full_time_range, event_id=event_id, project=project, reset=reset
    event_infos = azim_df_load_event_info(project=project)
    if ~event_infos.haskey(event_id) then begin
        errmsg = handle_error('Cannot find event_id: '+event_id+' ...')
        return
    endif
    event_info = get_var_data('event_info')
    save_vars = tnames('*')
    the_info = event_infos[event_id]
    if n_elements(event_time_range) ne 2 then event_time_range = the_info.time_range
    sorted_probes = event_info.sorted_probes
    cc_time_range = the_info.cc_time_range


;---Plot settings.
    margins = [12,3,2,2]
    fig_xsize = 6
    aspect_ratio = 0.3
    ytitle = '(deg)'
    xticklen = -0.06
    yticklen = xticklen*aspect_ratio
    half_ychsz = 0.35
    full_ychsz = 0.7
    line_skip = 0.35
    label_size = 0.7
    title_size = 1.0
    xpad = 10
    ypad = 0.2
    nvar = n_elements(sorted_probes)
    probe_colors = smkarthm(50,250,nvar, 'n')
    ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=ct)
    sgopen, 0, xsize=fig_xsize, ysize=fig_xsize
    tpos = sgcalcpos(1,2, xchsz=xchsz, ychsz=ychsz, ypad=ypad, xpad=xpad, margins=margins)
    sgclose, /wdelete
    pan_xsize = fig_xsize*(tpos[2,0]-tpos[0,0])
    pan_ysize = pan_xsize*aspect_ratio
    fig_ysize = pan_ysize*nvar+(ypad*(nvar-1)+(margins[1]+margins[3]))*ychsz*fig_xsize
    if keyword_set(test) then begin
        file = test
        magnify = 2
    endif else begin
        file = join_path([project.plot_dir,'time_lag','fig_theta_'+event_id+'.pdf'])
        magnify = 1
    endelse
    thick = (size(file,/type) eq 7)? 4:2
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify
    poss = sgcalcpos(nvar,2, margins=margins, xchsz=xchsz,ychsz=ychsz, ypad=ypad, xpad=xpad)

;---Plot B tilt angle.
    str_theta = tex2str('theta')
    str_delta = tex2str('Delta')
    plot_title = 'Detrended tilt angle '+str_theta+' in absolute time'

    xrange = event_time_range
    xstep = 30*60
    xtickv = smkarthm(xrange[0],xrange[1],xstep, 'dx')
    xticks = n_elements(xtickv)-1
    xminor = 6
    while xticks ge 4 do begin
        xstep *= 2
        xtickv = smkarthm(xrange[0],xrange[1],xstep, 'dx')
        xticks = n_elements(xtickv)-1
        xminor *= 2
    endwhile
    xtickn = strarr(xticks+1)
    for ii=0, xticks do xtickn[ii] = time_string(xtickv[ii],tformat='hh:mm')
    xtickn[0] = time_string(xtickv[0],tformat='YYYY-MM-DD') & for ii=0, strlen(xtickn[0])-1 do xtickn[0] += ' '

;---Loop through each sorted_probe.
    var_suffix = '_theta'
    probe_yranges = list()
    foreach probe, sorted_probes, ii do begin
        get_data, probe+var_suffix, times, dbtilt
        index = where_pro(times, event_time_range)
        times = times[index]
        dbtilt = dbtilt[index]

        tpos = reform(poss[*,0,ii])
        xtickformat = (ii eq nvar-1)? '': '(A1)'

        yrange = sg_autolim(minmax(dbtilt))
        yrange = [-1,1]*max(abs(yrange))
        probe_yranges.add, yrange
        yrange = round(yrange/2)*2
        yticks = 2
        ytickv = smkarthm(yrange[0],yrange[1],yticks+1,'n')
        ytickn = strarr(yticks+1)
        yminor = 5
        for jj=0, yticks do ytickn[jj] = sgnum2str(ytickv[jj])
        plot, times, dbtilt, $
            xstyle=1, xrange=xrange, xtickformat=xtickformat, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, $
            ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytickname=ytickn, $
            position=tpos, /noerase, /nodata
        oplot, times, dbtilt, color=sgcolor('black')
        tx = event_info[probe].arrival_time
        plots, tx+[0,0], yrange, linestyle=1

        if ii eq 0 then begin
            tx = (tpos[0]+tpos[2])*0.5
            ty = tpos[3]+ychsz*line_skip
            xyouts, tx,ty,/normal,alignment=0.5, plot_title, charsize=title_size
            ; add the time range for cc.
            plot, event_time_range, [0,1], /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
            ty = 0.1
            txs = cc_time_range
            tys = ty+[0,0]
            plots, txs,tys
            tmp = convert_coord(mean(txs),ty, /data, /to_normal)
            xyouts, tmp[0], tmp[1]+0.2*ychsz*label_size, alignment=0.5, /normal, 'data for c.c', charsize=label_size

            for jj=0,1 do begin
                tmp = convert_coord(txs[jj],tys[0], /data, /to_normal)
                tx = tmp[0]
                ty = tmp[1]
                plots, tx+[0,0], ty+[-1,1]*ychsz*0.15, /normal
            endfor
        endif
        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*half_ychsz
        xyouts, tx,ty,/normal, 'a-'+sgnum2str(ii+1)+'.'

        ; Add the error in c.c.
        if ii ne 0 then begin
            cc_info = event_info.cc_info
            plot, event_time_range, [0,1], /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
            ty = 0.1
            cc_error = cc_info[probe].time_lag_error
            txs = event_info[probe].arrival_time+[-1,1]*cc_error
            tys = ty+[0,0]
            plots, txs,tys, thick=thick
            tmp = convert_coord(max(txs),ty, /data, /to_normal)
            if ii eq 1 then begin
                xyouts, tmp[0]+xchsz*0.25, tmp[1]-0.2*ychsz*label_size, /normal, 'uncertainty of c.c', charsize=label_size
            endif

            for jj=0,1 do begin
                tmp = convert_coord(txs[jj],tys[0], /data, /to_normal)
                tx = tmp[0]
                ty = tmp[1]
                plots, tx+[0,0], ty+[-1,1]*ychsz*0.05, /normal
            endfor
        endif

        ; s/c name, MLT and R.
        tx = xchsz*5
        ty = mean(tpos[[1,3]])+ychsz*half_ychsz
        alignment = 1
        probe_info = resolve_probe(probe)
        short_name = probe_info.short_name
        xyouts, tx,ty,/normal, alignment=alignment, strupcase(short_name), color=probe_colors[ii]
        xyouts, tx,ty-ychsz*label_size*1,/normal, alignment=alignment, sgnum2str(event_info[probe].arrival_mlt,ndec=1)+' MLT', charsize=label_size
        arrival_rsm = event_info[probe].arrival_rsm
        xyouts, tx,ty-ychsz*label_size*2,/normal, alignment=alignment, sgnum2str(snorm(arrival_rsm[0:1]),ndec=1)+' Re', charsize=label_size
        if ii eq 0 then begin
            tx = (tpos[0]+tpos[2])*0.5
            ty = tpos[3]+ychsz*line_skip
            xyouts, tx,ty,/normal,alignment=0.5, plot_title, charsize=title_size
        endif
    endforeach


;---Time lag removed.
    plot_title = 'Time lag removed'
    foreach probe, sorted_probes, ii do begin
        get_data, probe+var_suffix, times, dbtilt
        tpos = reform(poss[*,1,ii])
        xtickformat = (ii eq nvar-1)? '': '(A1)'

        yrange = probe_yranges[ii]
        yrange = round(yrange/2)*2
        yticks = 2
        ytickv = smkarthm(yrange[0],yrange[1],yticks+1,'n')
        ytickn = strarr(yticks+1)
        yminor = 5
        for jj=0, yticks do ytickn[jj] = sgnum2str(ytickv[jj])
        plot, times, dbtilt, $
            xstyle=1, xrange=xrange, xtickformat=xtickformat, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, $
            ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytickname=ytickn, $
            position=tpos, /noerase, /nodata
        cc_info = event_info.cc_info
        time_lag = event_info[probe].time_lag
        max_corr = cc_info[probe].max_corr
        times = times-time_lag
        index = where_pro(times, event_time_range)
        times = times[index]
        dbtilt = dbtilt[index]
        oplot, times, dbtilt, color=sgcolor('gray')
        arrival_time = event_info[probe].arrival_time-time_lag
        tx = arrival_time
        plots, tx+[0,0], yrange, linestyle=1, color=sgcolor('gray')
        ; Add label.
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, str_delta+'t = '+sgnum2str(time_lag)+' sec', charsize=label_size
        if max_corr ne 0 then begin
            ;tx = tpos[0]+xchsz*8
            ty = tpos[3]-ychsz*full_ychsz*2
            xyouts, tx,ty,/normal, 'c.corr = '+sgnum2str(max_corr,ndec=2), charsize=label_size
        endif

        ; Add the error in c.c.
        if ii ne 0 then begin
            plot, event_time_range, [0,1], /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
            ty = 0.1
            cc_error = cc_info[probe].time_lag_error
            txs = arrival_time+[-1,1]*cc_error
            tys = ty+[0,0]
            plots, txs,tys, thick=thick
            tmp = convert_coord(max(txs),ty, /data, /to_normal)
            xyouts, tmp[0]+xchsz*0.25, tmp[1]-0.2*ychsz*label_size, /normal, tex2str('pm')+sgnum2str(round(cc_error))+' sec', charsize=label_size
            
            for jj=0,1 do begin
                tmp = convert_coord(txs[jj],tys[0], /data, /to_normal)
                tx = tmp[0]
                ty = tmp[1]
                plots, tx+[0,0], ty+[-1,1]*ychsz*0.05, /normal
            endfor
        endif

        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*half_ychsz
        xyouts, tx,ty,/normal, 'b-'+sgnum2str(ii+1)+'.'

        if ii eq 0 then begin
            tx = (tpos[0]+tpos[2])*0.5
            ty = tpos[3]+ychsz*line_skip
            xyouts, tx,ty,/normal,alignment=0.5, plot_title, charsize=title_size
        endif
    endforeach


    if keyword_set(test) then stop
    sgclose
end

time_range = time_double(['2014-08-28/09:30','2014-08-28/11:30'])
;time_range = time_double(['2017-03-28/02:00','2017-03-28/05:00'])
azim_df_plot_time_lag, time_range, project=project
end
