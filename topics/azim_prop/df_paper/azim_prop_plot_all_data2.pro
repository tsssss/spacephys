;+
; Plot all the data for all events, each panel has all spacecraft data shifted by the proper time lags.
;-

pro azim_prop_plot_all_data2, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    test = 0

    mlt_colors = list()
    mlt_values = list()
    ncolor = 175
    top_color = 255
    max_mlt = 7
    dawn_ct = 49
    dusk_ct = 50
    foreach key, ['dawn','dusk'] do begin
        ct = (key eq 'dawn')? dawn_ct: dusk_ct
        mlt_range = (key eq 'dawn')? [0.,max_mlt]: [0.,-max_mlt]
        mlt_values.add, smkarthm(mlt_range[0],mlt_range[1],ncolor,'n'), /extract
        colors = lonarr(ncolor)
        for ii=0, ncolor-1 do colors[ii] = sgcolor(top_color-ii, ct=ct)
        mlt_colors.add, colors, /extract
    endforeach

    mlt_colors = mlt_colors.toarray()
    mlt_values = mlt_values.toarray()
    index = sort(mlt_values)
    mlt_values = mlt_values[index]
    mlt_colors = mlt_colors[index]
    dusk_ct = 50    ; green.


;---Plot settings.
    re = project.constant.re
    deg = project.constant.deg
    rad = project.constant.rad
    str_pm = '!9'+string(177b)+'!X'
    str_theta = '!9'+string(113b)+'!X'
    str_delta = '!9'+string(68b)+'!X'
    str_omega = '!9'+string(119b)+'!X'
    r_kms = 6.6
    coef2kms = r_kms*re*rad/60
    data_type = 'db_tilt'   ; the data to be plotted.
    var_unit = 'deg'
    data_step = 2
    lmargin = 2
    label_size = 0.7

    ; Sort event by date.
    sorted_events = (events.keys()).toarray()
    index = sort(time_double(sorted_events,tformat='YYYY_MM_DD'))
    sorted_events = sorted_events[index]

    ; Figure out the positions of the panels.
    nevent = n_elements(sorted_events)
    panel_aspect_ratio = 0.2    ; ysize/xsize.
    xticklen = -0.06
    yticklen = xticklen*panel_aspect_ratio
    fig_xsize = 5.
    fig_ysize = panel_aspect_ratio*fig_xsize*nevent
    fig_xsize += fig_xsize*panel_aspect_ratio*1    ; add to the right a squre space.
    file = join_path([project.plot_dir,'fig_all_'+data_type+'2.pdf'])
    if keyword_set(test) then file = test

    sgopen, file, xsize=fig_xsize, ysize=fig_ysize
    poss = sgcalcpos(nevent,2, ypad=3, tmargin=5, bmargin=2.5, xchsz=xchsz,ychsz=ychsz, xpans=[5,1], xpad=8)

    panel_infos = hash()
    foreach event_id, sorted_events, ii do begin
        panel_infos[event_id] = dictionary($
            'position', poss[*,0,ii], $
            'orbit_position', poss[*,1,ii])
    endforeach

;---Overall label.
    tpos = poss[*,0,0]
    tx =  lmargin*xchsz
    ty = tpos[3]-ychsz*0.2
    xyouts, tx,ty,/normal,alignment=0, 'Tilt angle '+str_delta+str_theta+' time lags applied'
    tpos = poss[*,1,0]
    tx = 0.5*(tpos[0]+tpos[2])-xchsz*1.5
    ty = tpos[3]+ychsz*3
    xyouts, tx,ty,/normal,alignment=0.5, 'SM X-Y plane'


;---Color bar.
    tpos = poss[*,0,0]
    tpos = tpos[[0,3,2,3]]
    tpos[1] += ychsz*1.0
    tpos[3] = tpos[1]+ychsz*0.5

    device, decomposed=1
    ncolor = n_elements(colors)

    cbpos = tpos & cbpos[2] = mean(tpos[[0,2]])
    sgtv, reverse(ncolor-findgen(ncolor/2)) # (intarr(2)+1), position=cbpos, /resize, ct=dusk_ct
    cbpos = tpos & cbpos[0] = mean(tpos[[0,2]])
    sgtv, (ncolor-findgen(ncolor/2)) # (intarr(2)+1), position=cbpos, /resize, ct=dawn_ct
    xrange = [-1,1]*max_mlt
    plot, xrange, [-1,1], $
        xstyle=1, xtickformat='(A1)', xticks=1, xminor=0, $
        ystyle=1, ytickformat='(A1)', yticks=1, yminor=0, $
        position=tpos, /noerase, /nodata
    xminor = 2
    xtickv = make_bins(xrange, xminor)
    xticks = n_elements(xtickv)-1
    axis, xaxis=1, xrange=xrange, xstyle=1, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=-0.2, xtitle='MLT (hr)'



;---Panels.
    fig_id = 0
    foreach event_id, sorted_events do begin
        event = project.events[event_id]

        ; Load data for the current event.
        data_file = event.file
        if file_test(data_file) eq 0 then azim_prop_gen_survey_plot, /save_data
        if file_test(data_file) eq 0 then message, 'Something is wrong when loading the data ...'
        tplot_restore, file = data_file

        ; Prepare the x-axis info.
        event_time_range = event.time_range
        tpos = panel_infos[event_id].position
        probes = event.probes
        nprobe = n_elements(probes)
        probe_mlts = dblarr(nprobe)

        sorted_probes = event.sorted_probes
        xrange = event_time_range
        ; xrange = event.ref_time+[-1,1]*60*40 ; works but not significantly better than before.
        xstep = 30*60
        xtickv = smkarthm(xrange[0],xrange[1],3,'n')
        xticks = n_elements(xtickv)-1
        xminor = 0
        xtickn = strarr(xticks+1)
        for ii=0, xticks do xtickn[ii] = time_string(xtickv[ii],tformat='YYYY-MM-DD/hh:mm')
        for ii=0, strlen(xtickn[0]) do begin
            xtickn[0] = ' '+xtickn[0]
            xtickn[2] = xtickn[2]+' '
        endfor


        ; Prepare the y-axis info.
        yrange = []
        offsets = []
        foreach probe, sorted_probes do begin
            var = probe+'_'+data_type
            get_data, probe+'_'+data_type, times, data
            index = where_pro(times, event_time_range)
            data = data[index]
            offset = mean(data[where_pro(times, event_time_range[0]+[10,40]*60)])
            data = data-offset
            offsets = [offsets,offset]
            yrange = [yrange, minmax(data)]
        endforeach
        yrange = minmax(yrange)
        ytickv = make_bins(yrange, data_step)
        yrange = minmax(ytickv)
        yticks = n_elements(ytickv)-1

        plot, xrange, yrange, /noerase, /nodata, position=tpos, $
            xrange=xrange, xstyle=5, $
            yrange=yrange, ystyle=5
        oplot, event.ref_time+[0,0], yrange, linestyle=1

        ; Add y-labels.
        tx = tpos[0]-xchsz*1
        tys = mean(tpos[[1,3]])+[-1,1]*(tpos[3]-tpos[1])*0.25
        plots, tx+[0,0], tys, /normal
        foreach ty,tys do plots, tx+[-1,1]*ychsz*0.1/fig_xsize*fig_ysize, ty, /normal
        tx = tx-xchsz*1
        ty = mean(tys)-ychsz*0.3
        xyouts, tx,ty,/normal, alignment=1, string((yrange[1]-yrange[0])/2,format='(I02)')+' '+var_unit

        ; Add x-labels.
        ty = tpos[1]
        plots, tpos[[0,2]], ty+[0,0], /normal
        foreach tx, tpos[[0,2]] do plots, tx+[0,0],ty+[-1,1]*ychsz*0.1, /normal
        ty = tpos[1]-ychsz
        tx = tpos[0]
        msg = time_string(xrange[0],tformat='YYYY-MM-DD/hh:mm')
        xyouts, tx,ty,/normal, alignment=0.2, msg

        ; end time.
        tx = tpos[2]
        msg = time_string(xrange[1],tformat='YYYY-MM-DD/hh:mm')
        xyouts, tx,ty,/normal, alignment=0.8, msg

        fig_id += 1
        tx =  lmargin*xchsz
        xyouts, tx,ty,/normal, 'Event '+string(fig_id,format='(I0)')

        ; Add data and s/c label.
        ty = tpos[1]-ychsz
        probes = sorted_probes[sort(sorted_probes)]
        nprobe = n_elements(probes)
        strlength = nprobe*3
        tx = (tpos[0]+tpos[2])*0.5-xchsz*0.5*strlength
        foreach probe, sorted_probes, ii do begin
            ref_mlt = event[probe].ref_mlt
            tmp = min(mlt_values-ref_mlt, index, /absolute)
            color = mlt_colors[index]
            var = probe+'_'+data_type
            get_data, probe+'_'+data_type, times, data
            ;data = deriv(data)/5*10
            ;data = smooth(data,10)
            time_lag = event[probe].time_lag
            times = times+time_lag
            data = data-offsets[ii]
            oplot, times, data, color=color

            sc_short = strupcase(project[probe].short_name)
            xyouts, tx,ty,/normal, sc_short, color=color
            tx += strlen(sc_short)*xchsz
        endforeach

        txs = [xchsz*2,1-xchsz*4]
        tys = tpos[1]-ychsz*1.2+[0,0]
        plots, txs,tys,/normal, thick=0.2;, linestyle=1

        ; Add duration.
        tx = tpos[2]+xchsz*0.5
        ty = tpos[3]-ychsz*0.2
        duration = total(event_time_range*[-1,1])/3600.
        xyouts, tx,ty,/normal,alignment=1, 'dT = '+string(duration,format='(F3.1)')+' hr'


        ; Add ramp duration and scale.
        tilt_time_range = event.tilt_time_range
        duration = event.tilt_duration/60    ; in min.
        scale = event.tilt_scale             ; in Re.
        txs = tilt_time_range
        tmp = convert_coord(txs,[0,0], /data, /to_normal)
        txs = reform(tmp[0,*])
        tys = tpos[1]+ychsz*0.3+[0,0]
        red = sgcolor('red')
        plots, txs,tys,/normal, color=red
        foreach tx, txs do plots, tx+[0,0],tys[0]+[-1,1]*ychsz*0.2,/normal, color=red
        msg = string(duration,format='(I0)')+' min, '+string(scale,format='(F3.1)')+' Re'
        ;tx = mean(txs)
        ;ty = tys[0]+ychsz*0.1
        ;xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size, color=red
        tx = txs[1]+xchsz*1
        ty = tys[0]-ychsz*0.1
        xyouts, tx,ty,/normal, alignment=0, msg, charsize=label_size, color=red


    ;---Add orbit.
        tpos = panel_infos[event_id].orbit_position
        tpos[[1,3]] += [-1,1.5]*ychsz
        tpos[[0,2]] += [-3,5]*xchsz
        xrange = [15,-15]
        xticks = 2
        xminor = 5
        xticklen = -0.01
        yrange = [15,-15]
        yticks = 2
        yminor = 5
        yticklen = -0.01
        plot, xrange, yrange, /nodata, /noerase, $
            xstyle=5, xticks=xticks, xticklen=xticklen, xminor=xminor, xrange=xrange, $
            ystyle=5, yticks=yticks, yticklen=yticklen, yminor=yminor, yrange=yrange, $
            position=tpos, /isotropic, charsize=label_size

        ; Add earth.
        tmp = 50
        tmp = findgen(tmp)/(tmp-1)*2*!dpi
        xs = cos(tmp)
        ys = sin(tmp)
        polyfill, (xs*15)<0, ys*15, color=sgcolor('silver')
        plots, xs, ys
        foreach r, [5,10,15] do begin
            oplot, xs*r, ys*r, linestyle=1
            msg = string(r,format='(I0)')+' '
            tmp = -0*rad
            tx = cos(tmp)*r
            ty = sin(tmp)*r
            xyouts, tx,ty,/data, alignment=0.5, msg, charsize=label_size
        endforeach

        foreach probe, sorted_probes do begin
            ref_mlt = event[probe].ref_mlt
            tmp = min(mlt_values-ref_mlt, index, /absolute)
            color = mlt_colors[index]

            ref_rsm = event[probe].ref_rsm
            tx = ref_rsm[0]
            ty = ref_rsm[1]
            plots, tx,ty,/data, psym=1, symsize=label_size, color=color

            ;tmp = convert_coord([tx,ty], /data, /to_normal)
            ;tx = tmp[0]+xchsz*label_size*0.5
            ;ty = tmp[1]-ychsz*label_size*0.25
            ;short_name = strupcase(project[probe].short_name)
            ;xyouts, tx,ty,/normal, strupcase(short_name), color=color, charsize=label_size
        endforeach
        if keyword_set(test) then stop
    endforeach

    if keyword_set(test) then stop
    sgclose

end
