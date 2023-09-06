;+
; Plot all the data for all events, each panel has all spacecraft data shifted by the proper time lags.
;-

pro azim_prop_plot_all_data, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    test = 0
    
    mlt_colors = list()
    mlt_values = list()
    ncolor = 140
    top_color = 255
    foreach key, ['dawn','dusk'] do begin
        ct = (key eq 'dawn')? 49: 50
        mlt_range = (key eq 'dawn')? [0.,7]: [0.,-7]
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
    data_step = 2
    lmargin = 2

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
    file = join_path([project.plot_dir,'fig_all_'+data_type+'.pdf'])
    if keyword_set(test) then file = test

    sgopen, file, xsize=fig_xsize, ysize=fig_ysize
    poss = sgcalcpos(nevent, ypad=3, tmargin=2.5, bmargin=5, xchsz=xchsz,ychsz=ychsz)
    panel_infos = hash()
    foreach event_id, sorted_events, ii do begin
        panel_infos[event_id] = dictionary($
            'position', poss[*,ii])
    endforeach

    tpos = poss[*,0]
    tx = lmargin*xchsz
    ty = tpos[3]+ychsz*1
    xyouts, tx,ty,/normal,alignment=0, str_delta+str_theta+' of all events, after time lags are applied'

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
            ;xrange=xrange, xstyle=9, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickformat='(A1)', $
            xrange=xrange, xstyle=5, $
            yrange=yrange, ystyle=5
        oplot, event.ref_time+[0,0], yrange, color=color, linestyle=1

        ; Add y-labels.
        tx = tpos[0]-xchsz*1
        tys = mean(tpos[[1,3]])+[-1,1]*(tpos[3]-tpos[1])*0.25
        plots, tx+[0,0], tys, /normal
        foreach ty,tys do plots, tx+[-1,1]*ychsz*0.1/fig_xsize*fig_ysize, ty, /normal
        tx = tx-xchsz*1
        ty = mean(tys)-ychsz*0.3
        xyouts, tx,ty,/normal, alignment=1, string((yrange[1]-yrange[0])/2,format='(I02)')+' deg'

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
    endforeach


    if keyword_set(test) then stop
    sgclose

end
