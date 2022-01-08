;+
; Try to get the width of the ramp.
;-

pro azim_prop_calc_ramp_info, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    if ~project.haskey('done_calc_2d_vel') then azim_prop_calc_2d_vel, project


    test = 0


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
    panel_aspect_ratio = 0.25


    mlt_colors = list()
    mlt_values = list()
    ncolor = 140
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

;---Loop through each event.
    ; Sort event by date.
    events = project.events
    sorted_events = (events.keys()).toarray()
    index = sort(time_double(sorted_events,tformat='YYYY_MM_DD'))
    sorted_events = sorted_events[index]

    ; Common data rate.
    data_rate = project.common_data_rate
    ref_time_range = [-1,1]*30*60


    foreach event_id, sorted_events do begin
        event = project.events[event_id]
        omega_azim = event.omega_azim

        ; Load data for the current event.
        azim_prop_load_data, project=project, event_id=event_id
        sorted_probes = event.sorted_probes
        nprobe = n_elements(sorted_probes)

        ; Time range and common times.
        event_time_range = event.time_range
        common_times = smkarthm(event_time_range[0],event_time_range[1],data_rate,'dx')

        if keyword_set(test) then begin
            file = test
        endif else begin
            file = join_path([project.plot_dir,'tilt_ramp','fig_'+event_id+'_calc_ramp_info.pdf'])
        endelse
        sgopen, file, xsize=5, ysize=5*(nprobe+1)*panel_aspect_ratio
        poss = sgcalcpos(nprobe+1, ypad=4, tmargin=4,bmargin=4,rmargin=8, xchsz=xchsz,ychsz=ychsz)
        tx = 0.5
        ty = 1-ychsz*2
        xyouts, tx,ty,/normal, alignment=0.5, 'Event '+event_id

        ; Trim to the time around ref_time.
        foreach probe, sorted_probes, ii do begin
            ref_time = event[probe].ref_time
            search_time_range = ref_time_range+ref_time
            index = lazy_where(common_times, search_time_range)
            times = common_times[index]
            data = get_var_data(probe+'_db_tilt', at=times)
            sdespike, times, data, width=n_elements(times)/10

            var = probe+'_db_tilt_tmp'
            store_data, var, times, data, limits={$
                ytitle:'('+var_unit+')', labels:strupcase(project[probe].short_name)}

            ; Find the time of the ramp and the time range.
            max_val = max(data, max_id)
            min_val = min(data[0:max_id], min_id)
            del_val = 0.25*(max_val-min_val)
            mean_val = 0.5*(max_val+min_val)
            half_min_val = mean_val-del_val
            half_max_val = mean_val+del_val

            tdata = data[min_id:max_id]
            ttimes = times[min_id:max_id]
            index = where(tdata le half_min_val)
            min_time = ttimes[index[-1]]
            min_val = tdata[index[-1]]
            index = where(tdata ge half_max_val and ttimes gt min_time)
            max_time = ttimes[index[0]]
            max_val = tdata[index[0]]

            event[probe].ramp_time_range = [min_time, max_time]
            event[probe].ramp_value_range = [min_val, max_val]
            event[probe].ramp_duration = max_time-min_time
            event[probe].ramp_value_change = max_val-min_val

            ref_mlt = event[probe].ref_mlt
            ref_dis = snorm(event[probe].ref_rsm[0:1])
            tpos = poss[*,ii]

            xrange = search_time_range
            xtickv = smkarthm(xrange[0],xrange[1],3,'n')
            xticks = n_elements(xtickv)-1
            xtickn = time_string(xtickv,tformat='hh:mm')
            xminor = 10
            xticklen = -0.04

            yrange = minmax(make_bins(data,2))
            yticks = 2
            yminor = 5
            yticklen = -0.01
            plot, times, data, $
                xstyle=1, xticks=xticks, xrange=xrange, xtickv=xtickv, xtickname=xtickn, xticklen=xticklen, xminor=xminor, $
                ystyle=1, yticks=yticks, yrange=yrange, yticklen=yticklen, yminor=yminor, ytitle='('+var_unit+')', $
                /noerase, position=tpos
            foreach tx, [min_time, max_time] do plots, tx+[0,0], yrange
            plots, ref_time+[0,0], yrange, linestyle=1

            tx = tpos[0]
            ty = tpos[3]+ychsz*0.3
            duration = abs(max_time-min_time)   ; in sec.
            scale = duration*abs(omega_azim)*coef2kms/re
            msg = 'Ramp duration: '+$
                strtrim(string(duration/60,format='(F5.1)'),2)+' min, jump: '+$
                strtrim(string(max_val-min_val,format='(F5.1)'),2)+' '+var_unit+', scale: '+$
                strtrim(string(scale,format='(F5.1)'),2)+' Re'
            xyouts, tx,ty,/normal, msg

            tx = tpos[2]+xchsz*0.5
            ty = (tpos[1]+tpos[3])*0.5+ychsz*1
            xyouts, tx,ty,/normal, strupcase(project[probe].short_name)
            ty = ty-ychsz*1
            xyouts, tx,ty,/normal, strtrim(string(ref_mlt,format='(F5.1)'),2)+' MLT'
            ty = ty-ychsz*1
            xyouts, tx,ty,/normal, strtrim(string(ref_dis,format='(F5.1)'),2)+' Re'
        endforeach


    ;---Last panel, summarize all spacecraft.
        tpos = poss[*,nprobe]

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

        yrange = []
        offsets = []
        foreach probe, sorted_probes do begin
            var = probe+'_'+data_type
            get_data, probe+'_'+data_type, times, data
            index = lazy_where(times, event_time_range)
            data = data[index]
            offset = mean(data[lazy_where(times, event_time_range[0]+[10,40]*60)])
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
        oplot, event.ref_time+[0,0], yrange, color=color, linestyle=1

        ; Add y-labels.
        tx = tpos[0]-xchsz*1
        tys = mean(tpos[[1,3]])+[-1,1]*(tpos[3]-tpos[1])*0.25
        plots, tx+[0,0], tys, /normal
        foreach ty,tys do plots, tx+[-1,1]*xchsz*0.2, ty, /normal
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

        ; Add data and s/c label.
        ty = tpos[1]-ychsz
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


        ; Add duration.
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]+ychsz*0.3
        duration = total(event_time_range*[-1,1])/3600.
        xyouts, tx,ty,/normal,alignment=1, 'dT = '+string(duration,format='(F3.1)')+' hr'


    ;---Figure out the overall start and end times.
        time_ranges = fltarr(2,nprobe)
        foreach probe, sorted_probes, ii do begin
            ref_time = event[probe].ref_time
            time_ranges[*,ii] = event[probe].ramp_time_range-ref_time
        endforeach
        ref_time = event.ref_time

        ; Test 1.
        error = mean([stddev(time_ranges[0,*]),stddev(time_ranges[1,*])])
        min_time = mean(time_ranges[0,*])-0.5*error
        max_time = mean(time_ranges[1,*])+0.5*error

        ; Test 2.
        mid_time = mean((time_ranges[1,*]+time_ranges[0,*])*0.5)
        durations = time_ranges[1,*]-time_ranges[0,*]
        duration = mean(durations)
        error = stddev(durations)
        min_time = mid_time-duration*0.5-error*0.5
        max_time = mid_time+duration*0.5+error*0.5

;        ; center.
;        txs = [min_time,max_time]+ref_time
;        tmp = convert_coord(txs,[0,0], /data, /to_normal)
;        txs = reform(tmp[0,*])
;        tys = tpos[1]+ychsz*0.5+[0,0]
;        red = sgcolor('red')
;        plots, txs,tys,/normal, color=red
;        foreach tx, txs do plots, tx+[0,0],tys[0]+[-1,1]*ychsz*0.2,/normal, color=red


        tilt_time_range = [min_time,max_time]+ref_time
        txs = tilt_time_range
        tmp = convert_coord(txs,[0,0], /data, /to_normal)
        txs = reform(tmp[0,*])
        tys = tpos[1]+ychsz*0.25+[0,0]
        red = sgcolor('red')
        plots, txs,tys,/normal, color=red
        foreach tx, txs do plots, tx+[0,0],tys[0]+[-1,1]*ychsz*0.2,/normal, color=red


        duration = abs(tilt_time_range[1]-tilt_time_range[0])   ; in sec.
        scale = duration*abs(omega_azim)*coef2kms/re
        msg = 'Ramp duration: '+$
            strtrim(string(duration/60,format='(F5.1)'),2)+' min, scale: '+$
            strtrim(string(scale,format='(F5.1)'),2)+' Re'
        tx = tpos[0]
        ty = tpos[3]+ychsz*0.3
        xyouts, tx,ty,/normal, msg

        event.tilt_time_range = tilt_time_range
        event.tilt_duration = duration
        event.tilt_scale = scale

        project.events[event_id] = event
        if keyword_set(test) then stop
        sgclose
    endforeach


    scales = []
    foreach event_id, sorted_events do scales = [scales,project.events[event_id].tilt_scale]


    project.done_calc_ramp_info = 1
    if ~keyword_set(test) then azim_prop_update_project, project
end
