;+
; Plot dB tilt, time lagged, and MLT panel from the linear fit.
;-

pro azim_prop_plot_mlt_timing_only, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    test = 0

;---Plot settings.
    margins = [6,3,3,2]
    fig_xsize = 3
    aspect_ratio = 0.3
    yticklen = -0.02
    half_ychsz = 0.35
    full_ychsz = 0.7
    line_skip = 0.35
    label_size = 0.7
    xpad = 10
    ypad = 0.2
    mlt_panel = 2
    str_theta = '!9'+string(113b)+'!X'
    str_delta = '!9'+string(68b)+'!X'
    str_omega = '!9'+string(119b)+'!X'
    plot_title = str_delta+str_theta+', tilt angle variation from T89'
    psym = 6
    r_kms = 6.6
    re = project.constant.re
    rad = project.constant.rad
    coef2kms = r_kms*re*rad/60
    fig_labels = ['a','b','c','d','e','f','g','h','i','j','k']+'.'
    label_xpos = 5


    event_ids = (events.keys()).sort()
    foreach event_id, event_ids do begin
        data_file = join_path([project.data_dir,event_id+'_basic_data.tplot'])
        if file_test(data_file) eq 0 then azim_prop_gen_survey_plot, /save_data
        if file_test(data_file) eq 0 then message, 'Something is wrong when loading the data ...'
        events[event_id].file = data_file
        tplot_restore, file = data_file

        probes = events[event_id].probes

    ;---Sort probes.
        probe_mlts = list()
        foreach probe, probes do probe_mlts.add, (events[event_id])[probe].ref_mlt
        probe_mlts = probe_mlts.toarray()
        index = sort(probe_mlts)
        sorted_probes = probes[index]
        (events[event_id])[probe].sorted_probes = sorted_probes

    ;---Plot settings.
        event_time_range = events[event_id].time_range
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


        nvar = 1
        sgopen, 0, xsize=fig_xsize, ysize=fig_xsize
        tpos = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz, ypad=ypad, xpad=xpad, margins=margins)
        sgclose, /wdelete
        pan_xsize = fig_xsize*(tpos[2,0]-tpos[0,0])
        pan_ysize = pan_xsize*aspect_ratio
        fig_ysize = pan_ysize*mlt_panel+(margins[1]+margins[3])*ychsz*fig_xsize
        if keyword_set(test) then begin
            file = test
            magnify = 2
        endif else begin
            file = join_path([project.plot_dir,'time_lag','fig_'+event_id+'_plot_mlt_timing_only.pdf'])
            magnify = 1
        endelse

        thick = (size(file,/type) eq 7)? 4:2
        sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify
        poss = sgcalcpos(nvar, margins=margins, xchsz=xchsz,ychsz=ychsz, ypad=ypad, xpad=xpad)


        nprobe = n_elements(sorted_probes)
        ref_times = dblarr(nprobe)
        ref_mlts = dblarr(nprobe)
        probe_colors = lonarr(nprobe)
        foreach probe, sorted_probes, ii do begin
            ref_times[ii] = (events[event_id])[probe].ref_time
            ref_mlts[ii] = (events[event_id])[probe].ref_mlt
            probe_colors[ii] = project[probe].color
        endforeach
        yrange = minmax(ref_mlts)
        yrange = [floor(yrange[0]/2)*2,ceil(yrange[1]/2)*2]
        yticks = 2
        ytickv = smkarthm(yrange[0],yrange[1],yticks+1,'n')
        ytitle = 'MLT (hr)'
        yminor = abs(yrange[1]-yrange[0])
        xticklen = yticklen/aspect_ratio/mlt_panel
        tpos = poss[*,nvar-1]; & tpos[3] -= ychsz*1.5
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
            ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            position=tpos, /noerase, /nodata

        for ii=0, nprobe-1 do plots, ref_times[ii],ref_mlts[ii], psym=psym, symsize=label_size*0.5, color=probe_colors[ii]
        fit_result = linfit(ref_times, ref_mlts)
        ;omega_azim = fit_result[1]*15*60
        omega_azim = events[event_id].omega_azim
        txs = xrange
        tys = xrange*fit_result[1]+fit_result[0]
        oplot, txs, tys, linestyle=1
        string_dir = (omega_azim ge 0)? 'eastward': 'westward'
        alignment = (omega_azim ge 0)? 0: 1
        msg = ['|'+str_omega+'!Dazim!N| = '+sgnum2str(abs(omega_azim),nsgn=2)+' deg/min', $
            sgnum2str(round(abs(omega_azim*coef2kms)))+' km/s at '+sgnum2str(r_kms)+' Re', string_dir]
        
        rsqr = events[event_id].rsqr
        msg = [msg, '','r!U2!N = '+string(rsqr,format='(F4.2)')]
        
        tx = (omega_azim ge 0)? tpos[0]+xchsz*0.5: tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*label_size
        foreach tmp, msg, ii do $
            xyouts, tx,ty-ychsz*ii*label_size,/normal, alignment=alignment, msg[ii], charsize=label_size

        tx = tpos[0]
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,/normal, fig_labels[nvar-1]+' Rate of "drift" in MLT'
        if keyword_set(test) then stop
        sgclose
    endforeach

end
