;+
; Load B data and calculate tilt angle and plot the position of spacecraft.
;-
;


pro azim_prop_gen_survey_plot, project, event_list=event_list, save_data=save_data, event_ids=event_ids

;---Check input.
    if n_elements(project) eq 0 then project = azim_prop_load_project()
    if n_elements(event_list) eq 0 then event_list = join_path([project.data_dir,'event_list.txt'])
    if file_test(event_list) eq 0 then message, 'Event list does not exist ...'


;---Constants.
    deg = project.constant.deg
    model = 't89'
    trace_dir = -1
    orbit_data_rate = 60.   ; 1 min.
    test = 1


;---Read text file into a dictionary.
    lines = read_all_lines(event_list)
    events = hash()
    foreach line, lines do begin
        tinfo = strsplit(line,' ',/extract)
        event_id = tinfo[0]
        time_range = time_double(tinfo[1:2])
        probes = strsplit(tinfo[3],',',/extract)
        note = strjoin(tinfo[4:*],' ')
        events[event_id] = dictionary($
            'id', event_id, $
            'time_range', time_range, $
            'probes', probes, $
            'note', note)
    endforeach
    project.events = events
    azim_prop_update_project, project
    if n_elements(event_ids) eq 0 then event_ids = events.keys()



    foreach event_id, event_ids do begin
        time_range = events[event_id].time_range
        probes = events[event_id].probes
        nprobe = n_elements(probes)

    ;---Clear memory.
        vars = ['r_gsm','b_gsm','r_sm','b_tilt','b_mag','bmod_tilt','bmod_mag','db_tilt','db_mag']
        foreach var, vars do store_data, '*_'+var, /delete

    ;---Load s/c location, and B field.
        foreach probe, probes do begin
            mission = project[probe]
            call_procedure, mission.routine_name+'_read_bfield', time_range, probe=mission.probe, errmsg=errmsg
            call_procedure, mission.routine_name+'_read_orbit', time_range, probe=mission.probe
        endforeach

    ;---Caclulate xxx_r_sm, xxx_mlt, xxx_b_tilt.
        foreach probe, probes do begin
            mission = resolve_probe(probe)

            ; Calculate r_sm, and MLT.
            times = make_bins(time_range,orbit_data_rate)
            rvar = project[probe].prefix+'_r_gsm'
            rgsm = get_var_data(rvar, at=times)
            rsm = cotran(rgsm, times, 'gsm2sm')
            ;store_data, probe+'_r_gsm', times, rgsm
            store_data, probe+'_r_sm', times, rsm

            rmag = cotran(rgsm, times, 'gsm2mag')
            mlon = atan(rmag[*,1],rmag[*,0])*deg
            mlt = mlon2mlt(mlon, times)
            store_data, project[probe].prefix+'_mlt', times, mlt, limits={ytitle:'(deg)'}

            ; Calculate the model field.
            ntime = n_elements(times)
            bmodgsm = fltarr(ntime,3)
            for kk=0,ntime-1 do begin
                tilt = geopack_recalc(times[kk])

                ; in-situ position
                rx = rgsm[kk,0]
                ry = rgsm[kk,1]
                rz = rgsm[kk,2]

                ; in-situ B field.
                geopack_igrf_gsm, rx,ry,rz, bx,by,bz
                geopack_t89, 2, rx,ry,rz, dbx,dby,dbz
                bmodgsm[kk,*] = [bx,by,bz]+[dbx,dby,dbz]
            endfor
            store_data, probe+'_bmod_gsm_'+model, times, bmodgsm


            get_data, probe+'_b_gsm', times, bgsm
            bsm = cotran(bgsm, times, 'gsm2sm')
            b_tilt = atan(bsm[*,2],sqrt(bsm[*,0]^2+bsm[*,1]^2))*deg
            sdespike, times, b_tilt
            store_data, probe+'_b_tilt', times, b_tilt, limits={$
                ytitle: '(deg)', labels:strupcase(mission.short_name+mission.probe)}
            b_mag = snorm(bgsm)
            sdespike, times, b_mag
            store_data, probe+'_b_mag', times, b_mag, limits={$
                ytitle: '(nT)', labels:strupcase(mission.short_name+mission.probe)}

            bmodgsm = get_var_data(probe+'_bmod_gsm_'+model, at=times)
            bmodsm = cotran(bmodgsm, times, 'gsm2sm')
            bmod_tilt = atan(bmodsm[*,2],sqrt(bmodsm[*,0]^2+bmodsm[*,1]^2))*deg
            store_data, probe+'_bmod_tilt', times, bmod_tilt, limits={$
                ytitle: '(deg)', labels:strupcase(mission.short_name+mission.probe)}
            bmod_mag = snorm(bmodgsm)
            store_data, probe+'_bmod_mag', times, bmod_mag, limits={$
                ytitle: '(nT)', labels:strupcase(mission.short_name+mission.probe)}

            dis = snorm(get_var_data(probe+'_r_gsm', at=times))
            db_tilt = b_tilt-bmod_tilt
            index = where(dis ge 4 and db_tilt gt -100)
            yrange = minmax(db_tilt[index])
            store_data, probe+'_db_tilt', times, db_tilt, limits={$
                ytitle: '(deg)', labels:strupcase(mission.short_name+mission.probe), yrange: yrange}
            db_mag = b_mag-bmod_mag
            index = where((db_mag le 300) and (dis ge 4))
            yrange = minmax(db_mag[index])
            store_data, probe+'_db_mag', times, db_mag, limits={$
                ytitle: '(nT)', labels:strupcase(mission.short_name+mission.probe), yrange: yrange}
        endforeach


    ;---Save data.
        vars = ['r_gsm','b_gsm','bmod_gsm_'+model,'mlt','db_tilt','db_mag']
        saved_vars = []
        foreach probe, probes do saved_vars = [saved_vars,probe+'_'+vars]
        if keyword_set(save_data) then tplot_save, saved_vars, filename=join_path([project.data_dir,event_id+'_basic_data.tplot'])



    ;---Prepare to plot.
        fig_xsize = 4   ; inch.
        lmargin = 8
        rmargin = 6
        ticklen = -0.01
        sgopen, 0, xsize=fig_xsize, ysize=fig_xsize, /inch
        tpos = sgcalcpos(1,lmargin=lmargin, rmargin=rmargin, xchsz=xchsz, ychsz=ychsz)
        ypad = 0.5  ; y-pad space between panels.
        sgclose, /wdelete

        color_start = 50
        color_end = 250
        color_table = 40
        colors = round(smkarthm(color_start,color_end, nprobe, 'n'))
        for ii=0, nprobe-1 do colors[ii] = sgcolor(colors[ii], ct=color_table)

        ; TIlt angle panels.
        tilt_panel = dictionary()
        tilt_panel.aspect_ratio = 0.2
        tilt_panel_xsize = tpos[2]-tpos[0]
        tilt_panel_ysize = tilt_panel_xsize*tilt_panel.aspect_ratio*nprobe+ychsz*ypad*(nprobe-1)


        panel_vars = probes+'_b_tilt'
        panel_title = 'Tilt angle in SM'


        yticks = 2
        yminor = 5
        tilt_var_labels = strarr(nprobe)
        foreach tvar, panel_vars, ii do begin
            data = get_var_data(tvar)
            index = lazy_where(data-mean(data), [-1,1]*5*stddev(data))
            yrange = minmax(data[index])
            yrange = yrange-(yrange mod yticks)+[0,2]
            options, tvar, 'yrange', yrange
            options, tvar, 'yticks', yticks
            options, tvar, 'yminor', yminor
            options, tvar, 'colors', sgcolor('black')
            options, tvar, 'yticklen', ticklen
            options, tvar, 'xticklen', ticklen/tilt_panel.aspect_ratio
            options, tvar, 'labflag', -1
            options, tvar, 'ystyle', 1
            tilt_var_labels[ii] = get_setting(tvar, 'labels')
            if tilt_var_labels[ii] eq '' then tilt_var_labels[ii] = strupcase(strmid(tvar,0,strpos(tvar,'_')))
            options, tvar, 'labels', ''
        endforeach



        ; Orbit panel.
        pos_vars = probes+'_r_gsm'
        orbit_x_range = [-1,1]
        foreach tvar, pos_vars do orbit_x_range = minmax([orbit_x_range,minmax((get_var_data(tvar))[*,0])])
        orbit_y_range = [-1,1]
        foreach tvar, pos_vars do orbit_y_range = minmax([orbit_y_range,minmax((get_var_data(tvar))[*,1])])

        orbit_panel = dictionary()
        orbit_panel.xminor = 5
        orbit_panel.xrange = [ceil(max(orbit_x_range)),floor(min(orbit_x_range))]
        orbit_panel.xtickv = make_bins(orbit_panel.xrange, orbit_panel.xminor, /inner)
        orbit_panel.xticks = n_elements(orbit_panel.xtickv)-1
        orbit_panel.xtitle = 'SM X (Re)'
        orbit_panel.yminor = 5
        orbit_panel.yrange = [ceil(max(orbit_y_range)),floor(min(orbit_y_range))]
        orbit_panel.ytickv = make_bins(orbit_panel.yrange, orbit_panel.yminor, /inner)
        orbit_panel.yticks = n_elements(orbit_panel.ytickv)-1
        orbit_panel.ytitle = 'SM Y (Re)'
        orbit_panel.aspect_ratio = abs(double(orbit_panel.yrange[1]-orbit_panel.yrange[0])/(orbit_panel.xrange[1]-orbit_panel.xrange[0]))

        orbit_panel_xsize = tpos[2]-tpos[0]
        orbit_panel_ysize = orbit_panel_xsize*orbit_panel.aspect_ratio
        orbit_panel.xticklen = ticklen/orbit_panel.aspect_ratio
        orbit_panel.yticklen = ticklen

        tmargin = 2
        bmargin = 5
        fig_ysize = (tilt_panel_ysize+bmargin*2*ychsz+orbit_panel_ysize+tmargin*ychsz)*fig_xsize

    ;---Plot.
        if keyword_set(test) then begin
            file = test
            magnify = 1.2
        endif else begin
            file = join_path([project.plot_dir,'summary_plot',event_id+'_tilt_angle.pdf'])
            magnify = 1
        endelse

        sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify
        poss = sgcalcpos(2, lmargin=lmargin, rmargin=rmargin, tmargin=tmargin, bmargin=bmargin, ypans=[tilt_panel_ysize,orbit_panel_ysize], ypad=bmargin)

        tilt_panel.pos = poss[*,0]
        orbit_panel.pos = poss[*,1]

        middle_time = mean(time_range)
        psym = 6
        label_size = 0.8

        ; The tilt angle.
        poss = sgcalcpos(nprobe, position=tilt_panel.pos, ypad=ypad)
        ; Sort by MLT.
        mlts = fltarr(nprobe)
        for ii=0, nprobe-1 do mlts[ii] = get_var_data(probes[ii]+'_mlt',at=middle_time)
        index = sort(mlts)
        panel_vars = panel_vars[index]
        tilt_var_labels = tilt_var_labels[index]
        probes = probes[index]

        ;tplot, panel_vars, trange=time_range, position=poss, /novtitle
        for ii=0, nprobe-1 do begin
            tpos = poss[*,ii]
            nouttick = (ii ne nprobe-1)? 1: 0
            tplot, panel_vars[ii], trange=time_range, position=tpos, /novtitle, nouttick=nouttick, /noerase
            tx = tpos[2]+xchsz*1*label_size
            ty = 0.5*(tpos[3]+tpos[1])-ychsz*0.1*label_size
            xyouts, tx,ty,/normal, tilt_var_labels[ii], color=colors[ii], charsize=label_size
        endfor
        tpos = poss[*,0]
        tx = tpos[0]
        ty = tpos[3]+ychsz*0.3*label_size
        xyouts, tx,ty,/normal, panel_title

        ; The positions.
        plot, orbit_panel.xrange, orbit_panel.yrange, $
            xstyle=5, ystyle=5, $
            /nodata, /noerase, position=orbit_panel.pos, _extra=orbit_panel.tostruct()
        for ii=0, nprobe-1 do begin
            rsm = get_var_data(probes[ii]+'_r_sm')
            plots, rsm[*,0], rsm[*,1], color=colors[ii], linestyle=1
            rsm = get_var_data(probes[ii]+'_r_sm', at=middle_time)
            plots, rsm[0], rsm[1], color=colors[ii], psym=psym, symsize=label_size
            tmp = convert_coord(rsm[0:1], /data, /to_normal)
            tx = tmp[0]+xchsz*1*label_size
            ty = tmp[1]
            xyouts, tx,ty,/normal, tilt_var_labels[ii], color=colors[ii], charsize=label_size
        endfor

        ; Add earth and lines.
        tmp = 50
        tmp = findgen(tmp)/(tmp-1)*2*!dpi
        xs = cos(tmp)
        ys = sin(tmp)
        polyfill, xs<0, ys, /line_fill, orientation=45
        plots, xs, ys
        foreach r, [5,10] do oplot, xs*r, ys*r, linestyle=1

        plots, orbit_panel.xrange, [0,0], linestyle=1
        plots, [0,0], orbit_panel.yrange, linestyle=1

        plot, orbit_panel.xrange, orbit_panel.yrange, $
            xstyle=1, ystyle=1, $
            /nodata, /noerase, position=orbit_panel.pos, _extra=orbit_panel.tostruct()

        if keyword_set(test) then stop
        sgclose
    endforeach
end

azim_prop_gen_survey_plot, /save_data, event_ids='2016_1028_23'
end
