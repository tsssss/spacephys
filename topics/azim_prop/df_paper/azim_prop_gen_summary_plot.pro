;+
; Generate summary plot per storm, and generate 3 plots:
;   1. the tilt angle,
;   2. tilt angle difference from the model,
;   3. |B| difference from the model.
;
; Set start_after to start from the next storm. This is necessary since the
; program can be interrupted during downloading data from servers.
;
; stop_at is used for debugging, should be commented out in normal usage.
;-

pro azim_prop_gen_summary_plot, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()

;---Settings.
    summary_plot_setting = dictionary()
    summary_plot_setting.movie_fps = 10d         ; frame per sec.
    summary_plot_setting.movie_size = 800        ; pixel.
    summary_plot_setting.movie_time_step = 600.  ; 10 min.
    summary_plot_setting.model = 't89'
    summary_plot_setting.orbit_data_rate = 60.   ; sec.
    summary_plot_setting.log_file = join_path([project.data_dir,'gen_summary_plot.log'])
    if file_test(summary_plot_setting.log_file) eq 1 then file_delete, summary_plot_setting.log_file
    stouch, summary_plot_setting.log_file
    margins = [8,5,5,3]

    orbit = dictionary()
    orbit.position = [0.15,0.10,0.95,0.90]
    orbit.xrange = [15,-15]
    orbit.xminor = 5
    orbit.xtickv = make_bins(orbit.xrange,orbit.xminor)
    orbit.xticks = n_elements(orbit.xtickv)
    orbit.xticklen = -0.01
    orbit.xtitle = 'SM X (Re)'
    orbit.yrange = [15,-15]
    orbit.yminor = 5
    orbit.ytickv = make_bins(orbit.yrange,orbit.yminor)
    orbit.yticks = n_elements(orbit.ytickv)
    orbit.yticklen = -0.01
    orbit.ytitle = 'SM Y (Re)'

    constant = project.constant
    console = -1
    tab = '    '

    tplot_options, 'labflag', -1
    tplot_options, 'ynozero'
    tplot_options, 'ystyle', 1
    tplot_options, 'yticklen', -0.005
    tplot_options, 'xticklen', -0.03

    foreach dir, [project.root_dir,project.data_dir,project.plot_dir] do $
        if file_test(dir,/directory) eq 0 then file_mkdir, dir


;---For each storm,
;   1. generate the movie of spacecraft position in the SM X-Y plane,
;   2. a line plot of B tilt in the SM coordinate.

    if ~project.haskey('candidate_list') then azim_prop_filter_by_apogee, project
    event_list = project.candidate_list
    foreach time_range, event_list do begin
    ;---Settings.
        event_id = time_string(mean(time_range), tformat='YYYY_MMDD_hh')
        if n_elements(start_after) ne 0 then if time_double(event_id,tformat='YYYY_MMDD_hh') le start_after then continue

    ;---Load Dst and AE.
        omni = sread_omni(time_range)
        times = sfmepoch(omni.epoch,'unix')
        store_data, 'dst', times, omni.sym_h, limits={ytitle:'(nT)',labels:'Dst'}
        store_data, 'ae', times, omni.ae_index, limits={ytitle:'(nT)',labels:'AE'}

    ;---Clear memory.
        probes = project.probes
        nprobe = n_elements(probes)
        vars = ['r_gsm','b_gsm','r_sm','b_tilt','b_mag','bmod_tilt','bmod_mag','db_tilt','db_mag']
        foreach var, vars do store_data, probes+'_'+var, /delete


    ;---Load s/c location, and B field.
        foreach probe, probes do begin
            mission = project[probe]
            call_procedure, mission.routine_name+'_read_bfield', time_range, probe=mission.probe, errmsg=errmsg
            call_procedure, mission.routine_name+'_read_orbit', time_range, probe=mission.probe
        endforeach


    ;---Filter out the probes that both have position and B field data.
        probe_flags = bytarr(nprobe)
        foreach probe, probes, ii do begin
            if tnames(probe+'_r_gsm') eq '' then continue
            if tnames(probe+'_b_gsm') eq '' then continue
            probe_flags[ii] = 1
        endforeach
        index = where(probe_flags eq 1, nprobe)
        if nprobe eq 0 then continue    ; No data, skip this storm, although very unlikely.
        probes = probes[index]
        lprmsg, event_id+tab+strjoin(time_string(time_range),tab)+tab+strjoin(probes,','), summary_plot_setting.log_file


    ;---Caclulate xxx_r_sm, xxx_b_tilt.
        foreach probe, probes do begin
            mission = resolve_probe(probe)

            times = make_bins(time_range,summary_plot_setting.orbit_data_rate)
            rgsm = get_var_data(probe+'_r_gsm', at=times)
            rsm = cotran(rgsm, times, 'gsm2sm')
            ;store_data, probe+'_r_gsm', times, rgsm
            store_data, probe+'_r_sm', times, rsm

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
            store_data, probe+'_bmod_gsm_'+summary_plot_setting.model, times, bmodgsm


            get_data, probe+'_b_gsm', times, bgsm
            bsm = cotran(bgsm, times, 'gsm2sm')
            b_tilt = atan(bsm[*,2],sqrt(bsm[*,0]^2+bsm[*,1]^2))*constant.deg
            sdespike, times, b_tilt
            store_data, probe+'_b_tilt', times, b_tilt, limits={$
                ytitle: '(deg)', labels:strupcase(mission.short_name+mission.probe)}
            b_mag = snorm(bgsm)
            sdespike, times, b_mag
            store_data, probe+'_b_mag', times, b_mag, limits={$
                ytitle: '(nT)', labels:strupcase(mission.short_name+mission.probe)}

            bmodgsm = get_var_data(probe+'_bmod_gsm_'+summary_plot_setting.model, at=times)
            bmodsm = cotran(bmodgsm, times, 'gsm2sm')
            bmod_tilt = atan(bmodsm[*,2],sqrt(bmodsm[*,0]^2+bmodsm[*,1]^2))*constant.deg
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

    ;---Fill in the probes that do not have data?
        ; decided to not do this.


    ;---Plot B tilt angle.
        ;plot_vars = ['b_tilt','db_tilt','db_mag']
        ;plot_titles = ['SM B tilt angle','SM B-B_mod tilt angle','|B|-|B_mod|']
        plot_vars = ['db_tilt']
        plot_titles = ['SM B-B_'+strupcase(summary_plot_setting.model)+' tilt angle']
        foreach var, plot_vars, ii do begin
            if keyword_set(test) then begin
                file = 0
            endif else begin
                file = project.plot_dir+'/'+event_id+'/'+event_id+'_'+var+'.pdf'
            endelse
            temp_dir = file_dirname(file)
            if file_test(temp_dir,/directory) eq 0 then file_mkdir, temp_dir

            fig_ysize = (2+nprobe)*1
            fig_xsize = total(time_range*[-1,1])/constant.secofday*8
            sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch
            vars = ['dst','ae',probes+'_'+var]
            nvar = n_elements(vars)
            poss = sgcalcpos(nvar, margins=margins, xchsz=xchsz,ychsz=ychsz)
            tplot, vars, trange=time_range, /novtitle, position=poss
            tpos = poss[*,0]
            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            xyouts, tx,ty,/normal, event_id+', '+plot_titles[ii]
            sgclose
        endforeach


    ;---Make a movie of the spacecraft position.
        movfn = project.plot_dir+'/'+event_id+'/'+event_id+'_sm_orbit.mp4'
        temp_dir = file_dirname(movfn)
        if file_test(temp_dir,/directory) eq 0 then file_mkdir, temp_dir

        ovid = idlffvideowrite(movfn)
        vidstream = ovid.AddVideoStream(summary_plot_setting.movie_size, summary_plot_setting.movie_size, summary_plot_setting.movie_fps)
        temp_file = project.root_dir+'/tmp.png'

        times = make_bins(time_range, summary_plot_setting.movie_time_step)

        foreach time, times, i do begin
            printf, console, time_string(time)
            sgopen, temp_file, xsize=summary_plot_setting.movie_size, ysize=summary_plot_setting.movie_size
            tpos = sgcalcpos(1, position=orbit.position, xchsz=xchsz,ychsz=ychsz)

            ; Draw box.
            plot, orbit.xrange, orbit.yrange, $
                xstyle=1, ystyle=1, /nodata, _extra=orbit.tostruct()

            ; Add label.
            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            xyouts, tx,ty,/normal, time_string(time, tformat='YYYY-MM-DD/hh:mm:ss')
            tx = tpos[0]+xchsz*15
            xyouts, tx,ty,/normal, 'AE: '+string(get_var_data('ae',at=time),format='(I4)')+' nT'
            foreach probe, project.probes, jj do begin
                tx = tpos[0]+xchsz*(25+5*jj)
                mission = project[probe]
                xyouts, tx,ty,/normal, strupcase(mission.short_name), color=mission.color
            endforeach

            ; Add earth and lines.
            tmp = 50
            tmp = findgen(tmp)/(tmp-1)*2*!dpi
            xs = cos(tmp)
            ys = sin(tmp)
            polyfill, xs<0, ys, /line_fill, orientation=45
            plots, xs, ys
            foreach r, [5,10,15] do oplot, xs*r, ys*r, linestyle=1

            plots, orbit.xrange, [0,0], linestyle=1
            plots, [0,0], orbit.yrange, linestyle=1

            ; Add position of the probe.
            foreach probe, project.probes, jj do begin
                var = probe+'_r_sm'
                if tnames(var) eq '' then continue
                rsm = get_var_data(var, at=time)
                mission = project[probe]
                plots, rsm[0],rsm[1],/data, color=mission.color, psym=mission.psym
            endforeach

            sgclose
            read_png, temp_file, temp_image
            time = ovid.put(vidstream, temp_image)
        endforeach
        ovid.Cleanup
        file_delete, temp_file
    endforeach


end
