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

;---Settings.
    start_after = time_double('2017_0422',tformat='YYYY_MMDD')
    stop_at = time_double('2015_0608',tformat='YYYY_MMDD')
    stop_at = time_double('2014_0828',tformat='YYYY_MMDD')
    stop_at = time_double('0000_0000',tformat='YYYY_MMDD')
    
    project = dictionary()
    project.name = 'azim_prop'
    project.var = project.name+'_project_info'
    project.time_range = time_double(['2012-10-01','2017-10-01'])   ; omni dst only updated to 2017-11
    project.root_dir = join_path([shomedir(),'azim_prop'])
    project.data_dir = join_path([project.root_dir,'data'])
    project.plot_dir = join_path([project.root_dir,'plot'])
    project.file = join_path([project.data_dir,project.name+'_project_info.tplot'])
    project.log_file = join_path([project.data_dir,project.name+'_gen_b_tilt_summary_plot.log'])
    if file_test(project.log_file) eq 0 then stouch, project.log_file

    project.probes = ['rbspa','rbspb','g13','g14','g15','tha','thd','the']
    project.short_probes = ['rba','rbb','g13','g14','g15','tha','thd','the']
    project.probe_colors = sgcolor(['red','tomato','orange','yellow','lime','cyan','deep_sky_blue','purple'])
    project.probe_psyms = replicate(6,n_elements(project.probes))
    project.movie_fps = 10d         ; frame per sec.
    project.movie_size = 800        ; pixel.
    project.movie_time_step = 600.  ; 10 min.
    
    model = 't89'
    orbit_data_rate = 60.       ; sec.
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

    constant = dictionary()
    constant.rgb = sgcolor(['red','green','blue'])
    constant.xyz = ['x','y','z']
    constant.deg = 180d/!dpi
    constant.rad = !dpi/180d
    constant.re = 6378d
    constant.secofday = 86400d
    console = -1
    tab = '    '
    project.constant = constant

    tplot_options, 'labflag', -1
    tplot_options, 'ynozero'
    tplot_options, 'ystyle', 1
    tplot_options, 'yticklen', -0.005
    tplot_options, 'xticklen', -0.03

    foreach dir, [project.root_dir,project.data_dir,project.plot_dir] do $
        if file_test(dir,/directory) eq 0 then file_mkdir, dir


;---Get the list of storms.
    storm_list = join_path([project.data_dir,'storm_list.txt'])
    if file_test(storm_list) eq 0 then begin
        project.dst_threshold = -50d    ; nT.
        dst_var = 'dst'
        
        omni = sread_omni(project.time_range)
        store_data, dst_var, sfmepoch(omni.epoch,'unix'), omni.sym_h, $
            limits={ytitle:'(nT)',labels:'Dst'}

        storm_times = []

        get_data, dst_var, times, dst
        idx = where(dst le project.dst_threshold, cnt)
        for i=0, cnt-1 do begin
            ut0 = times[idx[i]]
            ut1 = times[idx[i]]
            while i lt cnt-1 do begin
                if (times[idx[i+1]]-times[idx[i]]) lt constant.secofday then begin
                    ut1 = times[idx[i+1]]
                    i += 1
                endif else break
            endwhile
            time_range = [ut0,ut1]+[-1,1]*0.5*constant.secofday
            ;time_range = time_range-(time_range mod constant.secofday)
            storm_times = [storm_times, time_range]
        endfor
        nstorm = n_elements(storm_times)/2
        storm_times = reform(storm_times, 2,nstorm)

        openw, lun, storm_list, /get_lun
        printf, lun, time_string(storm_times)
        free_lun, lun
    endif else begin
        nstorm = file_lines(storm_list)
        storm_times = dblarr(2,nstorm)
        lines = strarr(nstorm)
        openr, lun, storm_list, /get_lun
        readf, lun, lines
        free_lun, lun
        for ii=0, nstorm-1 do storm_times[*,ii] = time_double(strsplit(lines[ii],' ',/extract))
    endelse
    ; Convert to a list.
    storms = list(length=nstorm)
    for ii=0, nstorm-1 do storms[ii] = storm_times[*,ii]
    project.storm_list = storms
    ; Save the above info to disk.
    store_data, project.var, 0, project
    tplot_save, project.var, filename=project.file



;---For each storm,
;   1. generate the movie of spacecraft position in the SM X-Y plane,
;   2. a line plot of B tilt in the SM coordinate.
    foreach storm_time_range, storms do begin
    ;---Settings.
        event_id = time_string(mean(storm_time_range), tformat='YYYY_MMDD')
        if time_double(event_id,tformat='YYYY_MMDD') le start_after then continue
        ;if time_double(event_id,tformat='YYYY_MMDD') ne stop_at then continue



    ;---Load Dst and AE.
        omni = sread_omni(storm_time_range)
        times = sfmepoch(omni.epoch,'unix')
        store_data, 'dst', times, omni.sym_h, $
            limits={ytitle:'(nT)',labels:'Dst'}
        store_data, 'ae', times, omni.ae_index, $
            limits={ytitle:'(nT)',labels:'AE'}

    ;---Clear memory.
        probes = project.probes
        nprobe = n_elements(probes)
        vars = ['r_gsm','b_gsm','r_sm','b_tilt','b_mag','bmod_tilt','bmod_mag','db_tilt','db_mag']
        foreach var, vars do store_data, probes+'_'+var, /delete


    ;---Load s/c location, and B field.
        foreach probe, probes do begin
            mission = resolve_probe(probe)
            call_procedure, mission.routine_name+'_read_bfield', storm_time_range, probe=mission.probe, errmsg=errmsg
            call_procedure, mission.routine_name+'_read_orbit', storm_time_range, probe=mission.probe
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
        lprmsg, event_id+tab+strjoin(time_string(storm_time_range),tab)+tab+strjoin(probes,','), project.log_file
        

    ;---Caclulate xxx_r_sm, xxx_b_tilt.
        foreach probe, probes do begin
            mission = resolve_probe(probe)

            times = make_bins(storm_time_range,orbit_data_rate)
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
            store_data, probe+'_bmod_gsm_'+model, times, bmodgsm


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

            bmodgsm = get_var_data(probe+'_bmod_gsm_'+model, at=times)
            bmodsm = cotran(bmodgsm, times, 'gsm2sm')
            bmod_tilt = atan(bmodsm[*,2],sqrt(bmodsm[*,0]^2+bmodsm[*,1]^2))*constant.deg
            store_data, probe+'_bmod_tilt', times, bmod_tilt, limits={$
                ytitle: '(deg)', labels:strupcase(mission.short_name+mission.probe)}
            bmod_mag = snorm(bmodgsm)
            store_data, probe+'_bmod_mag', times, bmod_mag, limits={$
                ytitle: '(nT)', labels:strupcase(mission.short_name+mission.probe)}

            dis = snorm(get_var_data(probe+'_r_gsm', at=times))
            db_tilt = b_tilt-bmod_tilt
            index = where(dis ge 4)
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
        plot_vars = ['b_tilt','db_tilt','db_mag']
        plot_titles = ['SM B tilt angle','SM B-B_mod tilt angle','|B|-|B_mod|']
        foreach var, plot_vars, ii do begin
            if keyword_set(test) then begin
                file = 0
            endif else begin
                file = project.plot_dir+'/'+event_id+'/'+event_id+'_'+var+'.pdf'
            endelse
            temp_dir = file_dirname(file)
            if file_test(temp_dir,/directory) eq 0 then file_mkdir, temp_dir

            fig_ysize = (2+nprobe)*1
            fig_xsize = total(storm_time_range*[-1,1])/constant.secofday*8
            sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch
            vars = ['dst','ae',probes+'_'+var]
            nvar = n_elements(vars)
            poss = sgcalcpos(nvar, margins=margins, xchsz=xchsz,ychsz=ychsz)
            tplot, vars, trange=storm_time_range, /novtitle, position=poss
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
        vidstream = ovid.AddVideoStream(project.movie_size, project.movie_size, project.movie_fps)
        temp_file = project.root_dir+'/tmp.png'

        times = make_bins(storm_time_range, project.movie_time_step)

        foreach time, times, i do begin
            printf, console, time_string(time)
            sgopen, temp_file, xsize=project.movie_size, ysize=project.movie_size
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
            foreach probe, project.short_probes, jj do begin
                tx = tpos[0]+xchsz*(25+5*jj)
                xyouts, tx,ty,/normal, strupcase(project.short_probes[jj]), color=project.probe_colors[jj]
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
                plots, rsm[0],rsm[1],/data, color=project.probe_colors[jj], psym=project.probe_psyms[jj]
            endforeach

            sgclose
            read_png, temp_file, temp_image
            time = ovid.put(vidstream, temp_image)
        endforeach
        ovid.Cleanup
        file_delete, temp_file
        
        ;if time_double(event_id,tformat='YYYY_MMDD') eq stop_at then stop
    endforeach

;---Update info.
    lines = read_all_lines(project.log_file)
    event_info = hash()
    foreach line, lines do begin
        tinfo = strsplit(line,' ',/extract)
        event_id = tinfo[0]
        time_range = time_double(tinfo[0:1])
        probes = strsplit(tinfo[2],',',/extract)
        event_info[event_id] = dictionary($
            'time_range', time_range, $
            'probes', probes)
    endforeach
    project.event_info = event_info
    store_data, project.var, 0, project
    tplot_save, project.var, filename=project.file

end
