;+
; Use 3 spacecraft combinations to calculate the 2d velocity. Need to set reference time of the first spacecraft to get the MLT; and the spacecraft combinations.
;-


pro azim_prop_calc_2d_vel, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    if ~project.haskey('done_calc_time_lag') then azim_prop_calc_time_lag, project
    events = project.events

    dim_index = [0,1]   ; 2-d plane, the SM x-y plane.
    deg = project.constant.deg
    rad = project.constant.rad
    re = project.constant.re
    label_size = 0.8
    test = 0

;---Settings to ensure v_2d calculation is reliable.
    v2d_angle = 15.     ; triangle angles in [15,165] deg.
    v2d_dt = 120.       ; in sec, min time difference.
    project.v2d_angle = v2d_angle
    project.v2d_dt = v2d_dt


    foreach event_id, events.keys() do begin
        lprmsg, 'Processing '+event_id+' ...'

        event_info = events[event_id]
        probe_combos = choose_from(event_info.probes, 3)
        project.events[event_id].probe_combos = probe_combos

        data_file = join_path([project.data_dir,event_id+'_basic_data.tplot'])
        if file_test(data_file) eq 0 then azim_prop_gen_survey_plot, /save_data
        if file_test(data_file) eq 0 then message, 'Something is wrong when loading the data ...'
        project.events[event_id].file = data_file
        tplot_restore, file = data_file


    ;---Treat each combo.
        foreach probes, probe_combos do begin
            probe_combo = strjoin(probes[sort(probes)],'_')
            nprobe = n_elements(probes)
            times = dblarr(nprobe)
            rsms = fltarr(3,nprobe)
            mlts = fltarr(nprobe)   ; the MLT's of the probes.
            rxys = fltarr(nprobe)   ; the Rxy's of the probes.
            foreach probe, probes, ii do begin
                times[ii] = event_info[probe].ref_time
                rsms[*,ii] = event_info[probe].ref_rsm
                mlts[ii] = event_info[probe].ref_mlt
                rxys[ii] = snorm(rsms[0:1,ii])
            endforeach
            rsm_center = total(rsms[dim_index,*],2)/nprobe

            index = sort(times)
            times = times[index]
            rsms = rsms[*,index]
            mlts = mlts[index]
            rxys = rxys[index]
            angles = atan(rsms[2,*],rxys)*deg
            probes = probes[index]

            tt = times[1:nprobe-1]-times[0]
            rr = fltarr(3,nprobe-1) & for ii=0, nprobe-2 do rr[*,ii] = rsms[*,ii+1]-rsms[*,0]
            rr = rr[dim_index,*]
            vv = la_linear_equation(rr,tt)
            v_hat = sunitvec(vv)
            rr_normal = dblarr(nprobe-1)
            for ii=0, nprobe-2 do rr_normal[ii] = sdot(rr[*,ii],v_hat)
            fit_result = linfit(tt, rr_normal)
            v_mag = fit_result[1]*re
            angle_from_azimuth = sang(rsm_center, v_hat, /degree)-90
            if angle_from_azimuth gt 180 then angle_from_azimuth -= 360
            if angle_from_azimuth lt -180 then angle_from_azimuth += 360
            triangle_angles = [$
                sang(rr[*,0]-0,rr[*,1]-0,/degree), $
                sang(0-rr[*,0],rr[*,1]-rr[*,0],/degree), $
                sang(0-rr[*,1],rr[*,0]-rr[*,1],/degree)]

            if (project.events[event_id]).haskey('probe_combo') then (project.events[event_id]).remove, 'probe_combo'
            combo_info = dictionary($
                'times', times, $
                'mlts', mlts, $
                'rxys', rxys, $
                'rsms', rsms, $
                'angles', angles, $
                'mlt_range', minmax(mlts), $
                'rxy_range', minmax(rxys), $
                'z_range', minmax(rsms[1,*]), $
                'angle_range', minmax(angles), $
                'v_mag', v_mag, $
                'v_hat', v_hat, $
                'rsm_center', rsm_center, $
                'angle_from_azimuth', angle_from_azimuth)
            if keyword_set(test) then begin
                file = test
                magnify = 2
                hsize = 20
            endif else begin
                file = join_path([project.plot_dir,'v_2d','fig_'+event_id+'_'+strjoin(probes,'_')+'_v2d.pdf'])
                magnify = 1
                hsize = 160
            endelse

            lprmsg, probes
            lprmsg, 'v_mag (km/s): '+string(v_mag)
            lprmsg, 'v_hat (x,y) : '+strjoin(string(v_hat),',')

            sgopen, file, xsize=4, ysize=6, magnify=magnify
            tpos = sgcalcpos(1, position=[0.2,0.2,0.9,0.9], xchsz=xchsz, ychsz=ychsz)
            plot, [5,-15], [15,-15], xstyle=1, ystyle=1, position=tpos, $
                xrange=[5,-15], yrange=[15,-15], /nodata, /noerase, /isotropic, $
                xtitle='SM X (Re)', ytitle='SM Y (Re)', xticklen=-0.01, yticklen=-0.02

            ; Add earth.
            tmp = 50
            tmp = findgen(tmp)/(tmp-1)*2*!dpi
            xs = cos(tmp)
            ys = sin(tmp)
            polyfill, xs<0, ys, /line_fill, orientation=45
            plots, xs, ys
            foreach r, [5,10] do oplot, xs*r, ys*r, linestyle=1

            ; Add spacecraft and labels.
            for ii=0, nprobe-1 do begin
                tx = rsms[dim_index[0],ii]
                ty = rsms[dim_index[1],ii]
                plots, tx,ty, psym=6, symsize=label_size*0.5
                xyouts, tx,ty,/data, '  '+strupcase(probes[ii]), charsize=label_size
            endfor

            ; Add triangle.
            tx = rsms[dim_index[0],*]
            ty = rsms[dim_index[1],*]
            for ii=0, nprobe-1 do begin
                plots, (shift(tx,ii))[0:1], (shift(ty,ii))[0:1], linestyle=1
            endfor

            ; Add arrow.
            x0 = mean(rsms[dim_index[0],*])
            y0 = mean(rsms[dim_index[1],*])
            scale = 2d/20*v_mag
            x1 = x0+v_hat[0]*scale
            y1 = y0+v_hat[1]*scale
            arrow, x0,y0,x1,y1,/data, /solid, hsize=hsize

            ; Add triangle info.
            tx = tpos[0]+xchsz*1
            max_len = 5
            for ii=0, nprobe-1 do begin
                ty = tpos[3]-ychsz*(ii+2)
                xyouts, tx,ty,/normal, strupcase(probes[ii]), charsize=label_size
                xyouts, tx+xchsz*max_len,ty,/normal, charsize=label_size, ' SM R (Re): ('+strjoin(strtrim(string(rsms[*,ii],format='(F5.1)'),2),', ')+')'
            endfor

            triangles = triangle_angles
            triangles = triangles[sort(triangles)]
            combo_info.triangle_angles = triangles
            ty = tpos[3]-ychsz*5
            xyouts, tx,ty, /normal, 'Triangle (deg): '+strjoin(strtrim(string(triangles,format='(F5.1)'),2),', '), charsize=label_size

            ; Add time info.
            ty = tpos[3]-ychsz*6
            time_diffs = abs([times[1]-times[0],times[2]-times[1],times[0]-times[2]])
            time_diffs = time_diffs[sort(time_diffs)]
            combo_info.time_diffs = time_diffs
            xyouts, tx,ty, /normal, 'Time diff (sec): '+strjoin(strtrim(string(round(time_diffs),format='(I0)'),2),', '), charsize=label_size

            ; Add MLT info.
            ty = tpos[3]-ychsz*7
            mlt_diffs = abs([mlts[1]-mlts[0],mlts[2]-mlts[1],mlts[0]-mlts[2]])
            xyouts, tx,ty, /normal, 'MLT diff (hr): '+strjoin(strtrim(string(mlt_diffs,format='(F4.1)'),2),', '), charsize=label_size

            ; Add velocity info.
            ty = tpos[3]-ychsz*1
            xyouts, tx,ty, /normal, charsize=label_size, '|V| = '+string(v_mag,format='(F5.1)')+' km/s'
            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            xyouts, tx,ty, /normal, event_id+'. V_2D analysis among ['+strjoin(strupcase(probes),', ')+']', charsize=label_size
            if keyword_set(test) then stop
            sgclose

            ; Add omega_2d, in deg/min.
            omega_2d = v_mag/snorm(rsm_center)/re*deg*cos(angle_from_azimuth*rad)*60
            combo_info.omega_2d = omega_2d

            ; Apply geometry and timing criteria, add flag.
            good_combo = 1
            if max(triangle_angles) gt 180-project.v2d_angle then good_combo = 0
            if min(triangle_angles) lt project.v2d_angle then good_combo = 0
            if min(time_diffs) lt project.v2d_dt then good_combo = 0
            combo_info.good_combo = good_combo

            ; Add change to project.
            (project.events[event_id])[probe_combo] = combo_info
        endforeach

    ;---Get the overall omega_2d.
        omega_2ds = list()
        foreach probes, probe_combos do begin
            probe_combo = strjoin(probes[sort(probes)],'_')
            combo_info = (project.events[event_id])[probe_combo]
            if ~combo_info.good_combo then continue
            omega_2ds.add, combo_info.omega_2d
        endforeach
        omega_2ds = omega_2ds.toarray()
        project.events[event_id].omega_2ds = omega_2ds
        nomega_2d = n_elements(omega_2ds)
        if nomega_2d eq 0 then continue
        project.events[event_id].omega_2d = mean(omega_2ds)
        if nomega_2d eq 1 then continue
        project.events[event_id].domega_2d = stddev(omega_2ds)
    endforeach

    project.done_calc_2d_vel = 1
    azim_prop_update_project, project

end
