;+
; Plot the result of velocity from 2d time lag.
;-

pro azim_prop_plot_2d_timing, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    test = 0

;---Plot settings.
    margins = [8,4,5,2]
    margins = [6,4,3,2]
    fig_xsize = 3
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
    re = project.constant.re
    rad = project.constant.rad
    r_kms = 6.6
    coef2kms = r_kms*re*rad/60
    fig_labels = ['a','b','c','d','e','f','g','h','i','j','k']+'.'
    label_xpos = 5
    std_vmag = 40.
    std_length = 2.
    coef2length = std_length/std_vmag
    deg = project.constant.deg
    str_pm = '!9'+string(177b)+'!X'


    foreach event_id, events.keys() do begin
        data_file = join_path([project.data_dir,event_id+'_basic_data.tplot'])
        if file_test(data_file) eq 0 then azim_prop_gen_survey_plot, /save_data
        if file_test(data_file) eq 0 then message, 'Something is wrong when loading the data ...'
        events[event_id].file = data_file
        tplot_restore, file = data_file

        probes = events[event_id].probes


    ;---Plot settings.
        event_time_range = events[event_id].time_range
        xrange = [-1,1]
        yrange = [-1,1]
        foreach probe, probes do begin
            rgsm = get_var_data(probe+'_r_gsm', at=event_time_range)
            rsm = cotran(rgsm, event_time_range, 'gsm2sm')
            xrange = [xrange, rsm[*,0]]
            yrange = [yrange, rsm[*,1]]
        endforeach
        xrange = [floor(min(xrange)),ceil(max(xrange))]
        yrange = [floor(min(yrange)),ceil(max(yrange))]

        xstep = 5
        xtickv = make_bins(xrange, xstep)
        xticks = n_elements(xtickv)-1
        xminor = 5
        xtitle = 'SM X (Re)'
        xrange = reverse(minmax(xtickv))
        ystep = 5
        ytickv = make_bins(yrange, ystep)
        yticks = n_elements(ytickv)-1
        yminor = 5
        ytitle = 'SM Y (Re)'
        yrange = reverse(minmax(ytickv))
        aspect_ratio = total(abs(yrange*[-1,1]))/total(abs(xrange*[-1,1]))
        xticklen = yticklen/aspect_ratio

        if keyword_set(test) then begin
            file = test
            magnify = 2
            hsize = 10
        endif else begin
            file = join_path([project.plot_dir,'time_lag','fig_'+event_id+'_plot_2d_timing.pdf'])
            magnify = 1
            hsize = 120
        endelse


        ; Calculate ysize.
        sgopen, 0, xsize=fig_xsize, ysize=fig_xsize
        tpos = sgcalcpos(1, margins=margins, xchsz=xchsz,ychsz=ychsz)
        fig_ysize = (tpos[1]+1-tpos[3]+(tpos[2]-tpos[0])*aspect_ratio)*fig_xsize
        sgclose, /wdelete

    ;---Open the canvas.
        sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify
        tpos = sgcalcpos(1, margins=margins, xchsz=xchsz,ychsz=ychsz)

    ;---Create box, add earth and lines.
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtitle=xtitle, $
            ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            position=tpos, /nodata, /noerase, /isotropic

        ; Add earth.
        tmp = 50
        tmp = findgen(tmp)/(tmp-1)*2*!dpi
        xs = cos(tmp)
        ys = sin(tmp)
        polyfill, xs<0, ys, /line_fill, orientation=45
        plots, xs, ys
        foreach r, [5,10] do oplot, xs*r, ys*r, linestyle=1

        ; Add probes.
        probes = project.events[event_id].probes
        foreach probe, probes do begin
            ref_rsm = (project.events[event_id])[probe].ref_rsm
            label = project[probe].short_name
            color = project[probe].color
            tx = ref_rsm[0]
            ty = ref_rsm[1]
            plots, tx,ty,/data, psym=psym, symsize=label_size*0.5, color=color
            tmp = convert_coord([tx,ty], /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*full_ychsz
            xyouts, tx,ty,/normal, strupcase(label), alignment=0.5, charsize=label_size, color=color
        endforeach

    ;---Add legend.
        ty0 = tpos[3]-ychsz*label_size
        tmp = convert_coord([0,0],/data,/to_normal)
        x0 = tmp[0]
        tmp = convert_coord([std_length,0],/data,/to_normal)
        x1 = tmp[0]
        dx = abs(x1-x0)
        txs = tpos[2]-xchsz*8+[0,dx]
        tys = ty0+ychsz*(half_ychsz*0.5)+[0,0]
        plots, txs,tys,/normal
        foreach tx, txs, ii do plots, tx+[0,0],tys[ii]+[-1,1]*ychsz*0.1,/normal
        tx = max(txs)+xchsz*0.5
        ty = ty0
        xyouts, tx,ty,/normal, sgnum2str(std_vmag)+' km/s', charsize=label_size

        ;Add event id.
        ;tx = tpos[0]
        ;ty = tpos[3]+ychsz*line_skip
        ;xyouts, tx,ty,/normal, event_id, charsize=label_size


    ;---Add data for all combos.
        v2d_angle = project.v2d_angle
        v2d_dt = project.v2d_dt
        omega_2ds = []
        combos = project.events[event_id].probe_combos
        foreach combo, combos do begin
            key = strjoin(combo[sort(combo)],'_')
            combo_info = (project.events[event_id])[key]
            v_mag = combo_info.v_mag
            v_hat = combo_info.v_hat
            rsm_center = combo_info.rsm_center

            ; Filter by time difference.
            times = combo_info.times
            time_diffs = abs([times[0]-times[1],times[1]-times[2],times[2]-times[0]])
            if min(time_diffs) lt v2d_dt then continue

            ; Filter by geometry.
            triangle_angle = combo_info.triangle_angles
            if max(triangle_angle) gt 180-v2d_angle then continue
            if min(triangle_angle) lt v2d_angle then continue
            ;rsms = (combo_info.rsms)[0:1,*]
            ;diss = [snorm(rsms[*,0]-rsms[*,1]),snorm(rsms[*,1]-rsms[*,2]),snorm(rsms[*,2]-rsms[*,0])]
            ;if min(diss) lt small_dis then continue
            angle_from_azimuth = combo_info.angle_from_azimuth
            x0 = rsm_center[0]
            y0 = rsm_center[1]
            scale = coef2length*v_mag
            x1 = x0+v_hat[0]*scale
            y1 = y0+v_hat[1]*scale
            arrow, x0,y0,x1,y1,/data, hsize=hsize, /solid
            tmp = convert_coord([x0,y0], /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*1
            msg = strupcase(key)+$
                '!Cv_mag='+sgnum2str(v_mag,ndec=1)+$
                'km/s!Cazim='+sgnum2str(angle_from_azimuth,ndec=0)
            ;xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size*0.5
            lprmsg, msg
            print, 'MLT range=', combo_info.mlt_range
            print, 'Rxy range=', combo_info.rxy_range
            print, 'Z range=', combo_info.z_range
            ;stop

            ; angular speed in deg/min.
            omega_2d = v_mag/snorm(rsm_center)/re*deg*cos(angle_from_azimuth*rad)*60
            omega_2ds = [omega_2ds, omega_2d]
        endforeach

    ;---Add figure label.
        tx = tpos[0]
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,/normal, 'b. Permutations of 3-spacecraft timing'


    ;---Add the angular speed.
        if n_elements(omega_2ds) lt 1 then begin
            sgclose
            continue
        endif
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*label_size
        alignment = 0
        omega_2d = mean(omega_2ds)
        domega_2d = stddev(omega_2ds)
        msg = ['|'+str_omega+'!D2D!N| = '+sgnum2str(omega_2d,ndec=1)+str_pm+sgnum2str(domega_2d,ndec=1)+' deg/min', $
            sgnum2str(round(omega_2d*coef2kms))+str_pm+sgnum2str(round(domega_2d*coef2kms))+' km/s at '+sgnum2str(r_kms)+' Re']
        if n_elements(omega_2ds) eq 1 then $
            msg = ['|'+str_omega+'!D2D!N = '+sgnum2str(omega_2d,ndec=1)+' deg/min, ', $
                sgnum2str(round(omega_2d*coef2kms))+' km/s at '+sgnum2str(r_kms)+' Re']
        foreach tmp, msg, ii do $
            xyouts, tx,ty-ychsz*ii*label_size,/normal,alignment=alignment, charsize=label_size, tmp

    ;---Close the canvas.
        if keyword_set(test) then stop
        sgclose
    endforeach
end
