;+
; Plot all the data for all events, each panel has all spacecraft data shifted by the proper time lags.
; Not time-lagged.
; Color-coded with discrete colors.
;-

pro azim_prop_plot_all_data4, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    test = 0

;---Colors.
    ; Scheme 1.
    colors = list()
    colors.add, [38,115,77]
    colors.add, [45,134,134]
    colors.add, [51,102,153]
    colors.add, [57,57,172]
    colors.add, [127,64,191]
    colors.add, [134,45,134]
    
    ; Scheme 2.
    colors = list()
    colors.add, [85,64,191]
    colors.add, [64,149,191]
    colors.add, [64,191,128]
    colors.add, [191,128,64]
    colors.add, [179,64,64]
    colors.add, [179,77,110]
    
    ; Scheme 4.
    colors = list()
    colors.add, [51,51,204]
    colors.add, [51,178,204]
    colors.add, [170,191,64]
    colors.add, [242,166,13]
    colors.add, [217,68,38]    
    colors.add, [153,51,153]
    
;    ; Scheme 5, blue-green.
;    colors = list()
;    colors.add, [25,  0,153]
;    colors.add, [18, 65,161]
;    colors.add, [31,126,173]
;    colors.add, [46,184,184]
;    colors.add, [57,198,151]
;    colors.add, [83,198,121]
;    
;    ; Scheme 6, purple-blue-green.
;    colors = list()
;    colors.add, [102,  0,153]
;    colors.add, [ 42, 18,161]
;    colors.add, [ 31, 78,173]
;    colors.add, [ 46,161,184]
;    colors.add, [ 64,191,149]
;    colors.add, [ 94,186,110]


    ; switch order.    
    foreach rgb, colors, ii do colors[ii] = 256L*(256L*rgb[2]+rgb[1])+rgb[0] ; bbggrr.

    ; Scheme 3.
    ;colors = reverse(sgcolor(['tomato','orange','lime_green','deep_sky_blue','blue','purple']))

    ; Scheme 7, 6 colors. https://www.colorcombos.com/color-scheme-4291.html. Note: rgb is reversed in IDL.
    colors = ['4567F1'x,'5DC6FF'x,'A4C87B'x,'D9C34C'x,'8D6493'x,'404040'x]
    ;colors = colors[[5,3,4,2,0,1]] ; dark-light-dark-light.
    colors = colors[[5,4,3,2,1,0]]  ; reversed rainbow.
    ;colors = colors[[4,3,2,1,0,5]]  ; reversed rainbow, black last.
    
    ; Scheme 8, 6 colors, the most saturated colors.
    colors = list()
    colors.add, [  0,  0,255]   ; blue.
    colors.add, [  0,230,230]   ; cyan.
    colors.add, [230, 25, 25]   ; red.
    colors.add, [230,191,  0]   ; yellow.
    colors.add, [  0,204,  0]   ; green.
    colors.add, [204,  0,204]   ; purple.
    foreach rgb, colors, ii do colors[ii] = 256L*(256L*rgb[2]+rgb[1])+rgb[0] ; bbggrr.

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
    data_thick = 1
    ; For plotting sc.
    sc_psym = 8
    sc_size = 0.3
    sc_thick = 3
    ; For plotting arrows.
    std_vmag = 40.
    std_length = 2
    coef2length = std_length/std_vmag
    arrow_solid = 1
    arrow_hsize = keyword_set(test)? 4:40
    arrow_thick = 0.5
    arrow_color = sgcolor('black')

    
    tmp = findgen(21)/20*2*!dpi
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs,tys, thick=sc_thick

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
    file = join_path([project.plot_dir,'fig_all_'+data_type+'4.pdf'])
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
    tx =  (tpos[0]+tpos[2])*0.5
    ty = tpos[3]+ychsz*1
    xyouts, tx,ty,/normal,alignment=0.5, 'Tilt angle '+str_delta+str_theta+' time lag removed'
;    tpos = poss[*,1,0]
;    tx = 0.5*(tpos[0]+tpos[2])-xchsz*1.5
;    ty = tpos[3]+ychsz*2
;    xyouts, tx,ty,/normal,alignment=0.5, 'SM X-Y plane'


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
            color = colors[ii]
            var = probe+'_'+data_type
            get_data, probe+'_'+data_type, times, data
            ;data = deriv(data)/5*10
            ;data = smooth(data,10)
            time_lag = event[probe].time_lag
            times = times+time_lag
            data = data-offsets[ii]
            oplot, times, data, color=color, thick=data_thick

            sc_short = strupcase(project[probe].short_name)
            xyouts, tx,ty,/normal, sc_short, color=color
            tx += strlen(sc_short)*xchsz
        endforeach

        txs = [xchsz*2,1-xchsz*4]
        tys = tpos[1]-ychsz*1.2+[0,0]
        plots, txs,tys,/normal, thick=0.2;, linestyle=1

        ; Add duration.
        ;tx = tpos[0]+xchsz*0.5
        tx = lmargin*xchsz
        ty = tpos[3]-ychsz*0.1
        duration = total(event_time_range*[-1,1])/3600.
        xyouts, tx,ty,/normal,alignment=0, 'T = '+string(duration,format='(F3.1)')+' hr'


        ; Add ramp duration and scale.
        tilt_time_range = event.tilt_time_range
        duration = event.tilt_duration/60    ; in min.
        scale = event.tilt_scale             ; in Re.
        txs = tilt_time_range
        tmp = convert_coord(txs,[0,0], /data, /to_normal)
        txs = reform(tmp[0,*])
        tys = tpos[1]+ychsz*0.3+[0,0]
        scale_color = sgcolor('black')
        plots, txs,tys,/normal, color=scale_color
        foreach tx, txs do plots, tx+[0,0],tys[0]+[-1,1]*ychsz*0.2,/normal, color=scale_color
        msg = string(duration,format='(I0)')+' min, '+string(scale,format='(F3.1)')+' Re'
        tx = txs[1]+xchsz*1
        ty = tys[0]-ychsz*0.1
        xyouts, tx,ty,/normal, alignment=0, msg, charsize=label_size, color=scale_color


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

        ; Add sc.
        foreach probe, sorted_probes, ii do begin
            ref_mlt = event[probe].ref_mlt
            color = colors[ii]

            ref_rsm = event[probe].ref_rsm
            tx = ref_rsm[0]
            ty = ref_rsm[1]
            plots, tx,ty,/data, psym=sc_psym, symsize=sc_size, color=color
        endforeach
        
        ; Add arrows, from azim_prop_plot_2d_timing.
        v2d_angle = project.v2d_angle
        v2d_dt = project.v2d_dt
        omega_2ds = []
        combos = project.events[event_id].probe_combos
        combo_index = []
        foreach combo, combos, ii do begin
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
            
            combo_index = [combo_index,ii]
        endforeach
        
        if n_elements(combo_index) lt 2 then continue
        
        combos = combos[combo_index]
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

            x0 = rsm_center[0]
            y0 = rsm_center[1]
            scale = coef2length*v_mag
            x1 = x0+v_hat[0]*scale
            y1 = y0+v_hat[1]*scale
            arrow, x0,y0,x1,y1,/data, hsize=arrow_hsize, solid=arrow_solid, thick=arrow_thick, color=arrow_color
        endforeach

        ;if keyword_set(test) then stop
    endforeach

    if keyword_set(test) then stop
    sgclose

end
