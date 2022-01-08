;+
; Use 3 spacecraft combinations to calculate the 2d velocity. Need to set reference time of the first spacecraft to get the MLT; and the spacecraft combinations.
;-


pro fig_2014_0828_triangle_collage, project

    project = azim_prop_load_project()
    stop
    events = project.events

    dim_index = [0,1]   ; 2-d plane, the SM x-y plane.
    deg = project.constant.deg
    rad = project.constant.rad
    re = project.constant.re
    label_size = 0.8
    large_size = 1.5
    test = 0

;---Settings to ensure v_2d calculation is reliable.
    v2d_angle = project.v2d_angle   ; triangle angles in [15,165] deg.
    v2d_dt = project.v2d_dt         ; in sec, min time difference.
    black = sgcolor('black')
    bad_color = sgcolor('red')
    good_color = sgcolor('blue')

    event_id = '2014_0828_10'
    lprmsg, 'Processing '+event_id+' ...'
    data_file = join_path([project.data_dir,event_id+'_basic_data.tplot'])
    if file_test(data_file) eq 0 then azim_prop_gen_survey_plot, /save_data
    if file_test(data_file) eq 0 then message, 'Something is wrong when loading the data ...'
    tplot_restore, file = data_file

    event_info = events[event_id]
    probe_combos = event_info.probe_combos
    nprobe_combo = n_elements(probe_combos)

;---Plot settings.
    label_letters = letters(nprobe_combo)
    nxpanel = 5
    nypanel = ceil(nprobe_combo/nxpanel)
    panel_xrange = [5,-15]
    panel_yrange = [5,-15]
    panel_xsize = 2
    panel_ysize = panel_xsize*total(panel_yrange*[-1,1])/total(panel_xrange*[-1,1])
    sgopen, 0, xsize=panel_xsize, ysize=panel_ysize, /inch
    tmp = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz, position=[0,0,1,1])
    sgclose, /wdelete
    margins = [7.,5,2,1]
    xpad = 1*ychsz/xchsz
    ypad = 1
    fig_xsize = panel_xsize*nxpanel+(xpad*(nxpanel-1)+(margins[0]+margins[2]))*xchsz*panel_xsize
    fig_ysize = panel_ysize*nypanel+(ypad*(nypanel-1)+(margins[1]+margins[3]))*ychsz*panel_ysize

    if keyword_set(test) then begin
        file = test
        magnify = 1
        hsize = 2
    endif else begin
        file = join_path([project.plot_dir,'triangle_collage','fig_'+event_id+'_triangle_collage.pdf'])
        magnify = 1
        hsize = 160
    endelse
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify, /inch
    panel_poss = sgcalcpos(nypanel,nxpanel, xpad=xpad, ypad=ypad, margins=margins, xchsz=xchsz, ychsz=ychsz)


;---Treat each combo.
    foreach probes, probe_combos, panel_id do begin
        probe_combo = strjoin(probes[sort(probes)],'_')
        nprobe = n_elements(probes)
        ref_times = dblarr(nprobe)
        rsms = fltarr(3,nprobe)
        foreach probe, probes, ii do begin
            ref_times[ii] = event_info[probe].ref_time
            rsms[*,ii] = event_info[probe].ref_rsm
        endforeach
        index = sort(ref_times)
        ref_times = ref_times[index]
        rsms = rsms[*,index]
        probes = probes[index]  ; need to sort by ref_time, cannt read from event_info.

        combo_info = event_info[probe_combo]
        v_mag = combo_info.v_mag
        v_hat = combo_info.v_hat
        rsm_center = combo_info.rsm_center

        index = array_indices([nxpanel,nypanel], panel_id, /dimensions)
        tpos = panel_poss[*,index[0],index[1]]
        xtitle = (index[1] eq nypanel-1)? 'SM X (Re)': ''
        ytitle = (index[0] eq 0)? 'SM Y (Re)': ''
        xtickformat = (index[1] eq nypanel-1)? '(I0)': '(A1)'
        ytickformat = (index[0] eq 0)? '(I0)': '(A1)'
        plot, panel_xrange, panel_yrange, xstyle=1, ystyle=1, position=tpos, $
            xrange=panel_xrange, yrange=panel_yrange, /nodata, /noerase, /isotropic, $
            xtitle=xtitle, ytitle=ytitle, xticklen=-0.01, yticklen=-0.01, xtickformat=xtickformat, ytickformat=ytickformat

        ; Add figure label.
        tx = tpos[2]-xchsz*1
        ty = tpos[3]-ychsz*1
        label = +'('+sgnum2str(index[1]+1)+'-'+sgnum2str(index[0]+1)+')'
        xyouts, tx,ty,/normal,alignment=1, label

        ; Add geometry and timing criteria, add flag.
        good_combo = combo_info.good_combo
        msg = good_combo? 'Good': 'Bad'
        triad_color = good_combo? good_color: bad_color
        msg = msg+' triad'
        tx = tpos[2]-xchsz*1
        ty = tpos[1]+ychsz*0.5
        xyouts, tx,ty,/normal,alignment=1, msg, color=triad_color, charsize=large_size

        ; Add earth.
        tmp = 50
        tmp = findgen(tmp)/(tmp-1)*2*!dpi
        xs = cos(tmp)
        ys = sin(tmp)
        polyfill, xs<0, ys, /line_fill, orientation=45, spacing=ychsz*5
        plots, xs, ys
        foreach r, [5,10] do oplot, xs*r, ys*r, linestyle=1

        ; Add spacecraft and labels.
        for ii=0, nprobe-1 do begin
            tx = rsms[dim_index[0],ii]
            ty = rsms[dim_index[1],ii]
            plots, tx,ty, psym=6, symsize=label_size*0.5
            short_name = project[probes[ii]].short_name
            xyouts, tx,ty,/data, '  '+strupcase(short_name), charsize=label_size
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

        tx = tpos[0]+xchsz*1
        triangles = combo_info.triangle_angles
        triangles = triangles[sort(triangles)]
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty, /normal, 'Triad (deg):', charsize=label_size
        foreach triangle, triangles, jj do begin
            msg = strtrim(string(triangle,format='(F5.1)'),2)
            if jj ne 2 then msg += ','
            bad_flag = triangle le v2d_angle or triangle ge (180-v2d_angle)
            color = bad_flag? bad_color: black
            ttx = tx+(6+3*jj)*xchsz
            xyouts, ttx,ty,/normal, msg, charsize=label_size, color=color
        endforeach


        ; Add time info.
        ty = tpos[3]-ychsz*2
        time_diffs = abs([ref_times[1]-ref_times[0],ref_times[2]-ref_times[1],ref_times[0]-ref_times[2]])
        time_diffs = time_diffs[sort(time_diffs)]
        combo_info.time_diffs = time_diffs
        xyouts, tx,ty, /normal, 'T.diff (sec):', charsize=label_size
        foreach time_diff, time_diffs, jj do begin
            msg = strtrim(string(round(time_diff),format='(I0)'),2)
            if jj ne 2 then msg += ','
            bad_flag = time_diff le v2d_dt
            color = bad_flag? bad_color: black
            ttx = tx+(6+3*jj)*xchsz
            xyouts, ttx,ty,/normal, msg, charsize=label_size, color=color
        endforeach


        ; Add velocity info.
        ty = tpos[3]-ychsz*3
        xyouts, tx,ty, /normal, charsize=label_size, '|V| = '+string(v_mag,format='(F5.1)')+' km/s'

        if keyword_set(test) then stop
    endforeach

    if keyword_set(test) then stop
    sgclose
end
