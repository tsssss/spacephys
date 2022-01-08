;+
; Plot the result of velocity from 2d time lag.
;-

pro azim_prop_plot_2dv_vs_r, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    test = 0
    
    
    ; Scheme 7, 6 colors. https://www.colorcombos.com/color-scheme-4291.html. Note: rgb is reversed in IDL.
    colors = ['4567F1'x,'5DC6FF'x,'A4C87B'x,'D9C34C'x,'8D6493'x,'404040'x]
    colors = colors[[4,3,2,1,0,5]]  ; reversed rainbow.

;---Plot settings.
    margins = [8,4,5,2]
    margins = [6,4,3,2]
    fig_xsize = 3
    fig_ysize = 3
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

    min_ncombo = 1
    
    file = join_path([project.plot_dir,'fig_2dv_vs_r.pdf'])
    magnify = 1
    if keyword_set(test) then begin
        file = test
        magnify = 2
    endif
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify
    tpos = sgcalcpos(1, margins=margins, xchsz=xchsz,ychsz=ychsz)
    yrange = [0,80]
    yticklen = -0.01
    ytitle = '|V!D2D!N| (km/s)'
    ytickv = [0,30,60]
    yminor = 5
    yticks = n_elements(ytickv)-1
    xrange = [0,10]
    xticklen = -0.01
    xtitle = 'R!Dxy!N (Re)'
    plot, xrange, yrange, /nodata, $
        xstyle=1, xrange=xrange, xticklen=xticklen, xtitle=xtitle, $
        ystyle=1, yrange=yrange, yticklen=yticklen, ytitle=ytitle, yticks=yticks, ytickv=ytickv, yminor=yminor, $
        position=tpos

    sorted_events = (events.keys()).sort()
    color_index = 0
    foreach event_id, sorted_events do begin
        data_file = join_path([project.data_dir,event_id+'_basic_data.tplot'])
        if file_test(data_file) eq 0 then azim_prop_gen_survey_plot, /save_data
        if file_test(data_file) eq 0 then message, 'Something is wrong when loading the data ...'
        events[event_id].file = data_file
        tplot_restore, file = data_file

    ;---Add data for all combos.
        v2d_angle = project.v2d_angle
        v2d_dt = project.v2d_dt
        xs = list()
        ys = list()
        combos = project.events[event_id].probe_combos
        foreach combo, combos do begin
            key = strjoin(combo[sort(combo)],'_')
            combo_info = (project.events[event_id])[key]
            if ~combo_info.good_combo then continue

            v_mag = combo_info.v_mag
            angle_from_azimuth = combo_info.angle_from_azimuth
            dis_center = snorm(combo_info.rsm_center)
            v_azim = v_mag*cos(angle_from_azimuth*rad)
            omega_2d = v_mag/dis_center/re*deg*60
            
            tx = dis_center
            ty = v_mag
            xs.add, tx
            ys.add, ty
        endforeach
        print, event_id
        if n_elements(xs) le min_ncombo then continue
        xs = xs.toarray()
        ys = ys.toarray()
        
        color = colors[color_index]
        color_index += 1
        plots, xs, ys, psym=1, color=color, symsize=0.5
        fit_res = linfit(xs,ys)
        omega_fit = fit_res[1]/re*deg*60
        
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*(color_index)*label_size
        xyouts, tx,ty,/normal, strmid(event_id,0,9), color=color, charsize=label_size
    endforeach
    
;---Close the canvas.
    if keyword_set(test) then stop
    sgclose
end
