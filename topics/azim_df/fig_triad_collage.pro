;+
; Plot Triad collage to show good and bad triads.
;-

pro fig_triad_collage, event, project=project


test = 0

;---Check inputs.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    azim_df_load_basic_data, project=project


    deg = constant('deg')
    rad = constant('rad')
    re = constant('re')
    label_size = 0.8
    large_size = 1.5


    black = sgcolor('black')
    bad_color = sgcolor('red')
    good_color = sgcolor('blue')
    dim_index = [0,1]
    
    v2d_angle = 15.
    dis_scale = 4.
    vel_scale = 40


;---Collect info.
    triad_list = event.triad_list
    combo_keys = list()
    foreach triad, triad_list do combo_keys.add, strjoin(sort_uniq(triad.probes),'_')
    
    df_list = event.df_list
    df_dict = dictionary()
    foreach df, df_list do df_dict[df.probe] = df
    probe_xys = [-1,1]
    foreach df, df_list do probe_xys = [probe_xys,df.obs_rxy]
    xy_range = minmax(make_bins(probe_xys,5))
    
    probes = event.probes
    nvertex = 3
    probe_combos = choose_from(probes,nvertex)
    nprobe_combo = n_elements(probe_combos)
    region_name = (strsplit(event.region,'%',/extract))[0]


;---Plot settings.
    label_letters = letters(nprobe_combo)
    nxpanel = 4
    nypanel = ceil(nprobe_combo/nxpanel)
    panel_xrange = -xy_range
    panel_yrange = (region_name eq 'post_midn')? -xy_range: reverse(xy_range)
    panel_xsize = 2
    panel_ysize = panel_xsize*total(panel_yrange*[-1,1])/total(panel_xrange*[-1,1])
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    margins = [7.,4,2,2]
    xpad = 1*abs_ychsz/abs_xchsz
    ypad = 1
    fig_xsize = panel_xsize*nxpanel+(xpad*(nxpanel-1)+(margins[0]+margins[2]))*abs_xchsz
    fig_ysize = panel_ysize*nypanel+(ypad*(nypanel-1)+(margins[1]+margins[3]))*abs_ychsz
    arrow_hsize = keyword_set(test)? 4:80


    event_id = time_string(event.time_range[0],tformat='YYYY_MMDD_hh')
    file = join_path([project.plot_dir,'triangle_collage','fig_'+event_id+'_triangle_collage.pdf'])
    if keyword_set(test) then file = 0
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify, /inch
    panel_poss = sgcalcpos(nypanel,nxpanel, xpad=xpad, ypad=ypad, margins=margins, xchsz=xchsz, ychsz=ychsz)


;---Treat each combo.
    foreach probes, probe_combos, panel_id do begin
        probe_combo = strjoin(probes[sort(probes)],'_')
        nprobe = n_elements(probes)
        ref_times = dblarr(nprobe)
        rsms = fltarr(3,nprobe)
        foreach probe, probes, ii do begin
            ref_times[ii] = df_dict[probe].obs_time
            rsms[*,ii] = df_dict[probe].obs_r_sm
        endforeach
        index = sort(ref_times)
        ref_times = ref_times[index]
        rsms = rsms[*,index]
        probes = probes[index]  ; need to sort by ref_time, cannt read from event_info.
        rsm_center = total(rsms,2)/nvertex
        
        the_key = strjoin(sort_uniq(probes),'_')
        index = combo_keys.where(the_key)
        if index ne !null then begin
            good_combo = 1
            v_mag = triad_list[index[0]].vmag_obs_time
            v_hat = triad_list[index[0]].vhat_obs_time
        endif else begin
            good_combo = 0
            v_mag = 0
            v_hat = [0,0]
        endelse
        

        index = array_indices([nxpanel,nypanel], panel_id, /dimensions)
        tpos = panel_poss[*,index[0],index[1]]
        xtitle = (index[1] eq nypanel-1)? 'SM X (Re)': ''
        ytitle = (index[0] eq 0)? 'SM Y (Re)': ''
        xtickformat = (index[1] eq nypanel-1)? '(I0)': '(A1)'
        ytickformat = (index[0] eq 0)? '(I0)': '(A1)'
        plot, panel_xrange, panel_yrange, xstyle=1, ystyle=1, position=tpos, $
            xrange=panel_xrange, yrange=panel_yrange, /nodata, /noerase, /isotropic, $
            xtitle=xtitle, ytitle=ytitle, xticklen=-0.01, yticklen=-0.01, xtickformat=xtickformat, ytickformat=ytickformat

       

        ; Add geometry and timing criteria, add flag.
        msg = good_combo? 'Good': 'Bad'
        triad_color = good_combo? good_color: bad_color
        tx = tpos[2]-xchsz*1
        ty = tpos[1]+ychsz*0.5
        xyouts, tx,ty,/normal,alignment=1, msg, color=triad_color, charsize=large_size

        ; Add earth.
        nangle = 50
        angles = smkarthm(0,2*!dpi,nangle,'n')
        circle_x = cos(angles)
        circle_y = sin(angles)
        polyfill, circle_x<0, circle_y, color=sgcolor('silver')
        plots, circle_x, circle_y
        foreach r, make_bins(minmax(abs(xy_range)),5,/inner) do oplot, circle_x*r, circle_y*r, linestyle=1


        ; Add spacecraft and labels.
        foreach probe, probes, ii do begin
            tx = rsms[dim_index[0],ii]
            ty = rsms[dim_index[1],ii]
            plots, tx,ty,/data, psym=6, symsize=label_size*0.5

            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            if probe eq 'thd' then begin
                tx = tx-xchsz*1.5
                ty = ty+ychsz*0.3
            endif
            probe_info = resolve_probe(probe)
            short_name = probe_info.short_name
            xyouts, tx+xchsz*0.5,ty,/normal, strupcase(short_name), charsize=label_size
        endforeach

        ; Add triangle.
        tx = rsms[dim_index[0],*]
        ty = rsms[dim_index[1],*]
        for ii=0, nprobe-1 do begin
            plots, (shift(tx,ii))[0:1], (shift(ty,ii))[0:1], linestyle=0, color=sgcolor('silver')
        endfor

        ; Add arrow.
        x0 = mean(rsms[dim_index[0],*])
        y0 = mean(rsms[dim_index[1],*])
        scale = dis_scale/vel_scale*v_mag
        if scale ne 0 then begin
            x1 = x0+v_hat[0]*scale
            y1 = y0+v_hat[1]*scale
            arrow, x0,y0,x1,y1,/data, /solid, hsize=arrow_hsize*2
        endif
        

        ; Add figure label.
        ;fig_label = '('+sgnum2str(index[1]+1)+'-'+sgnum2str(index[0]+1)+')'
        fig_label = 'e-'+sgnum2str(index[0]+1)+'. '
        tx = tpos[0]+xchsz*0
        triangles = triangle_angles(reform(rsms,[1,3,3]))
        triangles = sort_uniq(triangles)
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty, /normal, fig_label+'Angle (deg):'
        foreach triangle, triangles, jj do begin
            msg = strtrim(string(triangle,format='(F5.1)'),2)
            if jj ne 2 then msg += ','
            bad_flag = triangle le v2d_angle or triangle ge (180-v2d_angle)
            color = bad_flag? bad_color: black
            ttx = tx+(11.5+3.5*jj)*xchsz
            xyouts, ttx,ty,/normal, msg, color=color
        endforeach


;        ; Add time info.
;        ty = tpos[3]-ychsz*2
;        time_diffs = abs([ref_times[1]-ref_times[0],ref_times[2]-ref_times[1],ref_times[0]-ref_times[2]])
;        time_diffs = sort_uniq(time_diffs)
;        xyouts, tx,ty, /normal, 'T.diff (sec):', charsize=label_size
;        foreach time_diff, time_diffs, jj do begin
;            msg = strtrim(string(round(time_diff),format='(I0)'),2)
;            if jj ne 2 then msg += ','
;            ;bad_flag = time_diff le v2d_dt
;            ;color = bad_flag? bad_color: black
;            color = black
;            ttx = tx+(6+3*jj)*xchsz
;            xyouts, ttx,ty,/normal, msg, charsize=label_size, color=color
;        endforeach


;        ; Add velocity info.
;        tx = tpos[2]+xchsz*0
;        ty = tpos[1]+ychsz*1
;        xyouts, tx,ty, alignment=0, /normal, charsize=label_size, '|V| = '+string(v_mag,format='(F5.1)')+' km/s'

        
        ; Add scale.
        tx = tpos[2]-xchsz*2
        ty0 = tpos[3]-ychsz*2.5
        ty = ty0
        tmp = convert_coord(tx,ty, /normal, /to_data)
        tmp = convert_coord(tmp[0]-dis_scale,tmp[1], /data, /to_normal)
        dx = abs(tmp[0]-tx)
        xs = tx+[-1,0]*dx
        plots, xs, ty+[0,0], /normal
        foreach tx, xs do plots, tx+[0,0], ty+[-1,1]*ychsz*0.1, /normal
        tx = mean(xs)
        ty = ty0+ychsz*0.5
        xyouts, tx,ty,/normal, alignment=0.5, sgnum2str(vel_scale)+' km/s'


        if keyword_set(test) then stop
    endforeach

    if keyword_set(test) then stop
    sgclose

end


;---Test events.
    test_event_times = time_double([$
        '2008-01-09/11:27:45' ])

    project = azim_df_load_project()
    candidates = azim_df_find_dfgroup(project=project)

;---Select test_list.
    test_candidates = list()
    foreach test_event_time, test_event_times do begin
        lprmsg, 'Find event '+time_string(test_event_time)+' ...'
        foreach candidate, candidates do begin
            if candidate.time_range[0] eq test_event_time then begin
                lprmsg, 'Found the wanted event ...'
                test_candidates.add, candidate
                break
            endif
        endforeach
    endforeach


    foreach candidate, test_candidates do begin
        fig_triad_collage, candidate
    endforeach

end
