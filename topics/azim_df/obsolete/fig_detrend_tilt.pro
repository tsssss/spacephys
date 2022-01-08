;+
; Create a plot to show the steps to detrend the tilt angle.
;-

pro fig_detrend_tilt_add_legend, position=tpos, colors=colors, legends=legends, $
    charsize=charsize

    if n_elements(xchsz) eq 0 then xchsz = double(!d.x_ch_size)/!d.x_size
    if n_elements(ychsz) eq 0 then ychsz = double(!d.y_ch_size)/!d.y_size
    if n_elements(charsize) eq 0 then charsize = 1.0
    if n_elements(lmargin) eq 0 then lmargin = 0.5
    if n_elements(tmargin) eq 0 then tmargin = 1.0

    nlegend = n_elements(legends)
    tx = tpos[0]+lmargin*xchsz
    ty0 = tpos[3]-tmargin*ychsz
    line_length = 3

    foreach legend, legends, ii do begin
        ty = ty0-ychsz*ii
        plots, tx+[0,line_length]*xchsz, ty+0.35*ychsz*charsize+[0,0], /normal, color=colors[ii]
        xyouts, tx+xchsz*(line_length+0.5),ty,/normal, legend, charsize=charsize, color=colors[ii]
    endforeach
end

pro fig_detrend_tilt, time_range, event_id=event_id, project=project, mission_probe=mission_probe

test = 0
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(time_range) eq 2 then begin
        event_id = time_string(mean(time_range),tformat='YYYY_MMDD_hh')
    endif else if n_elements(event_id) eq 0 then begin
        errmsg = handle_error('No input time_range or event_id ...')
        return
    endif
    if n_elements(mission_probe) ne 1 then begin
        errmsg = handle_error('Invalid mission_probe ...')
        return
    endif
    prefix = mission_probe+'_'

    data_file = join_path([project.data_dir,event_id+'_all_data.tplot'])
    if file_test(data_file) eq 0 then azim_df_load_basic_data, event_id=event_id, project=project
    tplot_restore, filename=data_file

    if n_elements(plot_time_range) ne 2 then begin
        vars = tnames('*_r_sm')
        get_data, vars[0], common_times
        plot_time_range = minmax(common_times)
    endif
    detrend_time = project['detrend_time']
    if n_elements(time_range) eq 2 then begin
        plot_time_range = time_range+[-1,1]*detrend_time*0.5
    endif


    file = join_path([project.plot_dir,'fig_detrend_tilt_'+event_id+'_'+mission_probe+'.pdf'])
    if keyword_set(test) then file = test

    sgopen, file, xsize=5, ysize=5, /inch
    nvar = 3
    poss = sgcalcpos(nvar, tmargin=2, rmargin=2, xchsz=xchsz, ychsz=ychsz)
    labels = letters(nvar)+'.'
    foreach label, labels, ii do begin
        tx = 2*xchsz
        ty = poss[3,ii]-ychsz*0.7
        xyouts, tx,ty,/normal, label
    endforeach
    tx = poss[0,0]
    ty = poss[3,0]+ychsz*0.35
    xyouts, tx,ty,/normal, 'Detrending '+strupcase(mission_probe)+' tilt angle data'
    

    ; The measured tilt angle.
    the_var = prefix+'alpha'
    stplot_merge, prefix+['b_tilt','bmod_tilt'], newname=the_var
    add_setting, the_var, /smart, {$
        display_type: 'vector', $
        short_name: '', $
        unit: 'deg', $
        coord: '', $
        coord_labels: ['',''], $
        colors: sgcolor(['black','red'])}
    yrange = minmax(get_var_data(the_var, in=plot_time_range))
    yrange = minmax(make_bins(yrange,2))
    yticks = 2
    yminor = total(yrange*[-1,1])/2
    foreach factor, [4,3,5,2] do if (yminor mod factor) eq 0 then begin
        yminor = factor
        break
    endif
    add_setting, the_var, {$
        yrange: yrange, $
        yticks: yticks, $
        ytickv: smkarthm(yrange[0],yrange[1],yticks+1,'n'), $
        yminor: yminor}
        
    tpos = poss[*,0]
    tplot, the_var, position=tpos, /noerase, /nouttick, trange=plot_time_range
    legends = tex2str('alpha')+'!D'+['measured','T89']
    colors = get_setting(the_var, 'colors')
    fig_detrend_tilt_add_legend, position=tpos, colors=colors, legends=legends


    ; Showing the detrending.
    the_var = prefix+'dalpha'
    sys_subtract, prefix+'db_tilt', prefix+'theta', to=the_var
    stplot_merge, prefix+['db_tilt','dalpha'], newname=the_var, labels=['','']
    add_setting, the_var, /smart, {$
        display_type: 'vector', $
        short_name: '', $
        unit: 'deg', $
        coord: '', $
        coord_labels: ['',''], $
        colors: sgcolor(['black','blue'])}
    yrange = minmax(get_var_data(the_var, in=plot_time_range))
    yrange = minmax(make_bins(yrange,2))
    yticks = 2
    yminor = total(yrange*[-1,1])/2
    foreach factor, [4,3,5,2] do if (yminor mod factor) eq 0 then begin
        yminor = factor
        break
    endif
    add_setting, the_var, {$
        yrange: yrange, $
        yticks: yticks, $
        ytickv: smkarthm(yrange[0],yrange[1],yticks+1,'n'), $
        yminor: yminor}
        
    tpos = poss[*,1]
    tplot, the_var, position=tpos, /noerase, /nouttick, trange=plot_time_range
    legends = [$
        tex2str('Delta')+tex2str('alpha')+' = '+tex2str('alpha')+'!Dmeasured!N - '+tex2str('alpha')+'!DT89', $
        tex2str('Delta')+tex2str('alpha')+' smoothed over '+sgnum2str(detrend_time/60)+' min']
    colors = get_setting(the_var, 'colors')
    fig_detrend_tilt_add_legend, position=tpos, colors=colors, legends=legends
    tx = mean(tpos[[0,2]])
    ty = tpos[1]+ychsz*0.35
    plot, plot_time_range, [0,1], /noerase, /nodata, position=tpos, xstyle=5, ystyle=5
    x1 = convert_coord(time_range[0],0, /data, /to_normal)
    x2 = convert_coord(time_range[0]+detrend_time,0, /data, /to_normal)
    xs = tx+[-1,1]*(x2[0]-x1[0])*0.5
    plots, xs, ty+[0,0],/normal, color=colors[1]
    foreach tmp, xs do plots, tmp+[0,0], ty+[-1,1]*ychsz*0.1, /normal, color=colors[1]
    ty = ty+ychsz*0.35
    xyouts, tx, ty, /normal, alignment=0.5, 'Window size = '+sgnum2str(detrend_time/60)+' min', color=colors[1]


    ; Showing the final result.
    tpos = poss[*,2]
    the_var = prefix+'theta'
    options, the_var, 'constant', 0
    options, the_var, 'labels', ''
    yrange = minmax(get_var_data(the_var, in=plot_time_range))
    yrange = [-1,1]*ceil(max(abs(yrange)))
    yticks = 2
    yminor = total(yrange*[-1,1])/2
    foreach factor, [4,3,5,2] do if (yminor mod factor) eq 0 then begin
        yminor = factor
        break
    endif
    add_setting, the_var, {$
        yrange: yrange, $
        yticks: yticks, $
        ytickv: smkarthm(yrange[0],yrange[1],yticks+1,'n'), $
        yminor: yminor}
    
    tplot, the_var, position=tpos, /noerase, trange=plot_time_range
    legends = tex2str('theta')+' the detrended tilt angle'
    colors = sgcolor('black')
    fig_detrend_tilt_add_legend, position=tpos, colors=colors, legends=legends

    if keyword_set(test) then stop
    sgclose

end

time_range = time_double(['2014-08-28/09:30','2014-08-28/11:30'])
mission_probe = 'rbspb'
fig_detrend_tilt, time_range, project=project, mission_probe=mission_probe
end
