pro rbsp_survey_on_ebfield_gen_plot, time_range, filename=plot_file, project=project, probe=probe
    
;---Internal do not check inputs.
    ; Panels: dEuv, |dE|, Vsc, dBxyz, |dB|
    ypans = [2.,1,1,2,1]
    nypanel = n_elements(ypans)
    margins = [8,3,8,2]
    panel_ypad = 0.5

    fig_xsize = 8.
    fig_ysize = 11.
    sgopen, 0, xsize=fig_xsize, ysize=fig_ysize, /inch
    poss = sgcalcpos(nypanel, ypans=ypans, ypad=panel_ypad, margins=margins, xchsz=xchsz, ychsz=ychsz)
    sgclose, /wdelete

    if n_elements(plot_file) eq 0 then plot_file = 0
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch


;---Settings for plot.
    ; Common x-axis.
    xrange = time_range
    xstep = constant('secofhour')*4
    xminor = 8
    xtickv = make_bins(xrange, xstep)
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    for ii=0, xticks do xtickn[ii] = time_string(xtickv[ii],tformat='hh:mm')
    xtickn[0] = time_string(time_range[0],tformat='YYYY-MM-DD')+'      '

    ; Other settings.
    xticklen_chsz = -0.10
    yticklen_chsz = -0.40
    full_ychsz = constant('full_ychsz')
    half_ychsz = constant('half_ychsz')
    panel_label_size = 1.0
    label_size = constant('label_size')
    label_xshift = 6
    bar_thick = (size(plot_file,/type) eq 7)? 8: 2
    prefix = 'rbsp'+probe+'_'



;---Add labels.
    panel_letters = letters(nypanel)+'.'
    panel_letter_xshift = 6
    for ii=0,nypanel-1 do begin
        tpos = poss[*,ii]
        tx = tpos[0]-panel_letter_xshift*xchsz
        ty = tpos[3]-full_ychsz*ychsz
        xyouts, tx,ty,/normal, panel_letters[ii]
    endfor

    ; Overall title.
    msg = 'RBSP'+strupcase(probe)+': E and B field on '+time_string(time_range[0],tformat='YYYY-MM-DD')
    tpos = poss[*,0]
    tx = tpos[0]
    ty = tpos[3]+ychsz*half_ychsz
    xyouts, tx,ty,/normal, msg

    ; Perigee.
    perigee_dis = 4.
    var = prefix+'r_gsm'
    orbit_time_step = 60.
    color = sgcolor('red')
    dis = snorm(get_var_data(var, in=time_range, times=times))
    index = where(dis le perigee_dis, count)
    perigee_time_ranges = time_to_range(times[index], time_step=orbit_time_step)
    tpos = poss[*,0]
    tpos[1] = min(poss[1,*])
    tpos[3] = max(poss[3,*])
    yrange = [0,1]
    plot, xrange, yrange, /nodata, /noerase, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos
    foreach time, perigee_time_ranges do begin
        oplot, time+[0,0], yrange, linestyle=1, color=color
    endforeach


;---Add flags.
    tpos = poss[*,1]
    tpos[1] = poss[3,1]
    tpos[3] = poss[1,0]
    tx = tpos[2]+xchsz*0.5
    ty = mean(tpos[[1,3]])
    color = sgcolor('blue')
    xyouts, tx,ty,/normal, 'SDT & eclipse', charsize=label_size, color=color
    var = prefix+'eclipse_flag'
    yys = get_var_data(var, in=time_range, times=xxs)
    var = prefix+'sdt_flag'
    yys += get_var_data(var, at=xxs)
    index = where(yys eq 1, count)
    flag_time_step = 60.
    if count ne 0 then begin
        yrange = [0,1]
        plot, xrange, yrange, /nodata, /noerase, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos
        flag_time_ranges = time_to_range(times[index], time_step=flag_time_step)
        nflag_time_range = n_elements(flag_time_ranges)/2
        for ii=0,nflag_time_range-1 do plots, flag_time_ranges[ii,*], mean(yrange)+[0,0], thick=bar_thick, color=color
    endif


;---Panel 1. dE_uv.
    var = prefix+'de_uv'
    xtickformat = '(A1)'
    efield_small_range = [-1,1]*10.
    efield_large_range = efield_small_range*30
    panel_pos = poss[*,0]
    subpanel_ypans = [1.,2,1]
    nsubpanel = n_elements(subpanel_ypans)
    subpanel_poss = sgcalcpos(nsubpanel, position=panel_pos, ypad=0, ypans=subpanel_ypans)
    panel_info = dictionary()

    ; Linear part.
    tpos = subpanel_poss[*,1]
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yrange = efield_small_range
    ystep = 4
    ytickv = make_bins(yrange, ystep, /inner)
    yminor = 2
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 0
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['linear_part'] = subpanel_info

    ; Positive log part.
    tpos = subpanel_poss[*,0]
    yrange = [efield_small_range[1],efield_large_range[1]]
    ytickv = efield_small_range[1]*[1,10]
    yminor = 10
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 1
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['positive_part'] = subpanel_info

    ; Negative log part.
    tpos = subpanel_poss[*,2]
    yrange = -[efield_large_range[0],efield_small_range[0]]
    ytickv = -efield_small_range[0]*[10,1]
    ytickn = '-'+string(ytickv,format='(I0)')
    yminor = 10
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 1
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['ytickname'] = ytickn
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['negative_part'] = subpanel_info

    ; Plot data.
    yys = get_var_data(var, in=xrange, times=xxs, limits=lim)
    if ~stagexist('colors',lim) then begin
        labels = lim.labels
        nlabel = n_elements(labels)
        case nlabel of
            1: colors = sgcolor('black')
            2: colors = sgcolor(['black','blue'])
            3: colors = constant('rgb')
            else: begin
                ct = 40
                colors = smkarthm(50,200,nlabel,'n')
                for ii=0, nlabel-1 do colors[ii] = sgcolor(colors[ii], ct=ct)
            end
        endcase
        lim = create_struct('colors', colors, lim)
    endif
    panel_keys = panel_info.keys()
    foreach key, panel_keys do begin
        subpanel_info = panel_info[key]
        plot, xrange, subpanel_info.yrange, /nodata, _extra=subpanel_info.tostruct()
        if get_setting(var, 'display_type') eq 'scalar' then begin
            the_yys = yys
            if key eq 'negative_part' then the_yys = -the_yys
            oplot, xxs, the_yys, color=lim.colors[0]
        endif else begin
            ndim = n_elements(lim.labels)
            for label_id=0,ndim-1 do begin
                the_yys = yys[*,label_id]
                if key eq 'negative_part' then the_yys = -the_yys
                oplot, xxs, the_yys, color=lim.colors[label_id]
            endfor
        endelse
    endforeach

    ; Draw other parts of the box.
    yrange = [0,1]
    ytitle = lim.ytitle
    plot, xrange, yrange, /nodata, /noerase, $
        xstyle=1, xrange=xrange, xtickformat=xtickformat, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
        ystyle=5, yrange=yrange, $
        position=panel_pos
    nsubpanel = n_elements(panel_keys)
    for panel_id=0,nsubpanel-2 do begin
        tpos = subpanel_poss[*,panel_id]
        plots, tpos[[0,2]], tpos[1]+[0,0], /normal, linestyle=1
    endfor

    ; Ytitle.
    tpos = panel_pos
    tx = tpos[0]-xchsz*4
    ty = total(tpos[[1,3]])*0.5-ychsz*half_ychsz
    xyouts, tx,ty,/normal, alignment=1, lim.ytitle, orientation=90

    ; Labels.
    labels = lim.labels
    colors = lim.colors
    nlabel = n_elements(labels)
    panel_ysize = total(tpos[[1,3]]*[-1,1])
    dy = panel_ysize/(nlabel+1)
    label_ys = tpos[1]+dy+findgen(nlabel)*dy-ychsz*panel_label_size*half_ychsz
    for ii=0, nlabel-1 do begin
        tx = tpos[2]+xchsz*1
        ty = label_ys[ii]
        xyouts, tx,ty,/normal, labels[ii], color=colors[ii], charsize=panel_label_size
    endfor


;---Panel 2. |dE|.
    var = prefix+'demag'
    xtickformat = '(A1)'
    efield_small_range = [-1,1]*10.
    efield_large_range = efield_small_range*30
    panel_pos = poss[*,1]
    subpanel_ypans = [1.,1]
    nsubpanel = n_elements(subpanel_ypans)
    subpanel_poss = sgcalcpos(nsubpanel, position=panel_pos, ypad=0, ypans=subpanel_ypans)
    panel_info = dictionary()
    ; Linear part.
    tpos = subpanel_poss[*,1]
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yrange = efield_small_range>0
    ystep = 4
    ytickv = make_bins(yrange, ystep, /inner)
    yminor = 2
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 0
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['linear_part'] = subpanel_info

    ; Positive log part.
    tpos = subpanel_poss[*,0]
    yrange = [efield_small_range[1],efield_large_range[1]]
    ytickv = efield_small_range[1]*[1,10]
    yminor = 10
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 1
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['positive_part'] = subpanel_info

    ; Plot data.
    yys = get_var_data(var, in=xrange, times=xxs, limits=lim)
    if ~stagexist('colors',lim) then begin
        labels = lim.labels
        nlabel = n_elements(labels)
        case nlabel of
            1: colors = sgcolor('black')
            2: colors = sgcolor(['black','blue'])
            3: colors = constant('rgb')
            else: begin
                ct = 40
                colors = smkarthm(50,200,nlabel,'n')
                for ii=0, nlabel-1 do colors[ii] = sgcolor(colors[ii], ct=ct)
            end
        endcase
        lim = create_struct('colors', colors, lim)
    endif
    panel_keys = panel_info.keys()
    foreach key, panel_keys do begin
        subpanel_info = panel_info[key]
        plot, xrange, subpanel_info.yrange, /nodata, _extra=subpanel_info.tostruct()
        if get_setting(var, 'display_type') eq 'scalar' then begin
            oplot, xxs, yys, color=lim.colors[0]
        endif else begin
            ndim = n_elements(lim.labels)
            for label_id=0,ndim-1 do oplot, xxs, yys[*,label_id], color=lim.colors[label_id]
        endelse
    endforeach

    ; Draw other parts of the box.
    yrange = [0,1]
    ytitle = lim.ytitle
    plot, xrange, yrange, /nodata, /noerase, $
        xstyle=1, xrange=xrange, xtickformat=xtickformat, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
        ystyle=5, yrange=yrange, $
        position=panel_pos
    nsubpanel = n_elements(panel_keys)
    for panel_id=0,nsubpanel-2 do begin
        tpos = subpanel_poss[*,panel_id]
        plots, tpos[[0,2]], tpos[1]+[0,0], /normal, linestyle=1
    endfor

    ; Ytitle.
    tpos = panel_pos
    tx = tpos[0]-xchsz*4
    ty = total(tpos[[1,3]])*0.5-ychsz*half_ychsz
    xyouts, tx,ty,/normal, alignment=1, lim.ytitle, orientation=90

    ; Labels.
    labels = lim.labels
    colors = lim.colors
    nlabel = n_elements(labels)
    panel_ysize = total(tpos[[1,3]]*[-1,1])
    dy = panel_ysize/(nlabel+1)
    label_ys = tpos[1]+dy+findgen(nlabel)*dy-ychsz*panel_label_size*half_ychsz
    for ii=0, nlabel-1 do begin
        tx = tpos[2]+xchsz*1
        ty = label_ys[ii]
        xyouts, tx,ty,/normal, labels[ii], color=colors[ii], charsize=panel_label_size
    endfor


;---Panel 3. Vsc.
    var = prefix+'vsc'
    xtickformat = '(A1)'
    vsc_small_range = [-1,1]*5
    vsc_large_range = vsc_small_range*20
    panel_pos = poss[*,2]
    subpanel_ypans = [1.,2,1]
    nsubpanel = n_elements(subpanel_ypans)
    subpanel_poss = sgcalcpos(nsubpanel, position=panel_pos, ypad=0, ypans=subpanel_ypans)
    panel_info = dictionary()

    ; Linear part.
    tpos = subpanel_poss[*,1]
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yrange = vsc_small_range
    ystep = 4
    ytickv = make_bins(yrange, ystep, /inner)
    yminor = 4
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 0
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['linear_part'] = subpanel_info

    ; Positive log part.
    tpos = subpanel_poss[*,0]
    yrange = [vsc_small_range[1],vsc_large_range[1]]
    ytickv = vsc_small_range[1]*[1,10]
    yminor = 10
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 1
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['positive_part'] = subpanel_info

    ; Negative log part.
    tpos = subpanel_poss[*,2]
    yrange = -[vsc_large_range[0],vsc_small_range[0]]
    ytickv = -vsc_small_range[0]*[10,1]
    ytickn = '-'+string(ytickv,format='(I0)')
    yminor = 10
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 1
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['ytickname'] = ytickn
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['negative_part'] = subpanel_info

    ; Plot data.
    yys = get_var_data(var, in=xrange, times=xxs, limits=lim)
    if ~stagexist('colors',lim) then begin
        labels = lim.labels
        nlabel = n_elements(labels)
        case nlabel of
            1: colors = sgcolor('black')
            2: colors = sgcolor(['black','blue'])
            3: colors = constant('rgb')
            else: begin
                ct = 40
                colors = smkarthm(50,200,nlabel,'n')
                for ii=0, nlabel-1 do colors[ii] = sgcolor(colors[ii], ct=ct)
            end
        endcase
        lim = create_struct('colors', colors, lim)
    endif
    panel_keys = panel_info.keys()
    foreach key, panel_keys do begin
        subpanel_info = panel_info[key]
        plot, xrange, subpanel_info.yrange, /nodata, _extra=subpanel_info.tostruct()
        if get_setting(var, 'display_type') eq 'scalar' then begin
            the_yys = yys
            if key eq 'negative_part' then the_yys = -the_yys
            oplot, xxs, the_yys, color=lim.colors[0]
        endif else begin
            ndim = n_elements(lim.labels)
            for label_id=0,ndim-1 do begin
                the_yys = yys[*,label_id]
                if key eq 'negative_part' then the_yys = -the_yys
                oplot, xxs, the_yys, color=lim.colors[label_id]
            endfor
        endelse
    endforeach

    ; Draw other parts of the box.
    yrange = [0,1]
    ytitle = lim.ytitle
    plot, xrange, yrange, /nodata, /noerase, $
        xstyle=1, xrange=xrange, xtickformat=xtickformat, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
        ystyle=5, yrange=yrange, $
        position=panel_pos
    nsubpanel = n_elements(panel_keys)
    for panel_id=0,nsubpanel-2 do begin
        tpos = subpanel_poss[*,panel_id]
        plots, tpos[[0,2]], tpos[1]+[0,0], /normal, linestyle=1
    endfor

    ; Ytitle.
    tpos = panel_pos
    tx = tpos[0]-xchsz*4
    ty = total(tpos[[1,3]])*0.5-ychsz*half_ychsz
    xyouts, tx,ty,/normal, alignment=1, lim.ytitle, orientation=90

    ; Labels.
    labels = lim.labels
    colors = lim.colors
    nlabel = n_elements(labels)
    panel_ysize = total(tpos[[1,3]]*[-1,1])
    dy = panel_ysize/(nlabel+1)
    label_ys = tpos[1]+dy+findgen(nlabel)*dy-ychsz*panel_label_size*half_ychsz
    for ii=0, nlabel-1 do begin
        tx = tpos[2]+xchsz*1
        ty = label_ys[ii]
        xyouts, tx,ty,/normal, labels[ii], color=colors[ii], charsize=panel_label_size
    endfor


;---Panel 4. dB_xyz.
    var = prefix+'db_gsm'
    xtickformat = '(A1)'
    small_range = [-1,1]*50.
    large_range = small_range*10
    panel_pos = poss[*,3]
    subpanel_ypans = [1.,2,1]
    nsubpanel = n_elements(subpanel_ypans)
    subpanel_poss = sgcalcpos(nsubpanel, position=panel_pos, ypad=0, ypans=subpanel_ypans)
    panel_info = dictionary()

    ; Linear part.
    tpos = subpanel_poss[*,1]
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yrange = small_range
    ystep = 20
    ytickv = make_bins(yrange, ystep, /inner)
    yminor = 2
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 0
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['linear_part'] = subpanel_info

    ; Positive log part.
    tpos = subpanel_poss[*,0]
    yrange = [small_range[1],large_range[1]]
    ytickv = small_range[1]*[1,10]
    yminor = 10
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 1
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['positive_part'] = subpanel_info

    ; Negative log part.
    tpos = subpanel_poss[*,2]
    yrange = -[large_range[0],small_range[0]]
    ytickv = -small_range[0]*[10,1]
    ytickn = '-'+string(ytickv,format='(I0)')
    yminor = 10
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 1
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['ytickname'] = ytickn
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['negative_part'] = subpanel_info

    ; Plot data.
    yys = get_var_data(var, in=xrange, times=xxs, limits=lim)
    if ~stagexist('colors',lim) then begin
        labels = lim.labels
        nlabel = n_elements(labels)
        case nlabel of
            1: colors = sgcolor('black')
            2: colors = sgcolor(['black','blue'])
            3: colors = constant('rgb')
            else: begin
                ct = 40
                colors = smkarthm(50,200,nlabel,'n')
                for ii=0, nlabel-1 do colors[ii] = sgcolor(colors[ii], ct=ct)
            end
        endcase
        lim = create_struct('colors', colors, lim)
    endif
    panel_keys = panel_info.keys()
    foreach key, panel_keys do begin
        subpanel_info = panel_info[key]
        plot, xrange, subpanel_info.yrange, /nodata, _extra=subpanel_info.tostruct()
        if get_setting(var, 'display_type') eq 'scalar' then begin
            the_yys = yys
            if key eq 'negative_part' then the_yys = -the_yys
            oplot, xxs, the_yys, color=lim.colors[0]
        endif else begin
            ndim = n_elements(lim.labels)
            for label_id=0,ndim-1 do begin
                the_yys = yys[*,label_id]
                if key eq 'negative_part' then the_yys = -the_yys
                oplot, xxs, the_yys, color=lim.colors[label_id]
            endfor
        endelse
    endforeach

    ; Draw other parts of the box.
    yrange = [0,1]
    ytitle = lim.ytitle
    plot, xrange, yrange, /nodata, /noerase, $
        xstyle=1, xrange=xrange, xtickformat=xtickformat, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
        ystyle=5, yrange=yrange, $
        position=panel_pos
    nsubpanel = n_elements(panel_keys)
    for panel_id=0,nsubpanel-2 do begin
        tpos = subpanel_poss[*,panel_id]
        plots, tpos[[0,2]], tpos[1]+[0,0], /normal, linestyle=1
    endfor

    ; Ytitle.
    tpos = panel_pos
    tx = tpos[0]-xchsz*4
    ty = total(tpos[[1,3]])*0.5-ychsz*half_ychsz
    xyouts, tx,ty,/normal, alignment=1, lim.ytitle, orientation=90

    ; Labels.
    labels = lim.labels
    colors = lim.colors
    nlabel = n_elements(labels)
    panel_ysize = total(tpos[[1,3]]*[-1,1])
    dy = panel_ysize/(nlabel+1)
    label_ys = tpos[1]+dy+findgen(nlabel)*dy-ychsz*panel_label_size*half_ychsz
    for ii=0, nlabel-1 do begin
        tx = tpos[2]+xchsz*1
        ty = label_ys[ii]
        xyouts, tx,ty,/normal, labels[ii], color=colors[ii], charsize=panel_label_size
    endfor



;---Panel 5. |dB|.
    var = prefix+'dbmag'
    xtickformat = ''
    small_range = [-1,1]*50.
    large_range = small_range*10
    panel_pos = poss[*,4]
    subpanel_ypans = [1.,1]
    nsubpanel = n_elements(subpanel_ypans)
    subpanel_poss = sgcalcpos(nsubpanel, position=panel_pos, ypad=0, ypans=subpanel_ypans)
    panel_info = dictionary()
    ; Linear part.
    tpos = subpanel_poss[*,1]
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yrange = small_range>0
    ystep = 20
    ytickv = make_bins(yrange, ystep, /inner)
    yminor = 2
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 0
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['linear_part'] = subpanel_info

    ; Positive log part.
    tpos = subpanel_poss[*,0]
    yrange = [small_range[1],large_range[1]]
    ytickv = small_range[1]*[1,10]
    yminor = 10
    subpanel_info = dictionary()
    subpanel_info['position'] = tpos
    subpanel_info['noerase'] = 1
    subpanel_info['ystyle'] = 1
    subpanel_info['yrange'] = yrange
    subpanel_info['ytickv'] = ytickv
    subpanel_info['yticks'] = n_elements(ytickv)-1
    subpanel_info['yminor'] = yminor
    subpanel_info['ylog'] = 1
    subpanel_info['ytitle'] = ' '
    subpanel_info['yticklen'] = yticklen
    subpanel_info['xticklen'] = xticklen
    subpanel_info['xtickformat'] = '(A1)'
    subpanel_info['xstyle'] = 5
    panel_info['positive_part'] = subpanel_info

    ; Plot data.
    yys = get_var_data(var, in=xrange, times=xxs, limits=lim)
    if ~stagexist('colors',lim) then begin
        labels = lim.labels
        nlabel = n_elements(labels)
        case nlabel of
            1: colors = sgcolor('black')
            2: colors = sgcolor(['black','blue'])
            3: colors = constant('rgb')
            else: begin
                ct = 40
                colors = smkarthm(50,200,nlabel,'n')
                for ii=0, nlabel-1 do colors[ii] = sgcolor(colors[ii], ct=ct)
            end
        endcase
        lim = create_struct('colors', colors, lim)
    endif
    panel_keys = panel_info.keys()
    foreach key, panel_keys do begin
        subpanel_info = panel_info[key]
        plot, xrange, subpanel_info.yrange, /nodata, _extra=subpanel_info.tostruct()
        if get_setting(var, 'display_type') eq 'scalar' then begin
            the_yys = yys
            if key eq 'negative_part' then the_yys = -the_yys
            oplot, xxs, the_yys, color=lim.colors[0]
        endif else begin
            ndim = n_elements(lim.labels)
            for label_id=0,ndim-1 do begin
                the_yys = yys[*,label_id]
                if key eq 'negative_part' then the_yys = -the_yys
                oplot, xxs, the_yys, color=lim.colors[label_id]
            endfor
        endelse
    endforeach

    ; Draw other parts of the box.
    yrange = [0,1]
    ytitle = lim.ytitle
    plot, xrange, yrange, /nodata, /noerase, $
        xstyle=1, xrange=xrange, xtickformat=xtickformat, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
        ystyle=5, yrange=yrange, $
        position=panel_pos
    nsubpanel = n_elements(panel_keys)
    for panel_id=0,nsubpanel-2 do begin
        tpos = subpanel_poss[*,panel_id]
        plots, tpos[[0,2]], tpos[1]+[0,0], /normal, linestyle=1
    endfor

    ; Ytitle.
    tpos = panel_pos
    tx = tpos[0]-xchsz*4
    ty = total(tpos[[1,3]])*0.5-ychsz*half_ychsz
    xyouts, tx,ty,/normal, alignment=1, lim.ytitle, orientation=90

    ; Labels.
    labels = lim.labels
    colors = lim.colors
    nlabel = n_elements(labels)
    panel_ysize = total(tpos[[1,3]]*[-1,1])
    dy = panel_ysize/(nlabel+1)
    label_ys = tpos[1]+dy+findgen(nlabel)*dy-ychsz*panel_label_size*half_ychsz
    for ii=0, nlabel-1 do begin
        tx = tpos[2]+xchsz*1
        ty = label_ys[ii]
        xyouts, tx,ty,/normal, labels[ii], color=colors[ii], charsize=panel_label_size
    endfor


    sgclose

end
