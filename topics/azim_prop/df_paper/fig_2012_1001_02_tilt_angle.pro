
pro fig_2012_1001_02_tilt_angle, sort_by=sort_by, sort=sort

    test = 0

;---Check inputs.
    if ~keyword_set(sort_by) then sort_by = 'mlt'


;---Load data.
    store_data, '*', /delete
    id = '2012_1001_02'
    df_paper_load_data, event=event, id=id
    event_info = event[id]
    time_range = event_info.time_range
    probes = event_info.probes
    nprobe = n_elements(probes)

    
;---Add settings.
    df_info = dictionary()
    foreach probe, probes do df_info[probe] = dictionary()
    df_info['g15'].time_range = time_double(['2012-10-01/02:52:30', '2012-10-01/03:14:39'])
    df_info['tha'].time_range = time_double(['2012-10-01/02:41:39', '2012-10-01/02:49:00'])
    df_info['the'].time_range = time_double(['2012-10-01/02:40:00', '2012-10-01/02:47:09'])
    df_info['thd'].time_range = time_double(['2012-10-01/02:36:51', '2012-10-01/02:46:09'])
    df_info['g13'].time_range = time_double(['2012-10-01/02:23:09', '2012-10-01/02:29:50'])
    foreach probe, probes do begin
        df_info[probe].time = mean(df_info[probe].time_range)
        df_info[probe].rgsm = get_var_data(probe+'_r_gsm', at=df_info[probe].time)
        df_info[probe].rsm = get_var_data(probe+'_r_sm', at=df_info[probe].time)
        df_info[probe].mlt = get_var_data(probe+'_mlt', at=df_info[probe].time)
        df_info[probe].mlon = get_var_data(probe+'_mlon', at=df_info[probe].time)
        df_info[probe].mlat = get_var_data(probe+'_mlat', at=df_info[probe].time)
    endforeach


;---Tune ranges.
    set_yaxis, 'g15_b_tilt', range=[10,90]
    set_yaxis, 'tha_b_tilt', range=[-30,90]
    set_yaxis, 'the_b_tilt', range=[-10,90]
    set_yaxis, 'thd_b_tilt', range=[10,90]
    set_yaxis, 'g13_b_tilt', range=[0,40]


;---Order the variables.
    sort_data = fltarr(nprobe)
    foreach probe, probes, ii do sort_data[ii] = (df_info[probe])[sort_by]
    sorted_probes = probes[sort(sort_data)]
    if n_elements(sort) eq nprobe then sorted_probes = sort


;---Plot settings.
    fig_xsize = 4   ; inch.
    fig_ysize = 8   ; inch.
    lmargin = 8
    rmargin = 6
    tmargin = 2
    bmargin = 5
    ticklen = -0.01

    color_start = 50
    color_end = 250
    color_table = 40
    colors = round(smkarthm(color_start,color_end, nprobe, 'n'))
    for ii=0, nprobe-1 do colors[ii] = sgcolor(colors[ii], ct=color_table)


    if keyword_set(test) then begin
        file = test
        magnify = 1.5
    endif else begin
        file = join_path([srootdir(),'fig_2012_1001_02_tilt_angle.pdf'])
        magnify = 1
    endelse

    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify

    ; Plot the variables.
    poss = sgcalcpos(nprobe+1, ypans=[intarr(nprobe)+1,1.5], $
        lmargin=lmargin,rmargin=rmargin,tmargin=tmargin,bmargin=bmargin, xchsz=xchsz, ychsz=ychsz)
    full_ysz = 0.8
    half_ysz = 0.3
    label_size = 0.8
    lineskip = 0.3
    psym = 8
    tmp = findgen(21)/20*2*!dpi
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs,tys, /fill
    hsize = (size(file,/type) eq 7)? 80: 8
    hr2deg = 15.

    xrange = time_range
    xstep = 20*60   ; sec.
    xminor = 4
    xtickv = smkarthm(xrange[0],xrange[1], xstep, 'dx')
    xticks = n_elements(xtickv)-1
    xtickn = time_string(xtickv,tformat='hhmm') & xtickn[0] = time_string(xtickv[0],tformat='YYYY-MM-DD          ')
    plot_vars = sorted_probes+'_b_tilt'
    panel_labels = ['a','b','c','d','e','f','g','h']
    for ii=0, nprobe-1 do begin
        get_data, plot_vars[ii], txs, tys, limits=lim
        tpos = poss[*,ii]
        probe = sorted_probes[ii]

        xtickformat = '(A1)'
        xticklen = (tpos[2]-tpos[0])/(tpos[3]-tpos[1])*fig_xsize/fig_ysize*ticklen
        yticklen = ticklen

        yrange = lim.yrange

        plot, xrange, yrange, /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
        oplot, txs, tys
        plot, xrange, yrange, /nodata, /noerase, position=tpos, $
            xstyle=1, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickn=xtickn, xtickformat=xtickformat, xticklen=xticklen, $
            ystyle=1, yticklen=yticklen, _extra=lim

        tx = tpos[2]+xchsz*1
        ty = mean(tpos[[1,3]])-ychsz*half_ysz
        mission = resolve_probe(probe)
        label = strupcase(mission.short_name+mission.probe)
        xyouts, tx,ty,/normal, label, color=colors[ii]

        if ii eq 0 then begin
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*full_ysz
            xyouts, tx,ty,/normal, 'B SM tilt angle', charsize=label_size
        endif

        ; Add arrow.
        ttxs = df_info[probe].time_range
        ttys = interpol(tys,txs, ttxs)
        for jj=0,1 do begin
            tmp = convert_coord(ttxs[jj],ttys[jj], /data,/to_normal)
            ttxs[jj] = tmp[0]
            ttys[jj] = tmp[1]
        endfor
        tx = mean(ttxs)
        ty = mean(ttys)-ychsz*full_ysz
        plots, ttxs,ty+[0,0],/normal, color=colors[ii]
        for jj=0,1 do plots, ttxs[jj]+[0,0],ty+[-1,1]*ychsz*0.5*half_ysz,/normal, color=colors[ii]
        plots, tx,ty,/normal, psym=psym, symsize=label_size*0.5, color=colors[ii]
    endfor


    tpos = poss[*,nprobe]

    yrange = [-6,0]
    yticks = 3
    ytickv = smkarthm(yrange[0],yrange[1],yticks+1,'n')
    yminor = 5
    ytitle = 'MLT (hr)'
    ytickn = strarr(yticks+1)
    for ii=0, yticks do ytickn[ii] = sgnum2str(ytickv[ii])

    xticklen = (tpos[2]-tpos[0])/(tpos[3]-tpos[1])*fig_xsize/fig_ysize*ticklen

    plot, xrange, yrange, /nodata, /noerase, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickn=xtickn, xticklen=xticklen, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickn=ytickn, yticklen=yticklen, ytitle=ytitle
    fit_xs = fltarr(nprobe)
    fit_ys = fltarr(nprobe)
    foreach probe, sorted_probes, ii do begin
        ttxs = df_info[probe].time_range
        ttys = df_info[probe].mlt+[0,0]

        fit_xs[ii] = mean(ttxs)
        fit_ys[ii] = mean(ttys)

        for jj=0,1 do begin
            tmp = convert_coord(ttxs[jj],ttys[jj], /data,/to_normal)
            ttxs[jj] = tmp[0]
            ttys[jj] = tmp[1]
        endfor
        tx = mean(ttxs)
        ty = mean(ttys)
        plots, ttxs,ty+[0,0],/normal, color=colors[ii]
        for jj=0,1 do plots, ttxs[jj]+[0,0],ty+[-1,1]*ychsz*0.5*half_ysz,/normal, color=colors[ii]
        plots, tx,ty,/normal, psym=psym, symsize=label_size*0.5, color=colors[ii]
    endforeach

    nline = 0
    foreach flag, ['gt','lt'] do begin
        index = lazy_where(fit_ys, flag, 0, count=count)
        if count eq 0 then continue
        nline += 1

        fit_coef = linfit(fit_xs[index]-xrange[0],fit_ys[index])
        ttxs = xrange
        ttys = (ttxs-xrange[0])*fit_coef[1]+fit_coef[0]
        oplot, ttxs,ttys, linestyle=1

        direction = (fit_coef[1] ge 0)? 'eastward': 'westward'
        label = 'DF '+direction+' @'+sgnum2str(abs(fit_coef[1])*hr2deg,nsgn=2)+' deg/sec'
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*full_ysz*nline
        xyouts, tx,ty,/normal, label, charsize=label_size
    endforeach



    ; Add label.
    tx = xchsz*3
    for ii=0, nprobe do begin
        tpos = poss[*,ii]
        ty = tpos[3]-ychsz*full_ysz
        xyouts, tx,ty,/normal, panel_labels[ii]+'.'
    endfor

    if keyword_set(test) then stop
    sgclose

;---Panel for position.
    lmargin = 8
    rmargin = 6
    tmargin = 2
    bmargin = 5

    psym = 6

    sgopen, 0, xsize=fig_xsize, ysize=fig_xsize
    tpos = sgcalcpos(1,lmargin=lmargin,rmargin=rmargin,tmargin=tmargin,bmargin=bmargin, xchsz=xchsz, ychsz=ychsz)
    sgclose, /wdelete

    pos_panel = dictionary()
    xsm = fltarr(nprobe)
    ysm = fltarr(nprobe)
    foreach probe, sorted_probes, ii do begin
        rsm = df_info[probe].rsm
        xsm[ii] = rsm[0]
        ysm[ii] = rsm[1]
    endforeach
    xrange = minmax([xsm,-1,1])
    pos_panel.xrange = [ceil(xrange[1]),floor(xrange[0])]
    pos_panel.xrange = [1,-11]
    pos_panel.xtitle = 'SM X (Re)'
    pos_panel.xminor = 5
    pos_panel.xtickv = make_bins(pos_panel.xrange, pos_panel.xminor, /inner)
    pos_panel.xticks = n_elements(pos_panel.xtickv)-1
    yrange = minmax([ysm,-1,1])
    pos_panel.yrange = [ceil(yrange[1]),floor(yrange[0])]
    pos_panel.ytitle = 'SM Y (Re)'
    pos_panel.yminor = 5
    pos_panel.ytickv = make_bins(pos_panel.yrange, pos_panel.yminor, /inner)
    pos_panel.yticks = n_elements(pos_panel.ytickv)-1
    pos_panel.aspect_ratio = abs(double(pos_panel.yrange[1]-pos_panel.yrange[0])/(pos_panel.xrange[1]-pos_panel.xrange[0]))
    pos_panel.xticklen = ticklen/pos_panel.aspect_ratio
    pos_panel.yticklen = ticklen

    pos_panel.xsize = tpos[2]-tpos[0]
    pos_panel.ysize = pos_panel.xsize*pos_panel.aspect_ratio
    fig_ysize = (tmargin*ychsz+pos_panel.ysize+bmargin*ychsz)*fig_xsize

    if keyword_set(test) then begin
        file = test+1
        magnify = 2
    endif else begin
        file = join_path([srootdir(),'fig_2012_1001_02_sc_pos.pdf'])
        magnify = 1
    endelse

    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify

    tpos = sgcalcpos(1,lmargin=lmargin,rmargin=rmargin,tmargin=tmargin,bmargin=bmargin)
    plot, pos_panel.xrange, pos_panel.yrange, $
        xstyle=5, ystyle=5, $
        /nodata, /noerase, position=tpos, /isotropic, _extra=pos_panel.tostruct()
    foreach probe, sorted_probes, ii do begin
        plots, xsm[ii], ysm[ii], color=colors[ii], psym=psym, symsize=label_size

        ; Add spacecraft short name.
        mission = resolve_probe(probe)
        label = strupcase(mission.short_name+mission.probe)

        tmp = convert_coord(xsm[ii], ysm[ii], /data, /to_normal)
        tx = tmp[0]+xchsz*1*label_size
        ty = tmp[1]
        xyouts, tx,ty,/normal, label, color=colors[ii], charsize=label_size
    endforeach

    ; Add earth and lines.
    tmp = 50
    tmp = findgen(tmp)/(tmp-1)*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    polyfill, xs<0, ys, /line_fill, orientation=45
    plots, xs, ys
    foreach r, [5,10] do oplot, xs*r, ys*r, linestyle=1

    plots, pos_panel.xrange, [0,0], linestyle=1
    plots, [0,0], pos_panel.yrange, linestyle=1

    plot, pos_panel.xrange, pos_panel.yrange, $
        xstyle=1, ystyle=1, /isotropic, $
        /nodata, /noerase, position=tpos, _extra=pos_panel.tostruct()

    if keyword_set(test) then stop
    sgclose

end

fig_2012_1001_02_tilt_angle
end
