

;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 0

    einfo = get_var_data(info_var)
    time_range = einfo.time_range_plot1

    scs = ['LANL-02A','LANL-04A','LANL-97A','1994-084','rbb','g13','1991-080','tha','the','g15','thd','LANL-01A']
    nsc = n_elements(scs)
    color_table1 = 40
    color_start = 50
    color_end = 250
    sc_colors = smkarthm(color_start,color_end, nsc, 'n')
    for ii=0, nsc-1 do sc_colors[ii] = sgcolor(sc_colors[ii],ct=color_table1)
    sc_colors[where(scs eq 'the')] = sgcolor('gold')    ; the yellow is too bright.
    

    xrange = [8,-13]
    xminor = 5
    xtickv = make_bins(xrange,xminor,/inner)
    xticks = n_elements(xtickv)-1
    yrange = [7,-8]
    yminor = 5
    ytickv = make_bins(yrange,yminor,/inner)
    yticks = n_elements(ytickv)-1
    ticklen = -0.01
    symsize = 0.5
    label_size = 0.8

    pos1 = [0.15,0.2,0.9,0.95]

    ofn = sparentdir(srootdir())+'/plot/fig_case2_sc_pos.pdf'
    if keyword_set(test) then begin
        ofn = test
        magnify = 2
    endif else begin
        magnify = 1
    endelse
    sgopen, ofn, xsize=4, ysize=3, magnify=magnify
    tpos = sgcalcpos(1, position=pos1, xchsz=xchsz, ychsz=ychsz)
    
    plot, xrange, yrange, /nodata, /noerase, position=tpos, /iso, $
        xstyle=1, xticks=xticks, xtickv=xtickv, xminor=xminor, xrange=xrange, xticklen=ticklen, xtitle='GSM X (Re)', $
        ystyle=1, yticks=yticks, ytickv=ytickv, yminor=yminor, yrange=yrange, yticklen=ticklen, ytitle='GSM Y (Re)'

    ; Add earth.
    tmp = 50
    tmp = findgen(tmp)/(tmp-1)*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    polyfill, xs<0, ys, /line_fill, orientation=45
    plots, xs, ys
    r = 5
    tx = xs*r
    ty = ys*r
    oplot, tx, ty, linestyle=1
    tmp = 2
    xyouts, tx[tmp], ty[tmp], '  '+string(r,format='(I0)')+' Re', charsize=label_size
    
    r = 10
    tx = xs*r
    ty = ys*r
    oplot, tx, ty, linestyle=1
    tmp = 21
    xyouts, tx[tmp], ty[tmp], '  '+string(r,format='(I0)')+' Re', charsize=label_size

    plots, xrange, [0,0], linestyle=1
    plots, [0,0], yrange, linestyle=1

    vars = scs+'_r_gsm'
    foreach tvar, vars, i do begin
        get_data, tvar, uts, rgsm
        idx = where(uts ge time_range[0] and uts le time_range[1])
        xs = rgsm[idx,0]
        ys = rgsm[idx,1]
        xs = interpol(rgsm[*,0], uts, mean(time_range))
        ys = interpol(rgsm[*,1], uts, mean(time_range))
        tmp = convert_coord(xs,ys, /data, /to_normal)
        tx = tmp[0]
        ty = tmp[1]
        plots, tx, ty, /normal, color=sc_colors[i], psym=6, symsize=symsize
        sc = strmid(tvar,0,strpos(tvar,'_'))
        if sc eq '1994-084' then begin
            alignment = 0.5
            ty += ychsz*0.5
            txt = strupcase(sc)
        endif else begin
            alignment = 0
            txt = '  '+strupcase(sc)
        endelse
        if sc eq 'rbb' then begin
            alignment = 0.5
            ty -= ychsz*1
            txt = strupcase(sc)
        endif
        if sc eq 'g15' then begin
            alignment = 1
            txt = strupcase(sc)+'  '
        endif
        xyouts, tx, ty, /normal, alignment=alignment, txt, color=sc_colors[i]
    endforeach

    sgclose
end
