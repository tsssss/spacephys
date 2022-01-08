

;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data
    
    test = 0
    
    einfo = get_var_data(info_var)
    time_range = einfo.time_range_plot1
    scs = einfo.scs
    sc_symbols = einfo.sc_symbols
    sc_colors = sgcolor(['red','orange','gold','lime','deep_sky_blue','medium_blue'])
    nsc = n_elements(scs)
    pres = scs+'_'
    
    x_rng = [5,-15]
    y_rng = [5,-10]
    ticklen = -0.01
    
    pos1 = [0.15,0.2,0.9,0.95]
    
    ofn = sparentdir(srootdir())+'/plot/fig_case2_sc_pos.pdf'
    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=4, ysize=3, /inch

    plot, x_rng, y_rng, /nodata, /noerase, position=pos1, $
        xstyle=1, xticks=4, xminor=5, xrange=x_rng, xticklen=ticklen, xtitle='GSM X (Re)', $
        ystyle=1, yticks=3, yminor=5, yrange=y_rng, yticklen=ticklen, ytitle='GSM Y (Re)'
    
    ; Add earth.
    tmp = findgen(21)/20*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    polyfill, xs<0, ys, /line_fill, orientation=45
    plots, xs, ys
    
    plots, x_rng, [0,0], linestyle=1
    plots, [0,0], y_rng, linestyle=1
    
    vars = scs+'_r_gsm'
    foreach tvar, vars, i do begin
        get_data, tvar, uts, rgsm
        idx = where(uts ge time_range[0] and uts le time_range[1])
        xs = rgsm[idx,0]
        ys = rgsm[idx,1]
        xs = interpol(rgsm[*,0], uts, mean(time_range))
        ys = interpol(rgsm[*,1], uts, mean(time_range))
        plots, xs, ys, color=sc_colors[i], psym=6
        xyouts, xs, ys, /data, alignment=0, '  '+strupcase(strmid(tvar,0,3)), color=sc_colors[i]
    endforeach
        
    sgclose
end
