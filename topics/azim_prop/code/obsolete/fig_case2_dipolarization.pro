;+
; A plot for dipolarization progagpation for case 2.
;-

;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 1

;---Settings.
    einfo = get_var_data(info_var)
    time_range = einfo.time_range_plot1
    scs = einfo.scs
    sc_symbols = einfo.sc_symbols
    sc_colors = sgcolor(['red','orange','gold','lime','deep_sky_blue','medium_blue'])
    nsc = n_elements(scs)
    pres = scs+'_'

    figfn = sparentdir(srootdir())+'/plot/fig_case2_dipolarization.pdf'
    if keyword_set(test) then figfn = 0

;---Constant.
    deg = 180d/!dpi
    rad = !dpi/180d


;---Get the time and location of the dipolarization.
    df_var = 'df_info'
    load = tnames(df_var) eq ''
    if load then begin
        df_infos = replicate({$
            sc:'', $
            time:0d, $
            mlon:0., $
            mlt:0.}, nsc)
        df_times = time_double([$
            '2014-08-28/10:06:30', $
            '2014-08-28/10:09:05', $
            '2014-08-28/10:11:05', $
            '2014-08-28/10:19:43', $
            '2014-08-28/10:43:13', $
            '2014-08-28/10:44:15'])

        foreach sc, scs, ii do begin
            df_infos[ii].sc = sc
            df_infos[ii].time = df_times[ii]
            df_infos[ii].mlon = get_var_data(sc+'_mlon', at=df_times[ii])
            df_infos[ii].mlt = mlon2mlt(df_infos[ii].mlon, df_times[ii])
        endforeach
        store_data, df_var, transpose(df_infos.time), transpose(df_infos.mlon), limits={$
            ytitle:'(deg)', psym:sc_symbols[0], colors:sc_colors}
    endif


;---Plot.
    vars = [pres+'dbmag','df_info']
    nvar = n_elements(vars)

    sgopen, figfn, xsize=5, ysize=6, /inch

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    poss = sgcalcpos(nvar, lmargin=9, rmargin=8, tmargin=2, bmargin=8)
    poss[1,nvar-1] = 4*ychsz    ; make the last panel bigger.
    
    ; Plot settings.
    yticklen = -0.01
    options, vars, 'yticklen', yticklen
    
    tvar = 'thd_dbmag'
    options, tvar, 'yrange', [-20,30]
    options, tvar, 'ytickv', [-20,5,30]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'labels', 'THD'
    
    tvar = 'g15_dbmag'
    options, tvar, 'yrange', [-20,30]
    options, tvar, 'ytickv', [-20,5,30]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'labels', 'G15'
    
    tvar = 'the_dbmag'
    options, tvar, 'yrange', [0,40]
    options, tvar, 'ytickv', [0,20,40]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'labels', 'THE'
    
    tvar = 'tha_dbmag'
    options, tvar, 'yrange', [0,30]
    options, tvar, 'ytickv', [0,15,30]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'labels', 'THA'
    
    tvar = 'rbb_dbmag'
    options, tvar, 'yrange', [-10,10]
    options, tvar, 'ytickv', [-10,0,10]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'labels', 'RBB'
    
    tvar = 'g13_dbmag'
    options, tvar, 'yrange', [-12,4]
    options, tvar, 'ytickv', [-12,-4,4]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 4
    options, tvar, 'labels', 'G13'
    
    tvar = df_var
    yrange = [-180,-90]+100
    ytickv = [-180,-150,-120,-90]+100
    yticks = n_elements(ytickv)-1
    yminor = 3
    options, tvar, 'yrange', yrange
    options, tvar, 'ytickv', ytickv
    options, tvar, 'yticks', yticks
    options, tvar, 'yminor', yminor
    options, tvar, 'symsize', 0.5
    options, tvar, 'colors', sc_colors
    options, tvar, 'ytitle', 'MLon (deg)'
    
    thick = (size(ofn,/type) eq 7)? 4:2
    

    ; Plot vars.
    tplot, vars, position=poss, /novtitle, figlab=figlab, trange=time_range

    ; Add labels.
    fig_labels = ['a','b','c','d','e','f','g','h']+'.'
    for ii=0, nvar-1 do begin
        tpos = poss[*,ii]
        xyouts, xchsz*2, tpos[3]-ychsz*0.8, /normal, fig_labels[ii]
    endfor
    
    for ii=0, nsc-1 do begin
        tpos = poss[*,ii]
        get_data, scs[ii]+'_dbmag', limits=lim
        plot, time_range, lim.yrange, position=tpos, xstyle=5, ystyle=5, /nodata, /noerase
        plots, df_infos[ii].time+[0,0], lim.yrange, color=sc_colors[ii]
    endfor
    
    tpos = poss[*,0]
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*1, /normal, '|B|-|B!Dmodel!N|'

    tpos = poss[*,nvar-1]
    get_data, df_var, limit=lim
    plot, time_range, lim.yrange, /noerase, /nodata, position=tpos, xstyle=5, ystyle=5
    res = linfit(df_infos.time, df_infos.mlon)
    oplot, time_range, time_range*res[1]+res[0], linestyle=1, color=sgcolor('black')
    
    ytickn = strarr(yticks+1)
    tut = mean(time_range)
    for ii=0, yticks do ytickn[ii] = string(mlon2mlt(ytickv[ii],tut),format='(F4.1)')
    axis, yaxis=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, yticklen=yticklen, ytitle='MLT (hr)'

    xyouts, tpos[0]+xchsz, tpos[3]-ychsz*1, /normal, alignment=0, $
        'DF eastward @ '+string(res[1],format='(F5.3)')+' deg/sec'
    sgclose


end