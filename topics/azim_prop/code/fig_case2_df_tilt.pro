;+
; A plot for tilt progagpation for case 2.
;-

;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 0

;---Settings.
    einfo = get_var_data(info_var)
    time_range = einfo.time_range_plot1

    ; All colors and spacecraft.
    all_scs = ['LANL-02A','LANL-04A','LANL-97A','1994-084','rbb','g13','1991-080','tha','the','g15','thd','LANL-01A']
    nsc = n_elements(all_scs)
    color_table1 = 40
    color_start = 50
    color_end = 250
    all_sc_colors = smkarthm(color_start,color_end, nsc, 'n')
    for ii=0, nsc-1 do all_sc_colors[ii] = sgcolor(all_sc_colors[ii],ct=color_table1)
    all_sc_colors[where(scs eq 'the')] = sgcolor('gold')    ; the yellow is too bright.

    scs = einfo.scs
    nsc = n_elements(scs)
    sc_symbols = einfo.sc_symbols
    sc_colors = lonarr(nsc)
    for ii=0, nsc-1 do sc_colors[ii] = all_sc_colors[where(all_scs eq scs[ii])]
    pres = scs+'_'

    figfn = sparentdir(srootdir())+'/plot/fig_case2_df_tilt.pdf'
    if keyword_set(test) then figfn = 0

;---Constant.
    deg = 180d/!dpi
    rad = !dpi/180d

    label_size = 0.7

;---Get the time and location of the tilt.
    df_var = 'df_info_tilt'
    load = tnames(df_var) eq ''
    load = 1
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
    df_times = time_double([$
        '2014-08-28/10:12:50', $
        '2014-08-28/10:16:50', $
        '2014-08-28/10:20:45', $
        '2014-08-28/10:23:10', $
        '2014-08-28/10:47:15', $
        '2014-08-28/10:50:25'])
    df_times = time_double([$
        '2014-08-28/10:12:50', $
        '2014-08-28/10:16:50', $
        '2014-08-28/10:20:45', $
        '2014-08-28/10:23:10', $
        '2014-08-28/10:50:00', $
        '2014-08-28/10:54:40'])

    foreach sc, scs, ii do begin
        df_infos[ii].sc = sc
        df_infos[ii].time = df_times[ii]
        df_infos[ii].mlon = get_var_data(sc+'_mlon', at=df_times[ii])
        df_infos[ii].mlt = mlon2mlt(df_infos[ii].mlon, df_times[ii])
    endforeach
    store_data, df_var, transpose(df_infos.time), transpose(df_infos.mlon), limits={$
        ytitle:'(deg)', psym:sc_symbols[0], colors:sc_colors}


;---Get the tilt angle.
    foreach sc, scs, ii do begin
        if tnames(sc+'_b_tilt') ne '' then continue
        get_data, sc+'_b_gsm', times, bgsm
        bsm = cotran(bgsm, times, 'gsm2sm')
        tilt = atan(bsm[*,2],sqrt(bsm[*,0]^2+bsm[*,1]^2))*deg
        store_data, sc+'_b_tilt', times, tilt, limits={ytitle:'(deg)'}
    endforeach


;---Plot.
    vars = [pres+'b_tilt','df_info_tilt']
    nvar = n_elements(vars)

    sgopen, figfn, xsize=4, ysize=6, /inch

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    poss = sgcalcpos(nvar, lmargin=9, rmargin=8, tmargin=2, bmargin=8)
    poss[1,nvar-1] = 4*ychsz    ; make the last panel bigger.

    ; Plot settings.
    yticklen = -0.01
    options, vars, 'yticklen', yticklen

    tvar = 'thd_b_tilt'
    options, tvar, 'yrange', [20,80]
    options, tvar, 'ytickv', [20,50,80]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'labels', 'THD'

    tvar = 'g15_b_tilt'
    options, tvar, 'yrange', [20,80]
    options, tvar, 'ytickv', [20,50,80]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'labels', 'G15'

    tvar = 'the_b_tilt'
    options, tvar, 'yrange', [10,60]
    options, tvar, 'ytickv', [10,35,60]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'labels', 'THE'

    tvar = 'tha_b_tilt'
    options, tvar, 'yrange', [0,50]
    options, tvar, 'ytickv', [0,25,50]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'labels', 'THA'

    tvar = 'rbb_b_tilt'
    options, tvar, 'yrange', [18,26]
    options, tvar, 'ytickv', [18,22,26]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 4
    options, tvar, 'labels', 'RBB'

    tvar = 'g13_b_tilt'
    options, tvar, 'yrange', [33,49]
    options, tvar, 'ytickv', [33,41,49]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 4
    options, tvar, 'labels', 'G13'

    tvar = df_var
    yrange = [-180,-90]+100
    yrange = reverse(minmax(yrange))
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

    thick = (size(figfn,/type) eq 7)? 4:2
    arrow_size = (size(figfn,/type) eq 7)? 160:8


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
        tvar = scs[ii]+'_b_tilt'
        get_data, tvar, limits=lim
        plot, time_range, lim.yrange, position=tpos, xstyle=5, ystyle=5, /nodata, /noerase
        tx = df_infos[ii].time
        tmp = convert_coord(tx, get_var_data(tvar,at=tx), /data, /to_normal)
        tx = tmp[0]
        ty = tmp[1]-1*ychsz
        arrow, tx,ty-1*ychsz, tx,ty, /normal, /solid, hsize=arrow_size, color=sc_colors[ii]
    endfor

    tpos = poss[*,0]
    xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*1, /normal, 'B SM tilt angle', charsize=label_size

    tpos = poss[*,nvar-1]
    get_data, df_var, limit=lim
    yrange = reverse(minmax(lim.yrange))
    plot, time_range, yrange, /noerase, /nodata, position=tpos, xstyle=5, ystyle=5, $
        yrange=yrange

    res = linfit(df_infos.time, df_infos.mlon)
    oplot, time_range, time_range*res[1]+res[0], linestyle=1, color=sgcolor('black')

    ytickn = strarr(yticks+1)
    tut = mean(time_range)
    for ii=0, yticks do ytickn[ii] = string(mlon2mlt(ytickv[ii],tut),format='(F4.1)')
    axis, yaxis=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, yticklen=yticklen, ytitle='MLT (hr)'

    xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.5, /normal, alignment=0, $
        'DF eastward @ '+string(res[1],format='(F5.3)')+' deg/sec', charsize=label_size
    sgclose


end
