;+
; Plot DF properties including:
;   a. DF height.
;   b. DF width.
;-

test = 1
    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_search_df_group_filter(project=project)
    ;
    ; z = (MLT+0.0)/4.2
    ;events = azim_df_search_event(project=project, /get_df_events) ; All DFs.
    ;
    ; z = (MLT+0.2)/3.8
    ;events = azim_df_search_event(project=project, /get_large_df_events) ; Removing small DFs using the scaling coef.




;---This part explores the error in timing, by comparing the arrival time and xcor time.
;    xcor_dts = list()
;    foreach event, events do foreach df, event.df_list do xcor_dts.add, df.xcor_time_lag
;    xcor_dts = xcor_dts.toarray()
;    xcor_yys = histogram(xcor_dts, binsize=10, locations=xcor_xxs)
;    index = where(xcor_xxs eq 0)
;    xcor_yys[index] = mean(xcor_yys[index+[-1,1]])
;
;    sgopen, 0, xsize=3, ysize=3, /inch, xchsz=xchsz, ychsz=ychsz
;    tpos = sgcalcpos(1)
;    plot, xcor_xxs, xcor_yys, psym=1, symsize=0.5, $
;        xstyle=1, xrange=[-1,1]*600, xtitle='dT (sec)', xticks=4, xminor=3, xticklen=-0.01, $
;        ystyle=1, yrange=[0,60], ytitle='Count (#)', yticks=3, yminor=4, yticklen=-0.01, position=tpos
;    yys = gaussfit(xcor_xxs, smooth(xcor_yys,3), nterm=3, coef)
;    oplot, xcor_xxs, yys, color=sgcolor('red')
;    tx = tpos[0]+xchsz*0.5
;    ty = tpos[3]-ychsz*1
;    hmfw = 2*sqrt(2*alog(2))*coef[2]
;    xyouts, tx,ty,/normal, 'FWHM='+sgnum2str(hmfw,ndec=0)+' sec'
;    ;fit_result = linfit(xcor_xxs^2, alog(xcor_yys))
;    ;a = exp(fit_result[0])
;    ;b = fit_result[1]
;    ;oplot, xcor_xxs, a*exp(b*xcor_xxs^2), color=sgcolor('red')



;;---This part calculates the spatial width in each event.
;    df_spatial_widths = list()
;    the_events = azim_df_load_event(project=project, /after_analysis)
;    rad = constant('rad')
;    foreach event, the_events do begin
;        the_temporal_widths = list()
;        the_spatial_widths = list()
;        the_rxys = list()
;        the_omega = abs(event.omega_linear_fit)
;
;        foreach df, event.df_list do begin
;            rxy = snorm(df.arrival_r_sm[0:1])
;            width = df.width/60
;            the_temporal_widths.add, width
;            the_rxys.add, rxy
;            the_spatial_widths.add, width*the_omega*rad*rxy
;        endforeach
;
;        the_temporal_widths = the_temporal_widths.toarray()
;        the_spatial_widths = the_spatial_widths.toarray()
;        the_rxys = the_rxys.toarray()
;        df_spatial_widths.add, mean(the_spatial_widths)
;    endforeach
;    df_spatial_widths = df_spatial_widths.toarray()
;    print, mean(df_spatial_widths)
;    stop


;---This part plots the width and height of large DFs in terms of MLT and Rxy.
    df_widths = list()
    df_heights = list()
    df_rxys = list()
    df_mlts = list()
    foreach event, events do foreach df, event.df_list do begin
        df_widths.add, df.width
        df_heights.add, df.height
        df_rxys.add, snorm(df.arrival_r_sm[0:1])
        df_mlts.add, azim_df_calc_pseudo_mlt(df.arrival_r_sm)
    endforeach
    df_widths = df_widths.toarray()
    df_heights = df_heights.toarray()
    df_rxys = df_rxys.toarray()
    df_mlts = df_mlts.toarray()


    file = join_path([project.plot_dir, 'fig_df_width_and_height.pdf'])
    if keyword_set(test) then file = 0
    fig_xsize = 4
    fig_ysize = 4
    magnify = keyword_set(test)? 2:1
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch, magnify=magnify

;---Plot settings.
    margins = [7,4,2,1]
    ypad = 1
    xpad = 1.5
    poss = sgcalcpos(2,2, margins=margins, xchsz=xchsz, ychsz=ychsz, ypad=ypad, xpad=xpad)

    mlt_range = [-1,1]*12
    mlt_step = 6
    mlt_minor = 3
    mlt_title = 'MLT (hr)'
    height_range = [1.,200]
    height_step = 30
    height_minor = 3
    height_title = tex2str('Delta')+tex2str('theta')+' (deg)'
    height_color = sgcolor('wheat')

    rxy_range = [3,33]
    rxy_step = 10
    rxy_minor = 5
    rxy_title = 'R!Dxy!N (Re)'

    width_range = [60.,3600]/60
    width_title = tex2str('Delta')+'T (min)'
    width_color = sgcolor('silver')

    xticklen = -0.02
    yticklen = -0.02
    label_size = constant('label_size')


;---Panel a. Height vs MLT.
    tpos = poss[*,0,0]
    xrange = mlt_range
    xstep = mlt_step
    xminor = mlt_minor
    xtickv = make_bins(xrange,xstep)
    xticks = n_elements(xtickv)-1
    xlog = 0
    xtitle = ''
    xtickformat = '(A1)'

    yrange = height_range
    ylog = 1
    ytitle = height_title
    ytickformat = ''

    bin_psym = 6
    bin_symsize = label_size
    symbol_color = height_color


    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, xlog=xlog, $
        ystyle=5, yrange=yrange, ylog=ylog, $
        /nodata, /noerase, position=tpos

    xxs = df_mlts
    yys = df_heights
    oplot, xxs, yys, psym=1, symsize=label_size, color=symbol_color
    ; Do binning.
    bins = make_bins(xrange,xstep/xminor)
    nbin = n_elements(bins)-1
    bin_xxs = (bins[1:nbin]+bins[0:nbin-1])*0.5
    bin_yys = fltarr(nbin)+!values.f_nan
    test_yys = fltarr(nbin)+!values.f_nan
    test_ratio = 0.7
    min_count = 10
    for ii=0, nbin-1 do begin
        index = where_pro(xxs, '[]', bins[ii:ii+1], count=count)
        if count le min_count then continue
        the_data = yys[index]
        bin_yys[ii] = median(the_data)
        sorted_data = the_data[sort(the_data)]
        test_yys[ii] = sorted_data[count*test_ratio]
    endfor
    plots, bin_xxs, bin_yys, psym=bin_psym, symsize=bin_symsize


    ; fit.
    index = where(finite(bin_yys))
    fit_xxs = bin_xxs[index]
    fit_yys = gaussfit(fit_xxs, bin_yys[index], coefs, nterm=3)
    oplot, fit_xxs, fit_yys, linestyle=fit_linestyle
    msg = 'Fit: '+tex2str('Delta')+tex2str('theta')+' = '+sgnum2str(coefs[0],ndec=1)+tex2str('times')+'exp(-z!U2!N/2),!Cz=(MLT'
    if coefs[1] ge 0 then msg += '-'+sgnum2str(coefs[1],ndec=1) else msg += '+'+sgnum2str(-coefs[1],ndec=1)
    msg += ')/'+sgnum2str(abs(coefs[2]),ndec=1)
    tx = tpos[2]-xchsz*1
    ty = tpos[3]-ychsz
    xyouts, tx,ty,/normal, msg, alignment=1, charsize=label_size

    label_color = sgcolor('salmon')
    txxs = bin_xxs
    tyys = 10*exp(-(txxs/3)^2/2)
    oplot, txxs, tyys, linestyle=1, color=label_color
    tx = (tpos[0]+tpos[2])*0.5
    ty = tpos[1]+ychsz*3
    msg = 'Lower!Cboundary!C10'+tex2str('times')+'exp(-z!U2!N/2),!Cz=MLT/3'
    xyouts, tx,ty,/normal, msg, alignment=0.5, charsize=label_size, color=label_color

    ; Plot axes.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xlog=xlog, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ylog=ylog, ytitle=ytitle, yticklen=yticklen, ytickformat=ytickformat, $
        /nodata, /noerase, position=tpos



;---Panel b.
    tpos = poss[*,0,1]

    yrange = width_range
    ylog = 1
    ytitle = width_title

    xtickformat = ''
    xtitle = mlt_title

    symbol_color = width_color


    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, xlog=xlog, $
        ystyle=5, yrange=yrange, ylog=ylog, $
        /nodata, /noerase, position=tpos

    xxs = df_mlts
    yys = df_widths/60
    oplot, xxs, yys, psym=1, symsize=label_size, color=symbol_color
    ; Do binning.
    bins = make_bins(xrange,xstep/xminor)
    nbin = n_elements(bins)-1
    bin_xxs = (bins[1:nbin]+bins[0:nbin-1])*0.5
    bin_yys = fltarr(nbin)+!values.f_nan
    min_count = 10
    for ii=0, nbin-1 do begin
        index = where_pro(xxs, '[]', bins[ii:ii+1], count=count)
        if count le min_count then continue
        the_data = yys[index]
        bin_yys[ii] = median(the_data)
    endfor
    plots, bin_xxs, bin_yys, psym=bin_psym, symsize=bin_symsize

    ; Plot axes.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xlog=xlog, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ylog=ylog, ytitle=ytitle, yticklen=yticklen, ytickformat=ytickformat, $
        /nodata, /noerase, position=tpos


;---Panel c.
    tpos = poss[*,1,0]

    yrange = height_range
    ylog = 1
    ytitle = ''
    ytickformat = '(A1)'

    xrange = rxy_range
    xstep = rxy_step
    xminor = rxy_minor
    xtickv = make_bins(xrange,xstep, /inner)
    xticks = n_elements(xtickv)-1
    xlog = 0
    xtickformat = '(A1)'
    xtitle = ''

    symbol_color = height_color

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, xlog=xlog, $
        ystyle=5, yrange=yrange, ylog=ylog, $
        /nodata, /noerase, position=tpos

    xxs = df_rxys
    yys = df_heights
    oplot, xxs, yys, psym=1, symsize=label_size, color=symbol_color
    ; Do binning.
    bins = make_bins(xrange,xstep/xminor)
    nbin = n_elements(bins)-1
    bin_xxs = (bins[1:nbin]+bins[0:nbin-1])*0.5
    bin_yys = fltarr(nbin)+!values.f_nan
    min_count = 10
    for ii=0, nbin-1 do begin
        index = where_pro(xxs, '[]', bins[ii:ii+1], count=count)
        if count le min_count then continue
        the_data = yys[index]
        bin_yys[ii] = median(the_data)
    endfor
    plots, bin_xxs, bin_yys, psym=bin_psym, symsize=bin_symsize

    ; fit.
    slope = mean(bin_yys/bin_xxs, /nan)
    index = where(finite(bin_yys))
    fit_xxs = bin_xxs[index]
    fit_yys = slope*fit_xxs
    plots, fit_xxs, fit_yys, linestyle=fit_linestyle
    msg = 'Fit: '+tex2str('Delta')+tex2str('theta')+' = '+sgnum2str(abs(slope),ndec=1)+tex2str('times')+'R!Dxy!N'
    tx = tpos[2]-xchsz*1
    ty = tpos[3]-ychsz
    xyouts, tx,ty,/normal, msg, alignment=1, charsize=label_size

    ; Plot axes.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xlog=xlog, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ylog=ylog, ytitle=ytitle, yticklen=yticklen, ytickformat=ytickformat, $
        /nodata, /noerase, position=tpos


;---Panel d.
    tpos = poss[*,1,1]

    yrange = width_range
    ylog = 1
    ytitle = ''
    ytickformat = '(A1)'

    xtitle = rxy_title
    xtickformat = ''

    symbol_color = width_color


    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, xlog=xlog, $
        ystyle=5, yrange=yrange, ylog=ylog, $
        /nodata, /noerase, position=tpos

    xxs = df_rxys
    yys = df_widths/60
    oplot, xxs, yys, psym=1, symsize=label_size, color=symbol_color
    ; Do binning.
    bins = make_bins(xrange,xstep/xminor)
    nbin = n_elements(bins)-1
    bin_xxs = (bins[1:nbin]+bins[0:nbin-1])*0.5
    bin_yys = fltarr(nbin)+!values.f_nan
    min_count = 10
    for ii=0, nbin-1 do begin
        index = where_pro(xxs, '[]', bins[ii:ii+1], count=count)
        if count le min_count then continue
        the_data = yys[index]
        bin_yys[ii] = median(the_data)
    endfor
    plots, bin_xxs, bin_yys, psym=bin_psym, symsize=bin_symsize

    ; Fit.
    slope = mean(bin_yys*bin_xxs, /nan)
    oplot, bin_xxs, slope/bin_xxs, linestyle=fit_linestyle
    msg = 'Fit: '+tex2str('Delta')+'T = '+sgnum2str(abs(slope),ndec=1)+'/R!Dxy!N'
    tx = tpos[2]-xchsz*1
    ty = tpos[3]-ychsz
    xyouts, tx,ty,/normal, msg, alignment=1, charsize=label_size

    ; Plot axes.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xlog=xlog, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ylog=ylog, ytitle=ytitle, yticklen=yticklen, ytickformat=ytickformat, $
        /nodata, /noerase, position=tpos


    npanel = 4
    labels = letters(npanel)+'.'
    for ii=0,npanel-1 do begin
        tmp = array_indices([2,2], ii, /dimensions)
        tpos = poss[*,tmp[0],tmp[1]]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,/normal, labels[ii]
    endfor

    if keyword_set(test) then stop
    sgclose


end
