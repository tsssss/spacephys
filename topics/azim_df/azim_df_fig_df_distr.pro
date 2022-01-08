;+
; Load all DF, and plot the distribution of height vs MLT, Rxy, R.
;-

pro azim_df_fig_df_distr_load_data, project=project, varname=var

    events = azim_df_search_event(project=project, /get_df_events)

    df_count = 0.
    foreach event, events do df_count += event.df_count
    ndim = 3
    df_r_sms = fltarr(df_count,ndim)
    df_heights = fltarr(df_count)
    df_widths = fltarr(df_count)

    df_index = 0
    foreach event, events do begin
        foreach df_info, event.df_list do begin
            df_heights[df_index] = df_info.height
            df_widths[df_index] = df_info.width
            df_r_sms[df_index,*] = df_info.arrival_r_sm
            df_index += 1
        endforeach
    endforeach

    df_info = dictionary($
        'count', df_count, $
        'heights', df_heights, $
        'widths', df_widths, $
        'r_sms', df_r_sms)

    store_data, var, 0, df_info
end



pro azim_df_fig_df_distr, project=project

test = 0
    if n_elements(project) eq 0 then project = azim_df_load_project()
    magnify = keyword_set(test)? 1.5: 1
    plot_file = join_path([project.plot_dir,'fig_df_distr.pdf'])
    if keyword_set(test) then plot_file = test


;---Collect data from events to array.
    df_distr_var = 'df_distr_var'
    if check_if_update(df_distr_var) then $
        azim_df_fig_df_distr_load_data, project=project, varname=df_distr_var
    df_info = get_var_data(df_distr_var)
    df_count = df_info.count
    df_r_sms = df_info.r_sms
    df_heights = df_info.heights
    df_widths = df_info.widths

;---Derive more position quantities.
    df_mlts = azim_df_calc_pseudo_mlt(df_r_sms)
    df_rxys = snorm(df_r_sms[*,0:1])
    df_diss = snorm(df_r_sms)
    df_xsms = df_r_sms[*,0]

;---Filter by height.
    df_min_height = 1.
    index = where(df_heights ge df_min_height)
    df_mlts = df_mlts[index]
    df_rxys = df_rxys[index]
    df_diss = df_diss[index]
    df_heights = df_heights[index]
    df_widths = df_widths[index]
    df_xsms = df_xsms[index]

;---Bin data.
    mlt_range = project.search_candidate.search_roi.mlt_range
    mlt_bin_size = 0.5
    mlt_bin_boundarys = make_bins(mlt_range, mlt_bin_size)
    nmlt_bin = n_elements(mlt_bin_boundarys)-1
    mlt_bin_centers = mlt_bin_boundarys[0:nmlt_bin-1]+mlt_bin_size*0.5

    rxy_range = project.search_candidate.search_roi.rxy_range
    rxy_bin_size = 1
    rxy_bin_boundarys = make_bins(rxy_range, rxy_bin_size)
    nrxy_bin = n_elements(rxy_bin_boundarys)-1
    rxy_bin_centers = rxy_bin_boundarys[0:nrxy_bin-1]+rxy_bin_size*0.5

    height_range = [0,max(df_heights)]
    width_range = [0,max(df_widths)]
    str_delta_theta = tex2str('Delta')+tex2str('theta')
    str_delta_t = tex2str('Delta')+'T'

;---Figure settings.
    fig_xsize = 6.
    fig_ysize = 6.
    xticklen_chsz = -0.15   ; in ychsz.
    yticklen_chsz = -0.30   ; in xchsz.
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify
    nxpanel = 2
    nypanel = 2
    margins = [8,5,2,2]
    poss = sgcalcpos(nxpanel, nypanel, margins=margins, $
        ypad=5, xpad=8, xchsz=xchsz, ychsz=ychsz)

    height_color = sgcolor('silver')
    point_psym = 3  ; dot.
    point_psym = 1  ; cross.
    point_symsize = 0.2
    width_color = sgcolor('wheat')

    fit_color = sgcolor('black')
    fit_psym = 1
    fit_symsize = 0.3
    fit_linestyle = 2

    stddev_color = sgcolor('red')

    full_ychsz = constant('full_ychsz')
    half_ychsz = constant('half_ychsz')
    label_size = constant('label_size')

;---Panel a: height vs MLT.
    tpos = poss[*,0,0]
    xrange = mlt_range
    xstep = 3
    xtickv = make_bins(xrange, xstep)
    xrange = minmax(xtickv)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    xtitle = 'MLT (hr)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])

    yrange = height_range
    ystep = 30.
    ytickv = make_bins(yrange, ystep)
    yrange = minmax(ytickv)
    yticks = n_elements(ytickv)-1
    yminor = 3
    ytitle = tex2str('Delta')+tex2str('theta')+' (deg)'
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ; Bin data.
    height_mlt_medians = fltarr(nmlt_bin)
    height_mlt_means = fltarr(nmlt_bin)
    height_mlt_stddevs = fltarr(nmlt_bin)
    for ii=0, nmlt_bin-1 do begin
        index = lazy_where(df_mlts, '[]', mlt_bin_boundarys[ii:ii+1], count=count)
        if count le 10 then continue
        the_heights = df_heights[index]
        the_heights = the_heights[sort(the_heights)]
        the_heights = the_heights[0.2*count:0.8*count]
        height_mlt_means[ii] = mean(the_heights)
        height_mlt_medians[ii] = median(the_heights)
        height_mlt_stddevs[ii] = stddev(the_heights)
    endfor

    ; Create coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /noerase, /nodata

    ; Add data.
    xxs = df_mlts
    yys = df_heights
    plots, xxs, yys, color=height_color, psym=point_psym, symsize=point_symsize
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*1
    plots, tx+xchsz*1.1,ty+ychsz*half_ychsz*label_size,/normal, color=height_color, psym=point_psym, symsize=point_symsize
    xyouts, tx+xchsz*3,ty,/normal, tex2str('Delta')+tex2str('theta')+', change in tilt angle of dipolarizations', charsize=label_size


    ; Add stddev.
    height_mag = 10.
    xxs = mlt_bin_centers
    yys = height_mlt_stddevs*height_mag
    plots, xxs, yys, color=stddev_color, psym=fit_psym, symsize=fit_symsize
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*2
    plots, tx+xchsz*1.1,ty+ychsz*half_ychsz*label_size,/normal, color=stddev_color, psym=fit_psym, symsize=fit_symsize
    xyouts, tx+xchsz*3,ty,/normal, sgnum2str(height_mag)+tex2str('times')+' standard deviation of '+tex2str('Delta')+tex2str('theta')+' in each MLT bin', charsize=label_size


    ; Add fit.
    xxs = mlt_bin_centers
    yys = gaussfit(mlt_bin_centers, height_mlt_stddevs, fit_coeff, nterms=3)*height_mag
    plots, xxs, yys, color=fit_color, linestyle=fit_linestyle
    msg = 'Gauss fit: '+tex2str('Delta')+tex2str('theta')+' = '+$
        sgnum2str(fit_coeff[0]*height_mag,ndec=1)+' exp(-z!U2!N/2), z='+$
        '(MLT+'+sgnum2str(abs(fit_coeff[1]),ndec=1)+')/'+sgnum2str(fit_coeff[2],ndec=1)
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*3
    plots, tx+xchsz*[0,2.2],ty+ychsz*half_ychsz*label_size+[0,0],/normal, color=fit_color, linestyle=fit_linestyle
    xyouts, tx+xchsz*3,ty,/normal, msg, charsize=label_size

    ; Draw axis.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xtickv=xtickv, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, ytickv=ytickv, ytitle=ytitle, yticklen=yticklen, $
        position=tpos, /noerase, /nodata


;---Panel c: width vs MLT.
    tpos = poss[*,0,1]

    yrange = width_range/60
    ystep = 20
    ytickv = make_bins(yrange, ystep)
    yrange = minmax(ytickv)
    yticks = n_elements(ytickv)-1
    yminor = 5
    ytitle = str_delta_t+' (min)'
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ; Bin data.
    width_mlt_stddevs = fltarr(nmlt_bin)
    for ii=0, nmlt_bin-1 do begin
        index = lazy_where(df_mlts, '[]', mlt_bin_boundarys[ii:ii+1], count=count)
        if count le 10 then continue
        the_widths = df_widths[index]
        the_widths = the_widths[sort(the_widths)]
        the_widths = the_widths[0.2*count:0.8*count]
        width_mlt_stddevs[ii] = stddev(the_widths)
    endfor

    ; Create coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /noerase, /nodata

    ; Add data.
    xxs = df_mlts
    yys = df_widths/60
    plots, xxs, yys, color=width_color, psym=point_psym, symsize=point_symsize
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*1
    plots, tx+xchsz*1.1,ty+ychsz*half_ychsz*label_size,/normal, color=width_color, psym=point_psym, symsize=point_symsize
    xyouts, tx+xchsz*3,ty,/normal, tex2str('Delta')+'T, duration of dipolarizations', charsize=label_size

    ; Add stddev.
    width_mag = 10.
    xxs = mlt_bin_centers
    yys = width_mlt_stddevs/60*width_mag
    plots, xxs, yys, color=stddev_color, psym=fit_psym, symsize=fit_symsize
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*2
    plots, tx+xchsz*1.1,ty+ychsz*half_ychsz*label_size,/normal, color=stddev_color, psym=fit_psym, symsize=fit_symsize
    xyouts, tx+xchsz*3,ty,/normal, sgnum2str(height_mag)+tex2str('times')+' standard deviation of '+str_delta_t+' in each MLT bin', charsize=label_size


    ; Draw axis.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xtickv=xtickv, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, ytickv=ytickv, ytitle=ytitle, yticklen=yticklen, $
        position=tpos, /noerase, /nodata


;---Panel b: heights vs Rxy.
    tpos = poss[*,1,0]

    xrange = rxy_range
    xstep = 5.
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    xtitle = 'Rxy (Re)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])

    yrange = height_range
    ystep = 30.
    ytickv = make_bins(yrange, ystep)
    yrange = minmax(ytickv)
    yticks = n_elements(ytickv)-1
    yminor = 3
    ytitle = str_delta_theta+' (deg)'
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ; Bin data.
    height_rxy_stddevs = fltarr(nrxy_bin)
    for ii=0, nrxy_bin-1 do begin
        index = lazy_where(df_rxys, '[]', rxy_bin_boundarys[ii:ii+1], count=count)
        if count le 10 then continue
        the_heights = df_heights[index]
        the_heights = the_heights[sort(the_heights)]
        the_heights = the_heights[0.2*count:0.8*count]
        height_rxy_stddevs[ii] = stddev(the_heights)
    endfor

    ; Create coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /noerase, /nodata

    ; Add data.
    xxs = df_rxys
    yys = df_heights
    plots, xxs, yys, color=height_color, psym=point_psym, symsize=point_symsize
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*1
    plots, tx+xchsz*1.1,ty+ychsz*half_ychsz*label_size,/normal, color=height_color, psym=point_psym, symsize=point_symsize
    xyouts, tx+xchsz*3,ty,/normal, tex2str('Delta')+tex2str('theta')+', change in tilt angle of dipolarizations', charsize=label_size


    ; Add stddev.
    xxs = rxy_bin_centers
    yys = height_rxy_stddevs*height_mag
    plots, xxs, yys, color=stddev_color, psym=fit_psym, symsize=fit_symsize
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*2
    plots, tx+xchsz*1.1,ty+ychsz*half_ychsz*label_size,/normal, color=stddev_color, psym=fit_psym, symsize=fit_symsize
    xyouts, tx+xchsz*3,ty,/normal, sgnum2str(height_mag)+tex2str('times')+' standard deviation of '+str_delta_t+' in each Rxy bin', charsize=label_size

    ; Draw axis.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xtickv=xtickv, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, ytickv=ytickv, ytitle=ytitle, yticklen=yticklen, $
        position=tpos, /noerase, /nodata


;---Panel d: width vs Rxy.
    tpos = poss[*,1,1]

    xrange = rxy_range
    xstep = 5.
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    xtitle = 'Rxy (Re)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])

    yrange = width_range/60
    ystep = 20
    ytickv = make_bins(yrange, ystep)
    yrange = minmax(ytickv)
    yticks = n_elements(ytickv)-1
    yminor = 5
    ytitle = tex2str('Delta')+'T (min)'
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ; Bin data.
    width_rxy_stddevs = fltarr(nrxy_bin)
    for ii=0, nrxy_bin-1 do begin
        index = lazy_where(df_rxys, '[]', rxy_bin_boundarys[ii:ii+1], count=count)
        if count le 10 then continue
        the_widths = df_widths[index]
        the_widths = the_widths[sort(the_widths)]
        the_widths = the_widths[0.2*count:0.8*count]
        width_rxy_stddevs[ii] = stddev(the_widths)
    endfor

    ; Create coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /noerase, /nodata

    ; Add data.
    xxs = df_rxys
    yys = df_widths/60
    plots, xxs, yys, color=width_color, psym=point_psym, symsize=point_symsize
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*1
    plots, tx+xchsz*1.1,ty+ychsz*half_ychsz*label_size,/normal, color=width_color, psym=point_psym, symsize=point_symsize
    xyouts, tx+xchsz*3,ty,/normal, tex2str('Delta')+'T, duration of dipolarizations', charsize=label_size


    ; Add stddev.
    xxs = rxy_bin_centers
    yys = width_rxy_stddevs/60*width_mag
    plots, xxs, yys, color=stddev_color, psym=fit_psym, symsize=fit_symsize
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*2
    plots, tx+xchsz*1.1,ty+ychsz*half_ychsz*label_size,/normal, color=stddev_color, psym=fit_psym, symsize=fit_symsize
    xyouts, tx+xchsz*3,ty,/normal, sgnum2str(width_mag)+tex2str('times')+' standard deviation of '+tex2str('Delta')+'T in each Rxy bin', charsize=label_size


    ; Add fit.
;    xxs = rxy_bin_centers
;    fit_coeff = linfit(rxy_bin_centers, 1/(width_rxy_stddevs/60*width_mag), yfit=yys)
;    yys = 1/yys
;    plots, xxs, yys, color=fit_color, linestyle=fit_linestyle
;    msg = 'Linear fit: '+str_delta_t+' = '+sgnum2str(1/fit_coeff[1],ndec=1)+'/('+$
;        sgnum2str(fit_coeff[0]/fit_coeff[1],ndec=1)+'+Rxy)'
    fit_coeff = mean(rxy_bin_centers*width_rxy_stddevs/60*width_mag)
    xxs = rxy_bin_centers
    yys = fit_coeff[0]/rxy_bin_centers
    plots, xxs, yys, color=fit_color, linestyle=fit_linestyle
    msg = 'Linear fit: '+str_delta_t+' = '+sgnum2str(fit_coeff[0],ndec=1)+'/Rxy'
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*full_ychsz*3
    plots, tx+xchsz*[0,2.2],ty+ychsz*half_ychsz*label_size+[0,0],/normal, color=fit_color, linestyle=fit_linestyle
    xyouts, tx+xchsz*3,ty,/normal, msg, charsize=label_size



    ; Draw axis.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xtickv=xtickv, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, ytickv=ytickv, ytitle=ytitle, yticklen=yticklen, $
        position=tpos, /noerase, /nodata


    npanel = nxpanel*nypanel
    fig_labels = letters(npanel)+'.'
    for ii=0, npanel-1 do begin
        index = array_indices([nxpanel,nypanel], ii, /dimension)
        tpos = poss[*,index[0],index[1]]
        tx = tpos[0]-xchsz*5
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, fig_labels[ii]
    endfor

    if keyword_set(test) then stop
    sgclose

end
