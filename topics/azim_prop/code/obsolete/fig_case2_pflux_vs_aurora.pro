;+
; A plot shows the mapping, correlation map, and pflux vs photon count.
;-

;+
; A figure shows the cross correlation of auroral brightness per pixel and Poynting flux.
;-
;

    test = 0
    all_models = 1
    add_grid = 1
    add_footpoint = 1
    model_chosen = 't89'


;---Get basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data
    einfo = get_var_data(info_var)
    time_range = time_double(['2013-06-07/04:54','2013-06-07/04:58'])


;---Settings.
    figfile = sparentdir(srootdir())+'/plot/fig_case2_pflux_vs_aurora.pdf'
    if keyword_set(test) then figfile = test

    probe = 'thd'
    probe_color = sgcolor(['magenta','red'])
    models = einfo.models
    model_symbols = [7,4,1,6]
    model_index = where(models eq model_chosen)
    if keyword_set(all_models) then model_index = indgen(n_elements(models))


    ; For drawing the snapshot.
    mlat_range = [61d,71]
    mlon_range = [-100d,-55]
    xrange = mlon_range
    yrange = mlat_range
    xminor = 5
    xtickv = smkarthm(xrange[0],xrange[1],xminor,'dx')
    xticks = n_elements(xtickv)-1
    xtickn = string(xtickv,format='(I0)') & xtickn[0] = ' '; & xtickn[xticks] = ' '
    yminor = 2
    ytickv = smkarthm(yrange[0],yrange[1],yminor,'dx')
    yticks = n_elements(ytickv)-1
    ytickn = string(ytickv,format='(I0)'); & ytickn[0] = ' ' & ytickn[yticks] = ' '
    xtitle = 'MLon (deg)'
    ytitle = 'MLat (deg)'
    ticklen = -0.01
    ct = 70
    nlevel = 21
    ztitle = 'Correlation'
    zrange = [-1,1]
    levels = smkarthm(zrange[0], zrange[1], nlevel,'n')
    colors = reverse(smkarthm(8,248,nlevel,'n'))

;---Crop the MLon images.
    mlonimg_var = pre0+'mlonimg'
    get_data, mlonimg_var, times, mlonimgs
    index = where(times ge time_range[0] and times le time_range[1])
    times = times[index]
    mlonimgs = mlonimgs[index,*,*]
    mlon_bins = get_var_data(pre0+'new_mlon_bins')
    mlat_bins = get_var_data(pre0+'new_mlat_bins')
    elev_2d = get_var_data(pre0+'new_elevs')
    edge_index = where(elev_2d le edge_elev)

    nmlon_bin = n_elements(mlon_bins)
    nmlat_bin = n_elements(mlat_bins)
    ntime = n_elements(times)
    mlon_index = where(mlon_bins ge mlon_range[0] and mlon_bins le mlon_range[1])
    mlat_index = where(mlat_bins ge mlat_range[0] and mlat_bins le mlat_range[1])
    cropped_images = mlonimgs[*,mlon_index,mlat_index]


;---Correlation map.
    dmlon = median(mlon_bins-shift(mlon_bins,1))
    dmlat = median(mlat_bins-shift(mlat_bins,1))
    mlon_err = 2d
    mlat_err = 1d
    dmlon_rec = mlon_err/dmlon*0.5
    dmlat_rec = mlat_err/dmlat*0.5

    pre1 = 'rb'+probe+'_'
    var = pre1+'pf_fac'
    if tnames(var) eq '' then _2013_0607_load_data2
    get_data, var, uts, pffac
    pfpara = interpol(pffac[*,0], uts, times)

    var = pre1+'corr_map'
    if test eq 1 then store_data, var, /delete
    if tnames(var) eq '' then begin
        fig_size = (size(cropped_images,/dimensions))[1:*]
        correlation_map = fltarr(fig_size)
        for ii=dmlon_rec, fig_size[0]-1-dmlon_rec do begin
            for jj=dmlat_rec, fig_size[1]-1-dmlat_rec do begin
                photon_count = fltarr(ntime)
                for kk=0, ntime-1 do photon_count[kk] = $
                    max(cropped_images[kk,ii-dmlon_rec:ii+dmlon_rec,jj-dmlat_rec:jj+dmlat_rec])
                photon_count = photon_count-(photon_count[0]+(photon_count[-1]-photon_count[0])*findgen(ntime)/(ntime-1))

                ;max_time_shift = 60
                ;lags = findgen(max_time_shift/3)    ; the lag for x, i.e., pflux.
                ;corrs = c_correlate(pfpara, photon_count, lags)
                ;max_corr = max(corrs, index)
                ;correlation_map[ii,jj] = corrs[index]

                correlation_map[ii,jj] = c_correlate(pfpara, photon_count, 0)
            endfor
        endfor

        store_data, var, 0, correlation_map
    endif


;---Plot correlation map.
    sgopen, figfile, xsize=4.5, ysize=5

    figpos = [0.15,0.10,0.85,0.95]
    poss = sgcalcpos(2, position=figpos)

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    pre1 = 'th'+probe+'_'

        get_data, pre1+'fpt_mlon', uts, fptmlon
        index = where(uts ge time_range[0] and uts le time_range[1], nrec)
        fptmlon = reform(total(fptmlon[index,*],1)/nrec)
        get_data, pre1+'fpt_mlat', uts, fptmlat
        index = where(uts ge time_range[0] and uts le time_range[1], nrec)
        fptmlat = reform(total(fptmlat[index,*],1)/nrec)


        tpos = poss[*,ii]

        get_data, pre1+'corr_map', tmp, zs
        xs = mlon_bins[mlon_index]
        ys = mlat_bins[mlat_index]

        if ii eq 0 then xtickformat= '(A1)' else xtickformat=''

        device, decomposed=0
        loadct, ct
        contour, zs, xs, ys, /noerase, $
            xstyle=5, xrange=mlon_range, xticklen=ticklen, xtitle=xtitle, $
            ystyle=5, yrange=mlat_range, yticklen=ticklen, ytitle=ytitle, $
            position=tpos, c_colors=colors, levels=levels, /fill
        device, decomposed=1
        ; Add grid.
        if keyword_set(add_grid) then begin
            plot, xrange, yrange, /nodata, position=tpos, $
                /noerase, ticklen=ticklen, $
                xstyle=1, xrange=xrange, xticks=xticks, xminor=0, xtickv=xtickv, xtitle='', xgridstyle=1, xtickformat='(A1)', xticklen=1, $
                ystyle=1, yrange=yrange, yticks=yticks, yminor=0, ytickv=ytickv, ytitle='', ygridstyle=1, ytickformat='(A1)', yticklen=1, $
                color=sgcolor('white')
        endif
        ; Add footpoint.
        foreach index, model_index do plots, fptmlon[index], fptmlat[index], psym=model_symbols[index], symsize=0.8, color=probe_colors[ii]
        ; Add ticks and axes.
        plot, xs, ys, /noerase, /nodata, $
            xstyle=1, xrange=mlon_range, xticklen=ticklen, xtitle='', xminor=xminor, xtickv=xtickv, xtickname=xtickn, xticks=xticks, xtickformat=xtickformat, $
            ystyle=1, yrange=mlat_range, yticklen=ticklen, ytitle='', yminor=yminor, ytickv=ytickv, ytickname=ytickn, yticks=yticks, $
            position=tpos

        if ii eq 1 then begin
            ; Add title.
            xyouts, tpos[0]-xchsz*6, tpos[1]-ychsz*1.6, /normal, xtitle, alignment=0
            xyouts, tpos[0]-xchsz*5, (figpos[1]+figpos[3])*0.5, /normal, ytitle, alignment=0.5, orientation=90
            ; Add label.
            foreach probe, probes, jj do xyouts, figpos[0]+xchsz*8*jj, figpos[3]+ychsz*0.5, /normal, color=probe_colors[jj], 'RBSP-'+strupcase(probe)
            foreach index, model_index, ii do begin
                tx = figpos[0]+xchsz*20+xchsz*ii*6
                ty = figpos[3]+ychsz*0.5
                xyouts, tx, ty, /normal, color=sgcolor('black'), strupcase(models[index])
                plots, tx-xchsz*1, ty+ychsz*0.30, /normal, psym=model_symbols[index], symsize=0.8
            endforeach
        endif
    endforeach



    cbpos = figpos
    cbpos[0] = tpos[2]+xchsz*1
    cbpos[2] = cbpos[0]+xchsz*0.8
    sgcolorbar, colors, zrange=zrange, ztitle=ztitle, position=cbpos, ct=ct
    if keyword_set(test) then stop
    sgclose
end
