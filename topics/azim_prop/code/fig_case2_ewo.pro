;+
; Generate EWOgram for case 2.
;-


;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 0

;---Settings.
    einfo = get_var_data(info_var)
    time_range_plot = time_double(['2014-08-28/10:10','2014-08-28/10:15'])
    count_range = [0,300]

    aurora_color_table = 49
    label_size = 0.7
    ticklen = -0.01


;---Calculate EWOgram.
    mlonimg_var = 'thg_mlonimg'
    ewo_var = mlonimg_var+'_ewo'

    ; Set settings.
    newo_info = 1
    ewo_infos = replicate({mlat_range:[0d,0], var_name:'', label:''},newo_info)
    ; the upper part, capture the expansion.
    ewo_infos[0].mlat_range = [66,69]
    ewo_infos[0].mlat_range = [67,68]+1
    ewo_infos[0].var_name = ewo_var+'_2'

    foreach tinfo, ewo_infos, ii do begin
        ;if tnames(tinfo.var_name) ne '' then continue
        themis_read_mlonimg_ewo, mlonimg_var, to=tinfo.var_name, mlat_range=tinfo.mlat_range
    endforeach

    ; The label.
    foreach tinfo, ewo_infos, ii do begin
        ewo_infos[ii].label = '['+sgnum2str(tinfo.mlat_range[0])+','+sgnum2str(tinfo.mlat_range[1])+']'
    endforeach


;---Prepare for plot.
    file = join_path([sparentdir(srootdir()),'plot','fig_case2_ewo.pdf'])
    magnify = 1
    if keyword_set(test) then begin
        file = 0
        magnify = 2
    endif
    sgopen, file, xsize=4, ysize=2.5, magnify=magnify

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    vars = ewo_infos.var_name
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, tmargin=4, bmargin=2.5, rmargin=3, lmargin=8)
    colorbar_pos = poss[[0,3,2,3],0]+[0,0.5,0,1]*ychsz

    options, vars, 'zrange', count_range
    options, vars, 'zposition', colorbar_pos
    options, vars, 'ztitle', 'Brightness (#)'
    options, vars, 'zhorizontal', 1
    options, vars, 'zticks', 3
    options, vars, 'zcharsize', label_size

    xrange = time_range_plot
    xminor = 4
    xtickv = smkarthm(xrange[0], xrange[1], 60, 'dx')
    xticks = n_elements(xtickv)-1
    xtickn = time_string(xtickv, tformat='hh:mm')
    xtickn[0] = 'UT (hh:mm)     '
    xticklen = ticklen*2

    yrange = [-60,-100]
    yticks = abs(yrange[1]-yrange[0])/20
    ytickv = smkarthm(yrange[0],yrange[1],yticks+1,'n')
    yminor = 4
    yticklen = ticklen
    ytitle = 'MLon (deg)'

    zrange = count_range
    ztitle = 'Brightness (#)'
    zticks = 3
    color_start = 0
    color_end = 254
    colors = round(smkarthm(color_start, color_end, 1, 'dx'))


    ; Plot.
    tpos = poss[*,0]
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        /nodata, /noerase, position=tpos
    get_data, ewo_infos[0].var_name, times, timg, mlon_bins
    index = lazy_where(mlon_bins, minmax(yrange))
    timg = timg[*,index]
    mlon_bins = mlon_bins[index]
    index = lazy_where(times, time_range_plot)
    timg = timg[index,*]
    timg = bytscl(timg, min=zrange[0], max=zrange[1], top=color_end)
    if yrange[1] lt yrange[0] then timg = reverse(timg,2)
    sgtv, timg, /resize, position=tpos, ct=aurora_color_table

    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickn=xtickn, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase, position=tpos


    ; The slope.
    get_data, ewo_infos[0].var_name, times, timg, mlon_bins
    time_range_slope = time_double(['2014-08-28/10:12:30','2014-08-28/10:13:10'])
    slope_brightness = 120
    slope_separator = -70.
    earth_rotation = 360d/86400
    slope_index = lazy_where(times, time_range_slope, count=npoint)
    for ii=0, 1 do begin
        txs = times[slope_index]
        tzs = timg[slope_index,*]
        tys = fltarr(npoint)+!values.f_nan
        if ii eq 0 then begin
            ; Eastward.
            tlabel = 'east'
            for jj=0, npoint-1 do begin
                index = where(tzs[jj,*] ge slope_brightness and mlon_bins gt slope_separator, count)
                if count eq 0 then continue
                tys[jj] = mlon_bins[index[count-1]]
            endfor
        endif else begin
            ; Westward.
            tlabel = 'west'
            for jj=0, npoint-1 do begin
                index = where(tzs[jj,*] ge slope_brightness and mlon_bins lt slope_separator, count)
                if count eq 0 then continue
                tys[jj] = mlon_bins[index[0]]
            endfor
        endelse

        index = where(finite(tys))
        txs = txs[index]
        tys = tys[index]
        plots, txs, tys, psym=1, symsize=0.1, color=sgcolor('white')
        fit_coef = linfit(txs, tys)
        linfit_yrange = yrange
        linfit_xrange = (linfit_yrange-fit_coef[0])/fit_coef[1]
        linfit_xrange = time_range_slope-(time_range_slope mod 60)+[-0.5,0.4]*60
        linfit_yrange = fit_coef[1]*linfit_xrange+fit_coef[0]
        oplot, linfit_xrange, linfit_yrange, linestyle=1

        tx = tpos[2]-xchsz*15
        ty = (tlabel eq 'east')? tpos[1]+ychsz*1: tpos[3]-ychsz*2
        xyouts, tx,ty, /normal, alignment=1, 'V!D'+tlabel+' !N: '+sgnum2str(abs(fit_coef[1]-earth_rotation), nsgn=2)+' deg/sec', charsize=label_size
    endfor



    ; Color bar.
    tpos = colorbar_pos
    sgcolorbar, colors=colors, ct=aurora_color_table, position=tpos, $
        /horizontal, ztitle=ztitle, zticks=zticks, charsize=label_size, zrange=zrange

    ; Add labels.
    fig_labels = ['a','b']+'.'
    for ii=0, nvar-1 do begin
        tpos = poss[*,ii]
        xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*1, /normal, 'Ewogram in '+ewo_infos[ii].label+' deg MLat', charsize=label_size
    endfor

    if keyword_set(test) then stop
    sgclose
end
