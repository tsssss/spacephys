;+
; Plot the aurora and draw spacecraft and magnetometers on it.
;-


;---Load basic data.
    _2014_0828_10_load_data, event_info=event_info

    test = 1


;---Settings.
    aurora = event_info['aurora']
    gmag_info = event_info['gmag_info']
    time_range = aurora['time_range']

    ; The magnetometers to be considered.
    mlonimg_var = aurora['var']
    mag_sites = ['whit','eagl','fsim']
    nmag_site = n_elements(mag_sites)
    mag_colors = sgcolor(['red','orange','gold'])
    time_colors = sgcolor(['brown','peru','dark_khaki'])
    site_infos = themis_read_mag_metadata(sites=mag_sites)

    ; For calculating the photon count.
    dmlon = 2d
    dmlat = 1d
    method = 'max'
    count_range = [0,300]

    ; For drawing snapshots of aurora.
    mlat_range = [61d,71]
    mlon_range = [-100d,-55]
    mlon_bins = get_setting(mlonimg_var, 'mlon_bins')
    if mlon_range[0] lt min(mlon_bins) then mlon_range[0] = min(mlon_bins)
    if mlon_range[1] gt max(mlon_bins) then mlon_range[1] = max(mlon_bins)
    mlat_bins = get_setting(mlonimg_var, 'mlat_bins')
    if mlat_range[0] lt min(mlat_bins) then mlat_range[0] = min(mlat_bins)
    if mlat_range[1] gt max(mlat_bins) then mlat_range[1] = max(mlat_bins)
    aspect_ratio = abs((mlat_range[1]-mlat_range[0])/(mlon_range[1]-mlon_range[0])*2.2)

    ; For drawing the magnetometer data.
    mag_xrange = time_range
    mag_xticks = 3
    mag_xtickv = smkarthm(mag_xrange[0], mag_xrange[1], mag_xticks+1, 'n')
    mag_xtickn = time_string(mag_xtickv, tformat='hh:mm') & mag_xtickn[0] = 'UT'
    mag_xminor = 5
    str_deg = '!9'+string(176b)+'!X'



;---Calculate the photon count near the magnetometers.
    mag_var = 'thg_dbh'
    foreach site, mag_sites, ii do begin
        var1 = 'thg_'+site+'_photon_count'
        mlon = (gmag_info[site])['mlon']
        mlat = (gmag_info[site])['mlat']
        themis_read_mlonimg_photon_count, mlonimg_var, to=var1, $
            mlon_center=mlon, mlat_center=mlat, dmlon=dmlon, dmlat=dmlat, method=method

        var2 = 'thg_'+site+'_mag_bh'
        get_data, mag_var, times, data, sites
        tdat = data[*,where(strlowcase(sites) eq site)]
        store_data, var2, times, tdat, limits={ytitle:'', labels:'', $
            ytickformat:'(A1)', colors:mag_colors[ii]}
        print, minmax(tdat)
        print, sg_autolim(tdat)
;        tplot, [var2,var1]
;        stop
    endforeach


;---Configure the photon count data.
    gamma_yrange = count_range*2
    gamma_yticks = 2
    gamma_yminor = 5
    gamma_ytickv = smkarthm(gamma_yrange[0], gamma_yrange[1], gamma_yticks+1, 'n')
    gamma_ytickn = string(gamma_ytickv, format='(I0)')

    vars = 'thg_'+mag_sites+'_photon_count'
    options, vars, 'labels', ''
    options, vars, 'ytitle', ''
    options, vars, 'ytickformat', '(A1)'
    options, vars, 'yrange', gamma_yrange
    options, vars, 'ytickv', gamma_ytickv
    options, vars, 'yticks', gamma_yticks
    options, vars, 'ytickn', gamma_ytickn
    options, vars, 'yminor', gamma_yminor



    mag_yticks = 2
    mag_yminor = 5
    foreach site, mag_sites do begin
        var = 'thg_'+site+'_mag_bh'
        get_data, var, times, data

        yrange = sg_autolim(data)
        ytickv = smkarthm(yrange[0], yrange[1], mag_yticks+1, 'n')
        yticks = n_elements(ytickv)-1
        ytickn = string(ytickv, format='(I0)')
        yminor = mag_yminor
        add_setting, var, {$
            yrange: yrange, $
            ytickv: ytickv, $
            yticks: yticks, $
            yminor: yminor}
    endforeach


;---Prepare for plotting.
    figure_size = [4.5,4.5]

    magnify = 2
    ticklen = -0.02
    aurora_color_table = 49
    top_color = 254
    colors = findgen(top_color+1)

    ztitle = 'Photon Count'
    label_size = 0.8

    fig_labels = ['a','b','c']


    wid = 0
    sgopen, wid, xsize=figure_size[0], ysize=figure_size[1], /inch, magnify=magnify
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    poss = sgcalcpos(nmag_site,2, lmargin=8, rmargin=6, tmargin=4, bmargin=4, xpad=1, ypad=1.5)
    xoffset = xchsz*5
    poss[2,0,*] += xoffset
    poss[0,1,*] += xoffset
    colorbar_pos = poss[[0,3,2,3],0,0]+[0,0.5,0,1]*ychsz
    sgclose
    wdelete, wid

;---Plot.
    get_data, mlonimg_var, mlonimg_times, mlonimgs
    figdir = sparentdir(srootdir())
    file = figdir+'/plot/fig_case2_mag_vs_aurora.pdf'
    if keyword_set(test) then file = 0

    sgopen, file, xsize=figure_size[0], ysize=figure_size[1], /inch, magnify=magnify


    foreach site, mag_sites, jj do begin
        fac_time = (gmag_info[site])['fac_time']

    ;---Draw the auroral snapshot.
        timg = get_var_data(mlonimg_var, at=fac_time)
        timg = bytscl(timg, min=count_range[0], max=count_range[1], top=top_color)
        index = where_pro(mlon_bins, mlon_range)
        timg = timg[index,*]
        index = where_pro(mlat_bins, mlat_range)
        timg = timg[*,index]

        figure_pos = poss[*,0,jj]
        tpos = figure_pos
        sgtv, timg, position=tpos, /resize, ct=aurora_color_table

        xtitle = 'MLon (deg)'
        xrange = mlon_range
        xticklen = ticklen
        xminor = 5
        dxmajor = 10
        xtickv = smkarthm(xrange[0],xrange[1],dxmajor,'dx')
        xticks = n_elements(xtickv)-1
        xtickn = string(xtickv,format='(I0)') & xtickn[0] = xtitle+'          '
        if jj ne nmag_site-1 then xtickn[*] = ' '   ; only draw x-axis for the last row.

        ytitle = 'MLat (deg)'
        yrange = mlat_range
        yticklen = ticklen*aspect_ratio
        yminor = 5
        dymajor = 5
        ytickv = make_bins(yrange, dymajor, /inner)
        yticks = n_elements(ytickv)-1
        ytickn = string(ytickv,format='(I0)')

        plot, xrange, yrange, /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
        if keyword_set(add_grid) then begin
            plot, xrange, yrange, /nodata, position=tpos, $
                /noerase, $
                xstyle=1, xrange=xrange, xticks=xticks, xminor=0, xtickv=xtickv, xtitle='', xgridstyle=1, xtickformat='(A1)', xticklen=1, $
                ystyle=1, yrange=yrange, yticks=yticks, yminor=0, ytickv=ytickv, ytitle='', ygridstyle=1, ytickformat='(A1)', yticklen=1, $
                color=sgcolor('white')
        endif

        ; Add axes.
        axis, yaxis=0, ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, ytickv=ytickv, yticklen=yticklen, ytitle='', ytickname=ytickn
        axis, yaxis=1, ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, ytickv=ytickv, yticklen=yticklen, ytitle='', ytickformat='(A1)'
        axis, xaxis=0, xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xtickv=xtickv, xticklen=xticklen, xtitle='', xtickname=xtickn
        axis, xaxis=1, xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xtickv=xtickv, xticklen=xticklen, xtitle='', xtickformat='(A1)'
        ; Add title.
        ;xyouts, tpos[0]-xchsz*6, tpos[1]-ychsz*1.8, /normal, xtitle, alignment=0
        xyouts, tpos[0]-xchsz*5, (tpos[1]+tpos[3])*0.5, /normal, ytitle, alignment=0.5, orientation=90
        ; Add label.
        xyouts, tpos[2]-xchsz*0.5*label_size, tpos[1]+ychsz*0.5*label_size, /normal, alignment=1, $
            time_string(fac_time,tformat='hh:mm:ss UT'), charsize=label_size, color=time_colors[jj]
        xyouts, tpos[2]-xchsz*0.5, alignment=1, tpos[3]-ychsz*1, /normal, fig_labels[jj]+'-1'

        ; Draw color bar.
        tpos = colorbar_pos
        zrange = count_range
        zticks = 3
        sgcolorbar, colors, zrange=zrange, ztitle=ztitle, zticks=zticks, $
            position=tpos, ct=aurora_color_table, /horizontal, zcharsize=label_size


    ;---Add magnetometers.
        psym = 6
        symsize = 0.5
        color = sgcolor('orange')
        color2 = sgcolor('red')
        tpos = figure_pos
        plot, xrange, yrange, /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
        foreach site, mag_sites, ii do begin
            mlon = site_infos[ii].mlon
            mlat = site_infos[ii].mlat
            plots, mlon, mlat, psym=psym, color=mag_colors[ii], symsize=symsize
        endforeach

        nouttick = (jj eq nmag_site-1)? 0: 1
        vars = 'thg_'+mag_sites[jj]+'_'+['mag_bh','photon_count']
        nvar = n_elements(vars)
        sub_pos = sgcalcpos(nvar, position=poss[*,1,jj])

    ;---Plot magnetometer data.
        tpos = sub_pos[*,0]
        get_data, vars[0], xs, ys, limits=lim
        xrange = mag_xrange
        xticks = mag_xticks
        xtickv = mag_xtickv
        xtickn = mag_xtickn
        xminor = mag_xminor
        xticklen = ticklen*2
        xticks = mag_xticks
        yrange = lim.yrange
        yticks = lim.yticks
        ytickv = lim.ytickv
        ytickn = string(ytickv,format='(I0)')
        yminor = lim.yminor
        yticklen = ticklen
        ytitle = '(nT)'
        yscale = 400

        plot, xrange, yrange, /noerase, /nodata, position=tpos, xstyle=5, ystyle=5
        oplot, xs, ys, color=mag_colors[jj]
        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickformat='(A1)', $
            ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytickformat='(A1)', $
            /noerase, /nodata
        ;axis, yaxis=1, ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytickn=ytickn, ytitle=ytitle

        xyouts, tpos[2]+xchsz*0.5, tpos[3]-ychsz*0.8, /normal, alignment=0, 'B!Deast!N'
        tx = tpos[2]+xchsz*0.5
        tys = yrange[0]+[0,yscale]
        foreach ty, tys, kk do begin
            tmp = convert_coord([xrange[1],ty], /data, /to_normal)
            tys[kk] = tmp[1]
        endforeach
        ynudge = (mean(tpos[[1,3]])-max(tys))*0.5
        if ynudge gt 0 then tys += ynudge
        plots, tx+[1,1]*xchsz*0.25, tys, /normal
        foreach ty, tys do plots, tx+[0,0.5]*xchsz, ty+[0,0], /normal
        xyouts, tx+xchsz*1, mean(tys)-ychsz*0.3*label_size, /normal, charsize=label_size, string(yscale,format='(I0)')+' nT'

        xyouts, tpos[2]-xchsz*0.5*label_size, tpos[3]-ychsz*1.0*label_size, /normal, alignment=1, charsize=label_size, strupcase(mag_sites[jj]), color=mag_colors[jj]
        xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.5, /normal, fig_labels[jj]+'-2'

    ;---Add notation.
        if jj eq 0 then begin
            tx = mean(tpos[[0,2]])
            ty = tpos[3]+ychsz*2.2
            xyouts, tx, ty, /normal, alignment=0.5, charsize=label_size, $
                'Magnetometer data vs!Cmax auroral brightness!C within 100 km x 100 km'
                ;'Magnetometer data vs!Cauroral brightness w/i '+strjoin(string([dmlon,dmlat],format='(I0)')+str_deg,'x')

        endif

    ;---Plot auroral brightness around the magnetometer.
        tpos = sub_pos[*,1]
        get_data, vars[1], xs, ys, limits=lim
        xtickformat = (jj eq nmag_site-1)? '': '(A1)'
        yrange = lim.yrange
        yticks = lim.yticks
        ytickv = lim.ytickv
        ytickn = lim.ytickn
        yminor = lim.yminor
        yticklen = ticklen
        ytitle = '(#)'

        plot, xrange, yrange, /noerase, /nodata, position=tpos, xstyle=5, ystyle=5
        oplot, xs, ys
        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickformat=xtickformat, xtickn=xtickn, $
            ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytickformat='(A1)', $
            /noerase, /nodata
        axis, yaxis=1, ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytickn=ytickn, ytitle=ytitle

        xyouts, tpos[2]-xchsz*0.5*label_size, tpos[3]-ychsz*1*label_size, /normal, 'Brightness', alignment=1, charsize=label_size
        xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*1, /normal, fig_labels[jj]+'-3'

    ;---Add the time bar.
        tmp = convert_coord([fac_time,0], /data, /to_normal)
        txs = tmp[0]+[0,0]
        tys = [sub_pos[1,1],sub_pos[3,0]]
        plots, txs,tys, /normal, color=time_colors[jj]
    endforeach

    sgclose

end
