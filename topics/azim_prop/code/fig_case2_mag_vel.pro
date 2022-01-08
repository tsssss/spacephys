
;---Get basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    einfo = get_var_data(info_var)
    time_range = einfo.time_range_plot1

    test = 0

;---Open plot.
    file = sparentdir(srootdir())+'/plot/fig_case2_mag_vel.pdf'
    if keyword_set(test) then file = 0
    xsize = 3
    ysize = 4
    sgopen, file, xsize=xsize, ysize=ysize, magnify=2

;---General settings.
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    pos_keo_mag = [0.2,0.45,0.95,0.97]   ; for mag data.
    pos_keo_vel = [0.2,0.15,0.95,0.44]   ; for ewogram.
    ;pos_mag_loc = [0.2,0.10,0.95,0.30]   ; for mag locations.

    ticklen = -0.01
    label_size = 0.8
    arrow_size = (size(file,/type) eq 7)? 160:8
    str_deg = '!9'+string(176b)+'!X'

    shift_val = -400         ; nT
    mag_yrange = [shift_val*9,500]    ; y-range for mag data.

    first_color = 200       ; start color for mag data.
    last_color = 50         ; last color for mag data.
    data_rate_asi = 3d      ; sec.
    data_rate_mag = 0.5d    ; sec.
    mag_var = 'thg_dbh'     ; name for the mag data.

;---Ewo/keogram settings.
    asi_top_color = 254     ; top color for scaling ewo/keogram.
    asi_color_table = 60    ; color table for ewo/keogram.

    ; Settings for ewogram.
    ewo_color_table = 65        ; color table for mag data.
    ewo_mlat_range = [67,69]        ; for calc ewogram.
    mag_ewo_range = [67,69]         ; for choosing mag sites.
    ewo_mlon_range = [-5,-125]     ; for plotting ewogram.
    ewo_count_range = [10,150]      ; for scaling ewogram.

    ; Settings for keogram.
    keo_color_table = 62
    keo_mlat_range = [60.,75]       ; for plotting keogram.
    keo_mlon_range = [-95,-80]      ; for calc keogram.
    mag_keo_range = [-95,-80]       ; for choosing mag sites.
    keo_count_range = [10,200]


;---Collect info for all available sites.
    get_data, mag_var, times, mag_data, sites
    sites = strlowcase(sites)
    nsite = n_elements(sites)

    tinfo = {name:'', time:0d, mlon:0., mlat:0., base_val:0.}
    site_infos = replicate(tinfo, nsite)
    tinfos = themis_read_mag_metadata(sites=sites)
    for ii=0, nsite-1 do begin
        site_infos[ii].name = tinfos[ii].id
        site_infos[ii].mlon = tinfos[ii].mlon
        site_infos[ii].mlat = tinfos[ii].mlat
    endfor


    test_var = 'ewo_test'
    foreach test_site, sites, ii do begin
        index = where(site_infos.name eq test_site)
        f0 = mag_data[*,index[0]]
        store_data, test_var, times, f0, limits={labels:test_site}
        f1 = deriv(f0)/data_rate_mag
        sdespike, times, f1
        store_data, test_var+'df', times, f1, limits={constant:0}


        dr0 = 0.5
        width = 30  ; sec.
        dn0 = width/dr0
        uts = make_bins(time_range, width)
        nut = n_elements(uts)
        data = interpol(f0,times, uts)
        slopes = (data[1:nut-1]-data[0:nut-2])/width
        uts = (uts[1:nut-1]+uts[0:nut-2])*0.5
        store_data, test_var+'slope', uts, slopes

        ; use data before 10:05 UT to derive the base slope
        base_time = time_double('2014-08-28/10:06')
        index = where(uts le base_time)
        base_slope = mean(abs(slopes[index]))

        ; find the time of negative slope that is large enough.
        foreach neg_slope, -base_slope*[6] do begin
            index = where(slopes le neg_slope and uts ge base_time, count)
            if count ne 0 then begin
                base_time2 = uts[index[0]]
                site_infos[ii].time = base_time2
                break
            endif
        endforeach

        ; find the time of value change that is large enough.
        base_value = sinterpol(f0, times, base_time2)
        neg_value = base_value+(min(f0[where(times gt base_time2)])-base_value)*0.05
        index = where(f0 le neg_value and times ge base_time2, count)
        if count ne 0 then site_infos[ii].time = times[index[0]]

;        tplot, test_var+'*'
;;        ;site_infos[ii].time = time_double(site_str_infos[1,ii])
;        timebar, site_infos[ii].time
;        stop
    endforeach
    for ii=0, nsite-1 do print, sites[ii], '    ', time_string(site_infos[ii].time)





;---Plot keogram and the associated mag data.
    ; Filter out the wanted mag sites.
    index = lazy_where(site_infos.mlon, 'in', mag_keo_range, count=nkeo_site)
    keo_infos = site_infos[index]
    keo_mag_data = mag_data[*,index]

    keo_str_infos = [$
        'inuv', '2014-08-28/10:17:20', $
        'kako', '2014-08-28/10:16:40', $
        'fykn', '2014-08-28/10:11:05', $
        'eagl', '2014-08-28/10:09:40', $
        'cigo', '2014-08-28/10:08:50', $
        'pokr', '']
        ;'cmo' , '']
        ;'whit', '2014-08-28/10:08:45']
        ;'trap', '2014-08-28/10:18:00']
    nkeo_str_site = n_elements(keo_str_infos)/2
    keo_str_infos = reform(keo_str_infos,2,nkeo_str_site)

    flags = bytarr(nkeo_site)
    for ii=0, nkeo_site-1 do begin
        index = where(keo_str_infos[0,*] eq keo_infos[ii].name, count)
        if count eq 0 then flags[ii] = 1
        ;keo_infos[ii].time = time_double(keo_str_infos[1,index[0]])
    endfor
    index = where(flags ne 1, nkeo_site)
    keo_infos = keo_infos[index]
    keo_mag_data = keo_mag_data[*,index]

    ; Sort by mlat.
    index = reverse(sort(keo_infos.mlat))
    keo_infos = keo_infos[index]
    keo_mag_data = keo_mag_data[*,index]

    for ii=1, nkeo_site-1 do begin
        base_val = (min(keo_mag_data[*,ii-1])-max(keo_mag_data[*,ii]))*0.9
        base_val <= shift_val*0.5   ; prevent from to small.
        base_val = shift_val*total(minmax(keo_mag_data[*,ii-1])*[-1,1])/550+keo_infos[ii-1].base_val
        keo_infos[ii].base_val = base_val
    endfor


    ;---Plot the data.
        xrange = time_range
        xtickv = smkarthm(xrange[0],xrange[1],1200,'dx')
        xticks = n_elements(xtickv)-1
        xminor = 5
        xtickn = time_string(xtickv,tformat='hh:mm')
        xtitle = 'Time (UT) on '+time_string(xrange[0],tformat='YYYY-MM-DD')
        yrange = mag_yrange
        tpos = pos_keo_mag
        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xticklen=ticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickformat='(A1)', $
            ystyle=1, yrange=yrange, yticklen=ticklen, ytickformat='(A1)', $
            /noerase, /nodata
        xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*1, /normal, 'a. B!Deast!N'
        xyouts, tpos[2]-xchsz*0.5, tpos[3]-ychsz*1, /normal, alignment=1, 'Ground magnetometers!Cin '+$
            '['+sgnum2str(mag_keo_range[0])+str_deg+','+sgnum2str(mag_keo_range[1])+str_deg+'] MLon', charsize=label_size

        ; Get the colors for each site.
        colors = smkarthm(first_color, last_color, nkeo_site, 'n')
        for ii=0, nkeo_site-1 do colors[ii] = sgcolor(colors[ii], ct=keo_color_table)

        for ii=0, nkeo_site-1 do begin
            ;base_value = shift_val*ii
            base_value = keo_infos[ii].base_val
            site = keo_infos[ii].name
            color = colors[ii]

            data = mag_data[*,where(sites eq site)]+base_value
            oplot, times, data, color=color
            ; Position of the tart point of the data.
            tmp = convert_coord(xrange[0],base_value, /data, /to_normal)
            tx1 = tmp[0]
            ty1 = tmp[1]
            ; Position of the site name.
            tx2 = tx1-xchsz*1
            ty2 = ty1-ychsz*0.3
            xyouts, tx2,ty2, /normal, strupcase(site), color=color, alignment=1, charsize=label_size

            ; Add mark for the start time of FAC.
            tmp = convert_coord(keo_infos[ii].time,interpol(data,times,keo_infos[ii].time), /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            ;plots, tx+[0,0], ty+[-1,1]*ychsz*0.5, /normal;, color=color
            arrow, tx, ty-ychsz*0.5, tx, ty-ychsz*0.1, /normal, /solid, hsize=arrow_size
        endfor

        ; Add legend.
        tmp = convert_coord(xrange[1],yrange[1], /data, /to_normal)
        ty1 = tmp[1]
        tmp = convert_coord(xrange[1],yrange[1]-shift_val, /data, /to_normal)
        ty2 = tmp[1]
        dy = abs(ty1-ty2)   ; the height of the legend.
        ; Position of the legend.
        tx = tpos[0]+(tpos[2]-tpos[0])*0.5-xchsz*2
        ty = tpos[1]+(tpos[3]-tpos[1])*0.5+ychsz*1
        tys = ty+[-1,1]*dy*0.5
        plots, tx+[0,0], tys, /normal
        foreach ty, tys do plots, tx+[-1,1]*xchsz*0.5, ty, /normal
        xyouts, tx+xchsz*1, mean(tys)-ychsz*0.3*label_size, /normal, string(abs(shift_val),format='(I0)')+' nT', charsize=label_size

    ;---Second panel.
        tpos = pos_keo_vel
        yrange = keo_mlat_range
        yticks = 3
        yminor = 5
        ytitle = 'MLat (deg)'

        xrange = [time_range[0],time_double('2014-08-28/10:46')]
        xtickv = smkarthm(xrange[0],xrange[1],1200,'dx')
        xticks = n_elements(xtickv)-1
        xminor = 5
        xtickn = time_string(xtickv,tformat='hh:mm')
        xtitle = 'Time (UT) on '+time_string(xrange[0],tformat='YYYY-MM-DD')
        tx = (xrange[1]-xrange[0])/(time_range[1]-time_range[0])*(tpos[2]-tpos[0])+tpos[0]
        colorbar_pos = tpos & colorbar_pos[0] = tx+xchsz*0.5 & colorbar_pos[2] = colorbar_pos[0]+xchsz*0.8
        tpos[2] = tx
        
        plot, xrange, yrange, position=tpos, $
            xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, /noerase, /nodata

        ; Keogram.
        keo_var = 'thg_mlonimg_keo'
        themis_read_mlonimg_keo, 'thg_mlonimg', to=keo_var, mlon_range=keo_mlon_range
        get_data, keo_var, uts, keo, mlat_bins
        timg = bytscl(keo, min=keo_count_range[0], max=keo_count_range[1], top=asi_top_color)
        ttpos = [(convert_coord(min(uts),min(mlat_bins), /data, /to_normal))[0:1], $
            (convert_coord(max(uts),max(mlat_bins), /data, /to_normal))[0:1]]
        sgtv, timg, position=ttpos, ct=asi_color_table, /resize
        linfit_xrange = minmax(uts)


        ; Magnetometers.
        psym = 6
        symsize = 0.5
        linestyle = 1
        txs = keo_infos.time
        tys = keo_infos.mlat
        for ii=0, nkeo_site-1 do plots, txs[ii], tys[ii], psym=psym, symsize=symsize, color=colors[ii]

        ; Fit upward motion.
        up_mlat = 65d
        index = where(keo_infos.mlat gt up_mlat)
        txs = (keo_infos.time)[index]
        tys = (keo_infos.mlat)[index]
        res = linfit(txs,tys)
        linfit_yrange = [60,70]+2.5
        linfit_xrange = (linfit_yrange-res[0])/res[1]
        oplot, linfit_xrange, linfit_yrange, linestyle=linestyle
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx, ty, /normal, alignment=1, charsize=label_size, 'V!Dnorth!N='+sgnum2str(abs(res[1]),nsgn=2)+' deg/sec'

        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xticklen=ticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtitle=xtitle, $
            ystyle=1, yrange=yrange, yticklen=ticklen, yticks=yticks, yminor=yminor, ytitle=ytitle, $
            /noerase, /nodata
        xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*1, /normal, 'b. Keogram'
        xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*2, /normal, $
            '['+sgnum2str(keo_mlon_range[0])+str_deg+','+sgnum2str(keo_mlon_range[1])+str_deg+'] MLon', charsize=label_size

    ;---Colorbar.
        sgcolorbar, findgen(asi_top_color), zrange=keo_count_range, ct=asi_color_table, ztitle='Photon count', position=colorbar_pos

;    ;---Magnetometer locations.
;        tpos = pos_mag_loc
;        yrange = [63,73]
;        xrange = total(yrange*[-1,1])*2*[-1,1]+mean(keo_mlon_range)
;        xtickv = [-100,-85,-70]
;        xticks = n_elements(xtickv)-1
;        ticklen = -0.01
;        xminor = 3
;        time_start = time_double('2014-08-28/10:36')
;        tmp = convert_coord(time_start,yrange[0], /data, /to_normal)
;        ty1 = tmp[1]
;        tmp = convert_coord(time_start,yrange[1], /data, /to_normal)
;        ty2 = tmp[1]
;        tx1 = tmp[0]
;        tx2 = tx1+(ty2-ty1)*ysize/xsize
;        tpos = [tx1,ty1,tx2,ty2]
;        plot, xrange, yrange, /nodata, position=tpos, $
;            xstyle=1, xtitle='MLon (deg)', xtickv=xtickv, xticks=xticks, xticklen=ticklen, xcharsize=label_size, $
;            ystyle=1, ytitle='MLat (deg)', $
;            /noerase
;        txs = keo_infos.mlon
;        tys = keo_infos.mlat
;        for ii=0, nkeo_site-1 do plots, txs[ii], tys[ii], psym=psym, symsize=symsize, color=colors[ii]


;---Clearn up
    if keyword_set(test) then stop
    sgclose

end
