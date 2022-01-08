;+
; Plot Mag H component for all loaded gmag, and plot the positions of them.
;-

pro fig_case2_mag_new, add_themis=add_themis

;---Load basic data.
    _2014_0828_10_load_data, event_info=event_info

    test = 0

;---Settings.
    time_range = time_double(['2014-08-28/09:50','2014-08-28/11:00'])


;---Load data and settings.
    mag_var = 'thg_dbh'
    get_data, mag_var, times, gmag_data
    ntime = n_elements(times)
    gmag_sites = event_info.gmag.sites      ; all gmag sites that are loaded.
    gmag_info = event_info.gmag_info        ; info for all gmag sites.
    gmag = event_info.gmag


;---Filter.
    mlon_range = gmag.mlon_range
    mlat_range = gmag.mlat_range
    mlons = []
    mlats = []
    foreach site, gmag_sites do begin
        mlons = [mlons, gmag_info[site].mlon]
        mlats = [mlats, gmag_info[site].mlat]
    endforeach

    used_sites = gmag_sites
    index = lazy_where(mlons, mlon_range)
    used_sites = used_sites[index]
    index = lazy_where(mlats, mlat_range)
    used_sites = used_sites[index]

    index = []
    foreach site, used_sites do index = [index, where(gmag_sites eq site)]
    gmag_data = gmag_data[*,index]
    mlons = mlons[index]
    mlats = mlats[index]

;---Sort.
    index = sort(mlons)
    gmag_data = gmag_data[*,index]
    mlons = mlons[index]
    mlats = mlats[index]
    used_sites = used_sites[index]


;---Combine data.
    used_gmag = dictionary()
    used_gmag.site = used_sites
    used_gmag.color_table = 40

    base_value_shift = 200d    ; nT. The basic shift when combining mag data.
    color_start = 100
    color_end = 250

    ; Get the sites, sort by mlon.
    sites = used_sites
    nsite = n_elements(sites)

    ; Collect data and shift to avoid overlapping.
    foreach site, sites, ii do gmag_data[*,ii] -= gmag_data[0,ii]

    base_values = fltarr(nsite)
    base_values[0] = 0
    for ii=1, nsite-1 do begin
        base_values[ii] = base_values[ii-1]-max(gmag_data[*,ii]-gmag_data[*,ii-1])-base_value_shift
    endfor
    used_gmag.base_value = base_values

    for ii=0, nsite-1 do gmag_data[*,ii] += base_values[ii]
    yrange = minmax(gmag_data)+[-1.5,0.5]*base_value_shift

    ; Save the combined var.
    used_gmag.var = 'gmag_combo'
    colors = round(smkarthm(color_start,color_end, nsite, 'n'))
    foreach color, colors, ii do colors[ii] = sgcolor(colors[ii], ct=used_gmag.color_table)
    store_data, used_gmag.var, times, gmag_data
    add_setting, used_gmag.var, /smart, {$
        display_type:'stack', $
        labels:strupcase(sites), $
        colors:colors, $
        yrange:yrange}


;---Plot the data.
    if keyword_set(test) then begin
        file = 0
        magnify = 1.5
    endif else begin
        file = sparentdir(srootdir())+'/plot/fig_case2_mag.pdf'
        if keyword_set(add_themis) then file = sparentdir(srootdir())+'/plot/fig_case2_mag_with_themis.pdf'
        magnify = 1
    endelse

    sgopen, file, xsize=4, ysize=8, magnify=magnify



    poss = sgcalcpos(2, ypans=[5,1], xchsz=xchsz, ychsz=ychsz, tmargin=2, rmargin=4, lmargin=8, bmargin=5, ypad=5)
    pos1 = poss[*,0]
    pos2 = poss[*,1]

    tpos = pos1

    tvar = used_gmag.var
    options, tvar, 'ytickformat', '(A1)'
    options, tvar, 'ytitle', ''
    options, tvar, 'labels', ''
    options, tvar, 'xticklen', -0.02/2
    tplot, used_gmag.var, trange=time_range, /novtitle, position=tpos

    xrange = time_range
    label_size = 0.8

    ; Add gmag name on the left.
    plot, xrange, yrange, /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
    for ii=0, nsite-1 do begin
        tx = xrange[0]
        ty = mean(gmag_data[0:100,ii])
        tmp = convert_coord(tx,ty, /data, /to_normal)
        tx = tmp[0]-xchsz*5*label_size
        ty = tmp[1]-ychsz*0.4*label_size
        xyouts, tx,ty, /normal, alignment=0, strupcase(sites[ii]), color=colors[ii], charsize=label_size
    endfor

    ; Add scale.
    gmag_scale = base_value_shift*2
    tx = xrange[0]
    ty = yrange[0]+[0,gmag_scale]
    tmp = convert_coord(tx,ty[0], /data, /to_normal)
    ty1 = tmp[1]
    tmp = convert_coord(tx,ty[1], /data, /to_normal)
    ty2 = tmp[1]
    tx = tmp[0]+xchsz*2
    tys = [ty1,ty2]+ychsz*0.5
    plots, tx+[0,0], tys, /normal
    foreach ty, tys do plots, tx+[-1,1]*0.2*xchsz, ty+[0,0], /normal
    xyouts, tx+xchsz*1, mean(tys)-ychsz*0.3, /normal, sgnum2str(gmag_scale)+' nT'
    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, 'a. dB!Deast!N (nT), ground magnetometer'


    ; Add position.
    tpos = pos2
    sym = 1
    xrange = mlon_range
    xminor = 6
    xticks = (xrange[1]-xrange[0])/(xminor*5)
    yrange = mlat_range
    yminor = 5
    yticks = (yrange[1]-yrange[0])/yminor
    plot, xrange, yrange, /nodata, /noerase, $
        xstyle=1, xticks=xticks, xticklen=-0.03, xminor=xminor, xtitle='MLon (deg)', $
        ystyle=1, yticks=yticks, yticklen=-0.01, yminor=yminor, ytitle='MLat (deg)', $
        position=tpos
    foreach site, sites, ii do begin
        tx = gmag_info[site].mlon
        ty = gmag_info[site].mlat
        plots, tx,ty, /data, psym=sym, symsize=label_size, color=colors[ii]
        tmp = convert_coord(tx,ty, /data, /to_normal)
        tx = tmp[0]+xchsz*1*label_size
        ty = tmp[1]-ychsz*0.3*label_size
        xyouts, tx,ty,/normal, strupcase(site), color=colors[ii], charsize=label_size
    endforeach
    xyouts, tpos[0], tpos[3]+ychsz*0.5, /normal, 'b. Ground magnetometer locations'

    ; Add themis.
    if keyword_set(add_themis) then begin
        themis = dictionary()
        themis_probes = ['a','d','e']
        themis_colors = sgcolor(['black','black','black'])
        themis_times = time_double('2014-08-28/'+['10:20','10:10','10:20'])
        model_chosen = 't89'
        sym = 6
        foreach probe, themis_probes, ii do begin
            var = 'th'+probe+'_fmlon_'+model_chosen
            tx = get_var_data(var, at=themis_times[ii])
            var = 'th'+probe+'_fmlat_'+model_chosen
            ty = get_var_data(var, at=themis_times[ii])
            plots, tx,ty,/data, psym=sym, symsize=label_size, color=themis_colors[ii]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]+xchsz*1*label_size
            ty = tmp[1]-ychsz*0.3*label_size
            xyouts, tx,ty,/normal, strupcase('th'+probe), color=themis_colors[ii], charsize=label_size
        endforeach
    endif


    if keyword_set(test) then stop
    sgclose
end

foreach flag, [0,1] do fig_case2_mag, add_themis=flag
end
