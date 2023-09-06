;+
; For THD, calc pflux, and time lag it, and calculate the cross correlation to auroral brightness for a given MLat and MLon range.
;-

pro fig_case2_pflux_corr_map, models=models, mlon_range=mlon_range, mlat_range=mlat_range

;---Load data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data
    einfo = get_var_data(info_var)

    test = 1


;---Settings.
    time_range_plot = time_double(['2014-08-28/10:00','2014-08-28/10:20'])

    ; The probe.
    probe = 'd'
    prefix = 'th'+probe
    pre0 = prefix+'_'

    ; The MLon/MLat range.
    if n_elements(mlon_range) ne 2 then mlon_range = [-90,-70]
    ;if n_elements(mlon_range) ne 2 then mlon_range = [-80,-70]
    ;if n_elements(mlon_range) ne 2 then mlon_range = [-77,-72]
    if n_elements(mlat_range) ne 2 then mlat_range = [60,70]
    
    mlon_range = [-75,-70]
    mlat_range = [65,68]

    if keyword_set(test) then begin
        figfile = test
        magnify = 2
    endif else begin
        figfile = join_path([sparentdir(srootdir()),'plot','fig_case2_pflux_vs_aurora.pdf'])
        magnify = 1
    endelse


;---Get the time lag.
    trace_dir = -1
    trace_height = 100. ; km.
    re = 6378.  ; km.
    r0 = trace_height/re+1
    model = 't89'
    va0 = 1000. ; km/s.

    rvar = pre0+'r_gsm'
    time = mean(time_range_plot)
    rgsm = get_var_data(rvar, at=time)
    tilt = geopack_recalc(time)
    par = 2
    geopack_trace, rgsm[0],rgsm[1],rgsm[2], trace_dir, par, tilt=tilt, r0=r0, /igrf, /t89, $
        fx,fy,fz, fline=fline
    ; fline in [n,3], in Re, in GSM.
    nstep = n_elements(fline)/3
    bmag = fltarr(nstep)
    for jj=0, nstep-1 do begin
        ; Need to recalc as time elapse???
        geopack_igrf_gsm, fline[jj,0],fline[jj,1],fline[jj,2], bx,by,bz
        geopack_t89, par, fline[jj,0],fline[jj,1],fline[jj,2], dbx,dby,dbz
        bmag[jj] = snorm([bx,by,bz]+[dbx,dby,dbz])
    endfor

    alfven_travel_time = 0
    field_line_length = 0
    for jj=0, nstep-2 do begin
        tva = va0*(bmag[jj]+bmag[jj+1])*0.5/bmag[0]
        tlength = snorm(fline[jj,*]-fline[jj+1,*])
        tduration = tlength/tva*re
        alfven_travel_time += tduration
        field_line_length += tlength
    endfor


;---Calculate pflux.
    scale_info = {s0:10d, s1:1000, dj:1d/8, ns:0d}

    fac = ['||','east','north']
    rgb = sgcolor(['red','green','blue'])
    probes = ['a','d','e']
    pf_unit = 'mW/m!U2!N'

    vars = pre0+['b0_gsm','db_gsm','r_gsm']
    foreach var, vars do interp_time, var, to=pre0+'e_gsm'

    ; Convert E/B fields to FAC.
    define_fac, pre0+'b0_gsm', pre0+'r_gsm'
    vars = pre0+['db_gsm','e_gsm']
    foreach var, vars do to_fac, var

    vars = pre0+['db','e']+'_fac'
    short_names = ['dB','dE']
    foreach var, vars, ii do begin
        add_setting, var, /smart, {$
            display_type: 'vector', $
            short_name: short_names[ii], $
            coord: 'FAC', $
            coord_labels: fac, $
            colors: rgb}
    endforeach

    ; Calculate pflux.
    stplot_calc_pflux_mor, pre0+'e_fac', pre0+'db_fac', pre0+'pf_fac', scaleinfo=scale_info;, /power

    ; Get parallel pflux and map to 100 km.
    pf_var = pre0+'pf_para_map'
    get_data, pre0+'pf_fac', times, data
    cmap = get_var_data(pre0+'cmap', at=times)
    pf_para = data[*,0]*cmap
    store_data, var, times, pf_para
    add_setting, var, /smart, {$
        display_type:'scalar', $
        short_name: 'S!D'+fac[0]+'!N', $
        unit:pf_unit}


;---Get auroral data.
    aurora_data_rate = 3.
    mlonimg_var = 'thg_mlonimg'
    get_data, mlonimg_var, times, mlonimgs

    ; Apply the time lag to aurora data, then crop in time.
    index = where_pro(times, time_range_plot+alfven_travel_time[0])
    mlonimgs = mlonimgs[index,*,*]
    times = times[index]

    ; Get the pflux para.
    pf_para = get_var_data(pf_var, at=times)

    ; Apply MLon/MLat range.
    mlon_bins = get_setting(mlonimg_var, 'mlon_bins')
    index = where_pro(mlon_bins, mlon_range, count=nmlon_bin)
    if nmlon_bin eq 0 then stop
    mlonimgs = mlonimgs[*,index,*]
    mlon_bins = mlon_bins[index]

    mlat_bins = get_setting(mlonimg_var, 'mlat_bins')
    index = where_pro(mlat_bins, mlat_range, count=nmlat_bin)
    if nmlat_bin eq 0 then stop
    mlonimgs = mlonimgs[*,*,index]
    mlat_bins = mlat_bins[index]

    fig_size = [nmlon_bin,nmlat_bin]
    corr_map = fltarr(fig_size)
    mlon_err = 2d
    mlat_err = 1d
    dmlon = median(mlon_bins-shift(mlon_bins,1))
    dmlat = median(mlat_bins-shift(mlat_bins,1))
    dmlon_rec = round(mlon_err/dmlon*0.5)
    dmlat_rec = round(mlat_err/dmlat*0.5)
    ntime = n_elements(times)
    for ii=dmlon_rec, fig_size[0]-1-dmlon_rec do begin
        for jj=dmlat_rec, fig_size[1]-1-dmlat_rec do begin
            photon_count = fltarr(ntime)
            for kk=0, ntime-1 do photon_count[kk] = $
                max(mlonimgs[kk,ii-dmlon_rec:ii+dmlon_rec,jj-dmlat_rec:jj+dmlat_rec])
            ; Remove background.
            photon_count = photon_count-(photon_count[0]+(photon_count[-1]-photon_count[0])*findgen(ntime)/(ntime-1))
            corr_map[ii,jj] = c_correlate(pf_para, photon_count, 0)
        endfor
    endfor

    corr_map_var = 'corr_map'
    store_data, corr_map_var, 0, corr_map
    add_setting, corr_map_var, {$
        mlon_bins: mlon_bins, $
        mlat_bins: mlat_bins}

    ; Get the photon count for the maximum correlation.
    max_corr = max(corr_map)
    index = where(corr_map eq max_corr)
    index = array_indices(fig_size, index, /dimensions)
    max_corr_mlon = mlon_bins[index[0]]
    max_corr_mlat = mlat_bins[index[1]]
max_corr_mlon = -75.2
max_corr_mlat = 65.8
tmp = min(mlon_bins-max_corr_mlon, idx1, /absolute)
max_corr_mlon = mlon_bins[idx1]
tmp = min(mlon_bins-max_corr_mlon, idx2, /absolute)
max_corr_mlon = mlon_bins[idx2]
index = [idx1,idx2]
    max_corr_photon_count = fltarr(ntime)
    for kk=0, ntime-1 do max_corr_photon_count[kk] = $
        max(mlonimgs[kk,index[0]-dmlon_rec:index[0]+dmlon_rec,index[1]-dmlat_rec:index[1]+dmlat_rec])
    max_corr_photon_count = max_corr_photon_count-(max_corr_photon_count[0]+(max_corr_photon_count[-1]-max_corr_photon_count[0])*findgen(ntime)/(ntime-1))
    photon_count_var = pre0+'photon_count'
    store_data, photon_count_var, times, max_corr_photon_count
    add_setting, photon_count_var, {$
        ytitle: '(Photon Count)', $
        mlon: max_corr_mlon, $
        mlat: max_corr_mlat}





;---Plot. Should have photon_count_var, corr_map_var, pf_var
    if keyword_set(test) then begin
        file = test
        magnify = 2
    endif else begin
        file = join_path([sparentdir(srootdir()),'plot','fig_case2_pflux_vs_aurora.pdf'])
        magnify = 1
    endelse

    sgopen, file, xsize=4, ysize=5, magnify=magnify
    poss = sgcalcpos(2, margin=[8,5,8,5], xchsz=xchsz, ychsz=ychsz, ypans=[3,2], ypad=5)

;---Panel 1.
    tpos = poss[*,0]

    zs = get_var_data(corr_map_var)
    xs = get_setting(corr_map_var, 'mlon_bins')
    ys = get_setting(corr_map_var, 'mlat_bins')

    zrange = [-1,1]*0.8
    ztitle = 'Cross Correlation (#)'
    nlevel = 41
    levels = smkarthm(zrange[0], zrange[1], nlevel, 'n')
    color_table = 70
    color_start = 250
    color_end = 5
    colors = smkarthm(color_start,color_end, nlevel, 'n')

    xrange = mlon_range
    xminor = 2
    xtickv = make_bins(xrange, xminor)
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    for ii=0, xticks do xtickn[ii] = sgnum2str(xtickv[ii])
    xtickn[0] = 'MLon (deg)'
    for ii=0, strlen(xtickn[0])-1 do xtickn[0] += ' '
    xticklen = -0.02
    xtitle = ''

    yrange = mlat_range
    yminor = 2
    ytickv = make_bins(yrange, yminor)
    yticks = n_elements(ytickv)-1
    ytickn = strarr(yticks+1)
    yticklen = -0.02
    ytitle = 'MLat (deg)'

    ; Add color table.
    cbpos = tpos[[0,3,2,3]]+ychsz*[0,0.5,0,1]
    sgcolorbar, colors, zrange=zrange, ztitle=ztitle, position=cbpos, /horizontal, ct=color_table

    ; Add contour.
    foreach color, colors, ii do colors[ii] = sgcolor(color, ct=color_table)
    contour, zs, xs, ys, /noerase, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickn=xtickn, xticklen=xticklen, xtitle=xtitle, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, ytickn=ytickn, yticklen=yticklen, ytitle=ytitle, $
        position=tpos, c_colors=colors, levels=levels, /fill

    ; Add model footpoint.
    if n_elements(models) eq 0 then models = ['t89','t96','t01','t04s']
    nmodel = n_elements(models)
    model_symbols = [7,4,1,6]
    model_colors = replicate(sgcolor('magenta'),nmodel)
    symsize = 0.5
    time = mean(time_range_plot)
    foreach model, models, ii do begin
        tx = get_var_data(pre0+'fmlon_'+model, at=time)
        ty = get_var_data(pre0+'fmlat_'+model, at=time)
        plots, tx, ty, /data, psym=model_symbols[ii], color=model_colors[ii], symsize=symsize
        tmp = convert_coord([tx,ty], /data, /to_normal)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.2
        xyouts, tx, ty, /normal, strupcase(model), charsize=symsize, alignment=0
    endforeach

    tmp = findgen(21)/20*2*!dpi
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs,tys
    psym = 8
    color = sgcolor('magenta')
    label = 'Max.Corr'
    tx = get_setting(photon_count_var, 'mlon')
    ty = get_setting(photon_count_var, 'mlat')
    plots, tx, ty, /data, psym=psym, color=color, symsize=symsize
    tmp = convert_coord([tx,ty], /data, /to_normal)
    tx = tmp[0]+xchsz*0.5
    ty = tmp[1]+ychsz*0.2
    xyouts, tx, ty, /normal, label, charsize=symsize, alignment=0



;---Panel 2.
    tpos = poss[*,1]
    get_data, photon_count_var, xs, ys
    pf_para = get_var_data(pf_var, at=times)

    log = 0

    xrange = time_range_plot
    xlog = log
    xminor = 5
    xstep = xminor*60    ; 5 min.
    xtickv = smkarthm(xrange[0],xrange[1],xstep,'dx')
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    for ii=0, xticks do xtickn[ii] = time_string(xtickv[ii],tformat='hhmm')
    xtickn[0] = time_string(xtickv[0],tformat='YYYY-MM-DD')
    for ii=0, strlen(xtickn[0])-1 do xtickn[0] += ' '

    ; For photon count.
    yrange = sg_autolim(ys)
    ylog = log
    ytickv = yrange
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytickn = strarr(yticks+1)
    for ii=0, yticks do ytickn[ii] = sgnum2str(ytickv[ii])
    ytitle = 'Photon Count (#)'

    plot, xs, ys, $
        xstyle=1, xlog=xlog, xrange=xrange, xticks=xticks, xtickv=xtickv, xtickn=xtickn, xminor=xminor, xticklen=xticklen, $
        ystyle=9, ylog=ylog, yrange=yrange, yticks=yticks, ytickv=ytickv, ytickn=ytickn, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /noerase, position=tpos

    ; For Poynting flux.
    ys = get_var_data(pf_var, at=xs)
    color = sgcolor('red')
    yrange = [1e0,1e2]
    ylog = log
    ytickv = [1e0,1e1,1e2]
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytickn = strarr(yticks+1)
    ytitle = 'S!D||!N@100km (mW/m!U!N)'
    for ii=0, yticks do ytickn[ii] = sgnum2str(ytickv[ii])
    axis, yaxis=1, /save, ylog=ylog, yrange=yrange, yticks=yticks, ytickv=ytickv, ytickn=ytickn, yminor=yminor, yticklen=yticklen, ytitle=ytitle, color=color
    oplot, xs, ys, color=color




    if keyword_set(test) then stop
    sgclose

end
