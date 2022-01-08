;+
; Plot gmag data in the format of EWOgram.
;-

    ;_2014_0828_10_load_data

    mlt_range = [-3,9]
    time_range = time_double(['2014-08-28/09:00','2014-08-28/14:00'])
    sort_time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
    sorted_probes = ['thd','the','g13','tha','rbspb','g15']
    panel_x_ratio = 0.8 ; inch per hour.
    gmag_zlog = 0
    gmag_range = [0,400]

    ; Add lines.
    line_color = sgcolor('grey')
    line_xx = list()
    line_yy = list()
    line_xx.add, time_double('2014-08-28/'+['10:05','11:10'])
    line_yy.add, [0.0,total(line_xx[0]*[-1,1])/60*2.1/15]
;    line_xx.add, time_double('2014-08-28/'+['16:30','17:10'])
;    line_yy.add, [0.0,9.0]
;    line_xx.add, time_double('2014-08-28/'+['15:40','16:30'])
;    line_yy.add, [0.0,9.0]


;    mlt_range = [-6,6]
;    time_range = time_double(['2013-06-07/04:00','2013-06-07/07:00'])
;    sort_time_range = time_double(['2013-06-07/04:50','2013-06-07/05:20'])
;    sorted_probes = ['thd','the','g13','rbspa','rbspb','g15']
;    panel_x_ratio = 0.8 ; inch per hour.
;    gmag_zlog = 0
;    gmag_range = [0,400]


    gmag_var = 'thg_dbh'
    smooth_window = 60*120
    long_time_range = time_range+[-1,1]*smooth_window
    if check_if_update(gmag_var, long_time_range) then themis_read_mag, long_time_range
    project = azim_df_load_project()
    azim_df_load_basic_data, project=project
    mlat_range = [60,80]


    line_linestyle = 1
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 1
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1
    rxy_range = [4.,30]
    root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
    ofn = join_path([root_dir,'fig_ewogram_of_dp_and_gmag'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hh'),'_')+'.pdf'])
test = 1

;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ; Get the sizes of the left panels.
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1,1,1]
    nypanel = n_elements(ypans)
    ypads = [0.4,0.4,2]
    panel_xsize = abs(total(time_range*[-1,1]))/3600*panel_x_ratio
    margins = [10.,4,8,1]
    fig_xsize = panel_xsize+(total(margins[[0,2]]))*abs_xchsz
    fig_ysize = panel_ysize*total(ypans)+(total(margins[[1,3]])+total(ypads))*abs_ychsz

;---Load data.
    if check_if_update('ae',time_range) then omni_read_index, time_range


;---Common x-axis.
    xrange = time_range
    xticklen_chsz = -0.2
    yticklen_chsz = -0.4
    xminor = 6 ; hr.
    xstep = 3600.
    xtickv = make_bins(xrange, xstep, /inner) ; make time line up.
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    secofday = 86400d
    for ii=0, xticks do begin
        the_time = xtickv[ii]
        xtickn[ii] = time_string(the_time,tformat='hh:mm')
        date = time_string(the_time,tformat='YYYY-MM-DD')
        if ii eq 0 then begin
            xtickn[ii] += '!C'+date
            continue
        endif
        if the_time mod secofday ne 0 then continue
        xtickn[ii] += '!C'+date
    endfor


    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz
    poss = sgcalcpos(nypanel, ypans=ypans, ypad=ypads, margins=margins)


    ; Label.
    fig_labels = letters(nypanel)+'. '+['AE','Dipolari-!C     zation','Gmag']
    for ii=0,nypanel-1 do begin
        tpos = poss[*,ii]
        tx = label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, fig_labels[ii]
    endfor



;---AE.
    the_var = 'ae'
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xtickformat = '(A1)'
    ystep = 1000.
    constants = [500]
    yys = get_var_data(the_var, in=time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    yminor = 5
    ytitle = '(nT)'

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Add data.
    oplot, xxs, yys

    ; Add grid.
    grid_color = sgcolor('gray')
    grid_linestyle = 1
    y_constants = [1000]
    x_constants = make_bins(time_range, 3600)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


;---Panel c. EWOgram for gmag.
    tpos = poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    ytitle = 'MLT (hr)'
    yrange = mlt_range
    yminor = 3
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0

    data_var_suffix = 'scaled_theta'
    ztitle = 'dB south'+' (nT)'
    if keyword_set(gmag_zlog) then begin
        zrange = alog10(gmag_range)
        log_ztickv = make_bins(zrange,1,/inner)
        ztickv = log_ztickv
        ztickn = '10!U'+string(log_ztickv,format='(I0)')
        index = where(log_ztickv eq 0, count)
        if count ne 0 then ztickn[index] = '1'
        index = where(log_ztickv eq 1, count)
        if count ne 0 then ztickn[index] = '10'
        index = where(log_ztickv eq -1, count)
        if count ne 0 then ztickn[index] = '0.1'
    endif else begin
        zrange = gmag_range
        zstep = abs(total(gmag_range*[-1,1])/2)
        ztickv = make_bins(zrange,zstep,/inner)
        ztickn = string(ztickv,format='(I0)')
    endelse
    zticks = n_elements(ztickv)-1

    ; Symbol.
    psym = 8
    symsize = 0.5
    usersym, [1,1,-1,-1,1]*0.1, [-1,1,1,-1,-1]*0.6
    ct = 49

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase

    pos_var_suffix = 'pseudo_mlt'
    flags = list()
    gmags = get_var_data(gmag_var, times=times, limits=lim)
    sites = lim.labels
    nsite = n_elements(sites)
    foreach site, sites, site_id do begin
        width = smooth_window/sdatarate(times)
        gmags[*,site_id] -= smooth(gmags[*,site_id], width, /edge_truncate, /nan)
    endforeach
    gmag_var2 = gmag_var+'_smooth'
    if keyword_set(gmag_zlog) then begin
        store_data, gmag_var2, times, alog10(-gmags), limits=lim
    endif else begin
        store_data, gmag_var2, times, -gmags, limits=lim
    endelse

    ; Site info.
    site_infos = themis_read_mag_metadata(sites=sites)
    mlons = fltarr(nsite)
    mlats = fltarr(nsite)
    foreach site_info, site_infos, site_id do begin
        mlons[site_id] = site_info.mlon
        mlats[site_id] = site_info.mlat
    endforeach

    gmags = get_var_data(gmag_var2, times=times, limits=lim, in=time_range)
    gmags = bytscl(gmags, max=zrange[1], min=zrange[0])
    ntime = n_elements(times)
    mlts = fltarr(ntime,nsite)+yrange[1]+10
    foreach site, sites, site_id do begin
        if mlats[site_id] ge mlat_range[1] then continue
        if mlats[site_id] le mlat_range[0] then continue
        mlts[*,site_id] = mlon2mlt(mlons[site_id], times)
    endforeach

    device, decomposed=0
    loadct, ct
    foreach time, times, time_id do begin
        yys = reform(mlts[time_id,*])
        xxs = dblarr(nsite)+time
        zzs = reform(gmags[time_id,*])
        ; Remove data outside yrange.
        index = lazy_where(yys,'[]', yrange, count=count)
        if count eq 0 then continue
        xxs = xxs[index]
        yys = yys[index]
        zzs = zzs[index]
        ; Sort according to color.
        index = sort(zzs)
        xxs = xxs[index]
        yys = yys[index]
        zzs = zzs[index]

        for ii=0, count-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach
    device, decomposed=1

    ; Add grid.
    y_constants = make_bins(mlt_range, 3, /inner)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle

    ; Add lines.
    for ii=0, line_xx.length-1 do begin
        oplot, line_xx[ii], line_yy[ii], color=line_color, linestyle=line_linestyle
    endfor

    ; Add box.
    xtickformat = ''
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Draw color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;---Panel b. UT-MLT long.
    tpos = poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    ytitle = 'MLT (hr)'
    yrange = mlt_range
    yminor = 3
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0

    data_var_suffix = 'scaled_theta'
    ztitle = 'Scaled '+tex2str('theta')+' (deg)'
    tilt_range = [-1,1]*64
    zrange = tilt_range
    ztickv = [-64,-16,-4,0,4,16,64]
    ztickn = string(ztickv,format='(I0)')
    ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
    zrange = alog10(abs(zrange))*[-1,1]
    zticks = n_elements(ztickv)-1

    ; Symbol.
    psym = 8
    symsize = 0.5
    usersym, [1,1,-1,-1,1]*0.1, [-1,1,1,-1,-1]*1
    ct = 70
    reverse_ct = 1

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase

    pos_var_suffix = 'pseudo_mlt'
    flags = list()
    foreach probe, reverse(sorted_probes) do begin
        prefix = probe+'_'
        pos_var = prefix+pos_var_suffix
        data_var = prefix+data_var_suffix
        zzs = get_var_data(data_var, in=time_range, times=xxs)
        yys = get_var_data(pos_var, at=xxs)
        mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
        zzs = azim_df_normalize_theta(zzs, zrange=tilt_range, ct=ct, /reverse_ct)

        ; Remove data out of range.
        r_sm = get_var_data(prefix+'r_sm', at=xxs)
        rxy = snorm(r_sm[*,0:1])
        ; Exclude data outside the MLT range.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = lazy_where(rxy, '][', rxy_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(r_sm)
        index = where(magn_flags eq 0, count)
        if count ne 0 then xxs[index] = !values.f_nan

        index = where(finite(xxs), count)
        the_flag = (count gt 0)
        flags.add, the_flag
        if count eq 0 then continue

        nxx = count
        xxs = xxs[index]
        yys = yys[index]
        zzs = zzs[index]
        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach

    ; Add grid.
    y_constants = make_bins(mlt_range, 3, /inner)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle

    ; Add lines.
    for ii=0, line_xx.length-1 do begin
        oplot, line_xx[ii], line_yy[ii], color=line_color, linestyle=line_linestyle
    endfor

    ; Add box.
    xtickformat = '(A1)'
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


    ; Draw color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks



    if keyword_set(test) then stop
    sgclose


end
