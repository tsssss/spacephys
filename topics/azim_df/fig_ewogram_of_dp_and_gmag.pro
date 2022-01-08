;+
; Plot gmag data in the format of EWOgram.
;-


test = 0

;    mlt_range = [-3,9]
;    time_range = time_double(['2014-08-28/09:00','2014-08-28/14:00'])
;    sorted_probes = ['thd','the','g13','tha','rbspb','g15']
;    panel_x_ratio = 0.8 ; inch per hour.
;    gmag_zlog = 0
;    gmag_range = [0,400]
;    mlat_range = [60,75]


    mlt_range = [-9,9]
    time_range = time_double(['2013-06-07/02:30','2013-06-07/07:00'])
    sorted_probes = ['tha','thd','the','g13','rbspa','rbspb','g15']
    panel_x_ratio = 0.8 ; inch per hour.
    gmag_zlog = 0
    gmag_range = [0,400]
    mlat_range = [60,80]


    gmag_var = 'thg_dbh'
    smooth_window = 120d*60
    long_time_range = time_range+[-1,1]*smooth_window
    if check_if_update(gmag_var, long_time_range) then themis_read_mag, long_time_range
    project = azim_df_load_project()
    azim_df_load_basic_data, project=project


    nprobe = n_elements(sorted_probes)

    line_linestyle = 2
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 1
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1
    rxy_range = [4.,30]
    root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
    ofn = join_path([root_dir,'fig_ewogram_of_dp_and_gmag_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hh'),'_')+'.pdf'])


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ; Get the sizes of the left panels.
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1,1]
    nypanel = n_elements(ypans)
    ypad = 0.4
    panel_xsize = abs(total(time_range*[-1,1]))/3600*panel_x_ratio
    margins = [10.,4,8,1]
    fig_xsize = panel_xsize+(total(margins[[0,2]]))*abs_xchsz
    fig_ysize = panel_ysize*total(ypans)+(total(margins[[1,3]])+ypad*(nypanel-1))*abs_ychsz

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
    poss = sgcalcpos(nypanel, ypans=ypans, ypad=ypad, margins=margins)


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
    x_constants = make_bins(time_range, 3600, /inner)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


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

    time_step = 10.
    times = make_bins(time_range+[0,-1]*time_step, time_step)
    ntime = n_elements(times)
    thetas = fltarr(ntime,nprobe)
    mlts = fltarr(ntime,nprobe)
    rxys = fltarr(ntime,nprobe)
    rsms = fltarr(ntime,3,nprobe)
    foreach probe, sorted_probes, probe_id do begin
        prefix = probe+'_'
        the_var = prefix+'scaled_theta'
        if tnames(the_var) eq '' then continue
        thetas[*,probe_id] = get_var_data(the_var, at=times)
        mlts[*,probe_id] = get_var_data(prefix+'pseudo_mlt', at=times)
        rsms[*,*,probe_id] = get_var_data(prefix+'r_sm', at=times)
        rxys[*,probe_id] = snorm(reform(rsms[*,0:1,probe_id]))
    endforeach

    foreach time, times, time_id do begin
        zzs = thetas[time_id,*]
        yys = mlts[time_id,*]
        rxy = rxys[time_id,*]
        rsm = rsms[time_id,*,*]

        index = where(finite(zzs),count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = azim_df_normalize_theta(zzs[index], zrange=tilt_range, ct=ct, /reverse_ct)
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Filter spatially.
        index = lazy_where(yys, '[]', yrange, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = lazy_where(rxy, '[]', rxy_range, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = where(check_if_in_magn(rsm, dynamic_pressure=pdyn) eq 1, count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Sort by tilt.
        index = reverse(sort(zzs))
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        for ii=0, count-1 do plots, time, yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach


    ; Add grid.
    y_constants = make_bins(mlt_range, 3, /inner)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle


    ; Add box.
    xtickformat = '(A1)'
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Draw color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks




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



    gmags = get_var_data(gmag_var, times=times, limits=lim)
    width = smooth_window/sdatarate(times)

    sites = lim.labels
    nsite = n_elements(sites)
    site_infos = themis_read_mag_metadata(sites=sites)
    mlons = fltarr(nsite)
    mlats = fltarr(nsite)
    foreach site_info, site_infos, site_id do begin
        mlons[site_id] = site_info.mlon
        mlats[site_id] = site_info.mlat

        gmags[*,site_id] -= smooth(gmags[*,site_id], width, /edge_truncate, /nan)
    endforeach
    gmag_var2 = gmag_var+'_smooth'
    if keyword_set(gmag_zlog) then begin
        store_data, gmag_var2, times, alog10(-gmags), limits=lim
    endif else begin
        store_data, gmag_var2, times, -gmags, limits=lim
    endelse


    gmags = get_var_data(gmag_var2, in=time_range, times=times)
    device, decomposed=0
    loadct, ct
    
    
    index = lazy_where(mlats, '[]', mlat_range, count=count)
    if count eq 0 then message, 'Invalid mlat_range ...'
    gmags = gmags[*,index]
    mlons = mlons[index]
    mlats = mlats[index]
    
    foreach time, times, time_id do begin
        yys = mlon2mlt(mlons, time)
        zzs = bytscl(gmags[time_id,*], min=zrange[0], max=zrange[1])

        index = lazy_where(yys, '[]', yrange, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]


        ; Sort by magnitude.
        index = sort(zzs)
        yys = yys[index]
        zzs = zzs[index]

        for ii=0, count-1 do plots, time, yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach
    device, decomposed=1

;    ; Sort sites.
;
;    gmags = get_var_data(gmag_var2, times=times, limits=lim, in=time_range)
;    avg_vals = fltarr(nsite)
;    foreach site, sites, site_id do begin
;        mlts = mlon2mlt(mlons[site_id], times)
;
;        index = lazy_where(mlts,'[]',yrange, count=nxx)
;        if nxx eq 0 then continue
;        xxs = times[index]
;        yys = mlts[index]
;        zzs = bytscl(gmags[index,site_id], min=zrange[0], max=zrange[1])
;        avg_vals[site_id] = mean(zzs)
;    endforeach
;    weight1 = avg_vals/255.         ; large current -> large value, plot last.
;    index = sort(weight1)
;
;    ; Do sort.
;    gmags = get_var_data(gmag_var2, times=times, limits=lim, in=time_range)
;    sites = sites[index]
;    gmags = gmags[*,index]
;    mlons = mlons[index]
;    mlats = mlats[index]
;
;    device, decomposed=0
;    loadct, ct
;    foreach site, sites, site_id do begin
;        if mlats[site_id] ge mlat_range[1] then continue
;        if mlats[site_id] le mlat_range[0] then continue
;        mlts = mlon2mlt(mlons[site_id], times)
;
;        index = lazy_where(mlts,'[]',yrange, count=nxx)
;        if nxx eq 0 then continue
;        xxs = times[index]
;        yys = mlts[index]
;        zzs = bytscl(gmags[index,site_id], max=zrange[1], min=zrange[0])
;        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize
;    endforeach
;    device, decomposed=1



    ; Add grid.
    y_constants = make_bins(mlt_range, 3, /inner)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle


    ; Add box.
    xtickformat = ''
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Draw color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks




    if keyword_set(test) then stop
    sgclose


end
