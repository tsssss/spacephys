;+
; Plot AE, ewograms of DP, gmag, upward current.
;-

pro fig_ewogram_of_dp_gmag_and_current, time_range, probes=probes, $
    filename=plot_file, test=test, $
    mlt_range=mlt_range, rxy_range=rxy_range, mlat_range=mlat_range, $
    fix_x_ratio=fix_x_ratio, panel_x_ratio=panel_x_ratio

    if n_elements(time_range) ne 2 then message, 'Invalid time_range ...'
    duration = total(time_range*[-1,1])
    if duration le 600 then message, 'time_range too short ...'

    if n_elements(mlt_range) ne 2 then mlt_range = [-9,9]
    if n_elements(rxy_range) ne 2 then rxy_range = [4.,30]
    if n_elements(mlat_range) ne 2 then mlat_range = [60,80]

    nprobe = n_elements(probes)
    if nprobe eq 0 then message, 'Invalid probes ...'
    if n_elements(plot_file) eq 0 then plot_file = join_path([homedir(),$
        'fig_ewogram_of_dp_gmag_and_current_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hh'),'_')+'_v01.pdf'])
    if keyword_set(test) then plot_file = 0

;---Read theta, mlt, r_sm, j_ewo_var, gmag_var.
    time_step = 10.
    common_times = make_bins(time_range+[0,-1]*time_step, time_step)
    foreach probe, probes do begin
        prefix = probe+'_'

        theta_var = prefix+'theta'
        if check_if_update(theta_var, time_range) then begin
            azim_dp_read_theta, time_range, probe=probe
            if tnames(theta_var) eq '' then continue
            interp_time, theta_var, common_times
        endif

        mlt_var = prefix+'mlt'
        if check_if_update(mlt_var, time_range) then begin
            azim_dp_read_mlt, time_range, probe=probe
            interp_time, mlt_var, common_times
        endif

        r_sm_var = prefix+'r_sm'
        if check_if_update(r_sm_var, time_range) then begin
            azim_dp_read_orbit, time_range, probe=probe
            interp_time, r_sm_var, common_times
        endif

        scaled_theta_var = prefix+'scaled_theta'
        if check_if_update(scaled_theta_var, time_range) then begin
            theta = get_var_data(theta_var)
            mlt = get_var_data(mlt_var)
            scaled_theta = azim_dp_scale_theta(theta, mlt, width=scale_width)
            store_data, scaled_theta_var, common_times, scaled_theta
        endif
    endforeach

    if check_if_update('ae', time_range) then omni_read_index, time_range

    j_ewo_var = 'thg_j_ver_up_ewo'
    if check_if_update(j_ewo_var, time_range) then begin
        themis_read_upward_current_ewo, time_range, mlt_range=mlt_range, mlat_range=mlat_range
    endif


    gmag_var = 'thg_dbh'
    smooth_window = 120d*60
    long_time_range = time_range+[-1,1]*smooth_window
    gmag_var2 = gmag_var+'_smooth'
    if check_if_update(gmag_var2, long_time_range) then begin
        if check_if_update(gmag_var, long_time_range) then begin
            themis_read_mag, long_time_range
        endif
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
        store_data, gmag_var2, times, -gmags, limits=lim
        options, gmag_var2, 'mlons', mlons
        options, gmag_var2, 'mlats', mlats
        options, gmag_var2, 'sites', sites
    endif


;---Settings.
    probe_xys = [-1,1]
    line_linestyle = 2
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 1
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1
    panel_xsize_max = 10
    panel_xsize_min = 4


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ; Get the sizes of the left panels.
    margins = [10.,4,8,1]
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1,1,1]
    nypanel = n_elements(ypans)
    fig_labels = letters(nypanel)+'. '+['AE','Dipolari-!C     zation','Gmag','Upward!C    current']
    ypad = 0.4
    ; xsize of the main panel.
    if n_elements(panel_x_ratio) eq 0 then panel_x_ratio = 1 ; inch per hour.
    panel_xsize = abs(total(time_range*[-1,1]))/3600*panel_x_ratio
    if ~keyword_set(fix_x_ratio) then begin
        panel_xsize <= panel_xsize_max
        panel_xsize >= panel_xsize_min
    endif
    fig_xsize = panel_xsize+(total(margins[[0,2]]))*abs_xchsz
    fig_ysize = panel_ysize*total(ypans)+(total(margins[[1,3]])+ypad*(nypanel-1))*abs_ychsz

;---Common x-axis.
    xrange = time_range
    xticklen_chsz = -0.2
    yticklen_chsz = -0.4
    xminor = 6 ; hr.
    if duration ge 7200. then begin
        xstep = 3600.
    endif else if duration ge 3600 then begin
        xstep = 1800.
    endif else if duration ge 1800 then begin
        xstep = 600.
    endif else step = 60.
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


    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz
    poss = sgcalcpos(nypanel, ypans=ypans, ypad=ypad, margins=margins)


    ; Label.
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
    yys = float(get_var_data(the_var, in=time_range, times=xxs))
    index = where(abs(yys) ge 1e4, count)
    if count ne 0 then yys[index] = !values.f_nan
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
    x_constants = make_bins(time_range, xstep, /inner)
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
    ystep = 3
    yminor = ystep
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


    ntime = n_elements(common_times)
    thetas = fltarr(ntime,nprobe)
    mlts = fltarr(ntime,nprobe)
    rxys = fltarr(ntime,nprobe)
    rsms = fltarr(ntime,3,nprobe)
    foreach probe, probes, probe_id do begin
        prefix = probe+'_'
        the_var = prefix+'scaled_theta'
        if tnames(the_var) eq '' then continue
        thetas[*,probe_id] = get_var_data(the_var)
        mlts[*,probe_id] = get_var_data(prefix+'mlt')
        rsms[*,*,probe_id] = get_var_data(prefix+'r_sm')
        rxys[*,probe_id] = snorm(reform(rsms[*,0:1,probe_id]))
    endforeach


    foreach time, common_times, time_id do begin
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
    y_constants = make_bins(yrange, ystep, /inner)
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


;---Panel b, EWOgram of gmag.
    tpos = poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    yrange = mlt_range
    ystep = 3
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1

    ztitle = 'dB south'+' (nT)'
    gmag_range = [0,400]
    zrange = gmag_range
    zstep = abs(total(gmag_range*[-1,1])/2)
    ztickv = make_bins(zrange,zstep,/inner)
    ztickn = string(ztickv,format='(I0)')
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

    ; Load data.
    gmags = get_var_data(gmag_var2, in=time_range, times=times, limits=lim)
    mlons = lim.mlons
    mlats = lim.mlats
    index = lazy_where(mlats, '[]', mlat_range, count=count)
    if count eq 0 then message, 'Invalid mlat_range ...'
    gmags = gmags[*,index]
    mlons = mlons[index]
    mlats = mlats[index]

    device, decomposed=0
    loadct, ct
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

    ; Add grid.
    y_constants = make_bins(yrange, ystep, /inner)
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
    sgcolorbar, findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks



;---Panel c, EWOgram of current.
    tpos = poss[*,3]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    yrange = mlt_range
    ystep = 3
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1


    j_ewo_var1 = j_ewo_var+'1'
    copy_data, j_ewo_var, j_ewo_var1
    get_data, j_ewo_var1, xx,yy,zz
    yy *= 1e-3
    store_data, j_ewo_var1, xx,yy,zz
    options, j_ewo_var1, 'ztitle', 'Upward current (kA)'
    options, j_ewo_var1, 'zrange', [0,50]
    options, j_ewo_var1, 'zcharsize', 0.8
    options, j_ewo_var1, 'zposition', cbpos
    j_ewo_color_table = 62
    options, j_ewo_var1, 'xticklen', xticklen
    options, j_ewo_var1, 'yticklen', yticklen
    options, j_ewo_var1, 'color_table', j_ewo_color_table
    options, j_ewo_var1, 'yrange', yrange
    options, j_ewo_var1, 'yticks', yticks
    options, j_ewo_var1, 'ytickv', ytickv
    tplot, j_ewo_var1, trange=time_range, position=tpos, /noerase

    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, noerase=1, nodata=1

    ; Add grid.
    grid_color = sgcolor('gray')
    grid_linestyle = 1
    y_constants = make_bins(yrange, ystep, /inner)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle



    if keyword_set(test) then stop
    sgclose
end




time_range = time_double(['2013-06-07/02:30','2013-06-07/07:00'])
probes = ['tha','thd','the','g13','rbspa','rbspb','g15']
mlt_range = [-9,7]
mlat_range = [60,75]

;time_range = time_double(['2014-08-28/09:50','2014-08-28/11:10'])
;probes = ['tha','thd','the','g13','rbspb','g15']
;mlt_range = [-3,7]
;mlat_range = [60,75]
;
time_range = time_double(['2008-02-29/07:50','2008-02-29/09:00'])
probes = ['th'+['a','c','d','e'],'g'+['11','12']]
mlt_range = [-1,1]*6
mlat_range = [60,70]


;; No current data.
;time_range = time_double(['2008-03-20/11:30','2008-03-20/13:00'])
;probes = ['th'+['a','b','c','d','e']]
;mlt_range = [-1,1]*6
;mlat_range = [60,70]

;
;; Current not obvious.
;time_range = time_double(['2014-12-26/00:30','2014-12-26/02:00'])
;probes = ['g13','g15','tha','the']
;mlt_range = [-12,0]
;mlat_range = [60,80]
;
;; Current does not go along with DP, but start at the same time.
;time_range = time_double(['2016-06-15/03:30','2016-06-15/06:00'])
;probes = ['g13','g14','g15','tha','mms1']
;mlt_range = [-9,6]
;mlat_range = [60,75]
;
; Good example.
time_range = time_double(['2016-10-13/12:00','2016-10-13/13:30'])
probes = ['rbspa','rbspb','g13','g14','g15','thd']
mlt_range = [0,9]
mlat_range = [60,80]
;
; Interesting example.
;time_range = time_double(['2017-03-31/02:00','2017-03-31/03:30'])
;probes = ['g13','g15','tha','thd','the','rbspa','rbspb']
;mlt_range = [-9,3]
;mlat_range = [60,65]

; Runov's example. No significant FAC.
time_range = time_double(['2009-02-27/07:00','2009-02-27/09:00'])
probes = 'th'+letters('e')
mlt_range = [-1,1]*6
mlat_range = [60,80]


; Mary Hudson's event.
time_range = time_double(['2017-03-21/07:00','2017-03-21/09:00'])
probes = ['tha','thd','the','g13','g15']
mlt_range = [-6,9]
mlat_range = [70,80]
test = 1

time_range = time_double(['2019-08-05/23:00','2019-08-06/00:30'])
probes = ['tha','thd','the','g14','g16']
mlt_range = [-9,3]
mlat_range = [55,70]
test = 1



fig_ewogram_of_dp_gmag_and_current, time_range, probes=probes, test=test, $
    mlt_range=mlt_range, mlat_range=mlat_range
end
