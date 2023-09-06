;+
; Generate the summary plot for the 2014-08-28 event.
; fig_2014_0828_0950_v01.pdf
;-


;---Settings.
    time_range = time_double(['2014-08-28/09:50','2014-08-28/11:10'])
    probes = ['tha','thd','g15','the','rbspb','g13']
    mlt_range = [-2,7]
    rxy_range = [4.,30]
    mlat_range = [60,80]
    ewo_mlat_range=  [60,75]
test = 0

;---Settings for the start location and slopes (omega).
    the_line_style = 2
    the_line_color = sgcolor('dim_grey')
    the_start = dictionary($
        'times', time_double(['2014-08-28/10:08']), $
        'mlts', [0], $
        'omega1', [-2.1], $
        'omega2', [3] )


;---Adopted from fig_ewogram_of_dp_and_up_down_current.pro.
    duration = total(time_range*[-1,1])

    nprobe = n_elements(probes)
    plot_file = join_path([homedir(),'Dropbox','mypapers','dp_vs_fac','plot',$
        'fig_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmm'),'_')+'_v01.pdf'])
    if keyword_set(test) then plot_file = 0


;---Read theta, mlt, r_sm, j_up_ewo_var.
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

    j_up_ewo_var = 'thg_j_up_ewo'
    if check_if_update(j_up_ewo_var, time_range) then begin
        themis_read_upward_current_ewo, time_range, mlt_range=mlt_range, mlat_range=ewo_mlat_range
    endif

    j_down_ewo_var = 'thg_j_down_ewo'
    if check_if_update(j_down_ewo_var, time_range) then begin
        themis_read_downward_current_ewo, time_range, mlt_range=mlt_range, mlat_range=ewo_mlat_range
    endif


    ; Settings.
    probe_xys = [-1,1]
    line_linestyle = 2
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 1
    label_yshift = full_ychsz
    label_size = 0.8
    constant_linestyle = 1


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ; Get the sizes of the left panels.
    margins = [10.5,4,7,1]
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1,1,1,0.8]
    nypanel = n_elements(ypans)
    fig_labels = letters(nypanel)+'. '+$
        ['AE','DP','J up','J down',$
        'J vert']
    fig_labels = letters(nypanel)+'2. '+$
        ['AE','Scaled '+tex2str('theta'),'Upwd!Ccurrent','Downwd!Ccurrent',$
        'Vertical!Ccurrent']
    ypad = 0.4
    ; xsize of the main panel.
    panel_x_ratio = 1.5 ; inch per hour.
    panel_xsize = abs(total(time_range*[-1,1]))/3600*panel_x_ratio
    fig_xsize = panel_xsize+(total(margins[[0,2]]))*abs_xchsz
    fig_ysize = panel_ysize*total(ypans)+(total(margins[[1,3]])+ypad*(nypanel-1))*abs_ychsz

;---Common x-axis.
    xrange = time_range
    xticklen_chsz = -0.2
    yticklen_chsz = -0.4
    xminor = 3 ; hr.
    xstep = 1800.
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



;---Panel a. AE.
    the_var = 'ae'
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xtickformat = '(A1)'
    ystep = 500.
    yys = float(get_var_data(the_var, in=time_range, times=xxs))
    index = where(abs(yys) ge 1e5, count)
    if count ne 0 then yys[index] = !values.f_nan
    yrange = [100,1200]
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
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
    y_constants = ytickv
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

;    tx = tpos[0]+xchsz*0.5
;    ty = tpos[3]-ychsz*0.8
;    msg = 'Auroral electrojet (AE) index'
;    xyouts, tx,ty,/normal, msg, charsize=label_size


;---Panel b. UT-MLT long.
    tpos = poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    ytitle = 'MLT!C(hour)'
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
        index = where_pro(yys, '[]', yrange, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = where_pro(rxy, '[]', rxy_range, count=count)
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


;---Panel c, EWOgram of upward current.
    tpos = poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    yrange = mlt_range
    ystep = 3
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    ytitle = 'MLT!C(hour)'

    j_up_ewo_color_table = 62
    j_down_ewo_color_table = 57
    j_ver_keo_color_table = 70

    j_up_ewo_var1 = j_up_ewo_var+'1'
    copy_data, j_up_ewo_var, j_up_ewo_var1
    get_data, j_up_ewo_var1, xx,yy,zz
    yy *= -1e-3
    store_data, j_up_ewo_var1, xx,yy,zz
    options, j_up_ewo_var1, 'ztitle', 'Upward current (kA)'
    options, j_up_ewo_var1, 'zrange', [-50,0]
    options, j_up_ewo_var1, 'zcharsize', 0.8
    options, j_up_ewo_var1, 'zposition', cbpos
    options, j_up_ewo_var1, 'xticklen', xticklen
    options, j_up_ewo_var1, 'yticklen', yticklen
    options, j_up_ewo_var1, 'color_table', j_up_ewo_color_table
    options, j_up_ewo_var1, 'reverse_color_table', 1
    options, j_up_ewo_var1, 'no_color_scale', 1
    options, j_up_ewo_var1, 'yrange', yrange
    options, j_up_ewo_var1, 'yticks', yticks
    options, j_up_ewo_var1, 'ytickv', ytickv
    options, j_up_ewo_var1, 'ytitle', ytitle
    tplot, j_up_ewo_var1, trange=time_range, position=tpos, /noerase, /nouttick

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


;---Panel d, EWOgram of downward current.
    tpos = poss[*,3]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    yrange = mlt_range
    ystep = 3
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    ytitle = 'MLT!C(hour)'

    j_down_ewo_var1 = j_down_ewo_var+'1'
    copy_data, j_down_ewo_var, j_down_ewo_var1
    get_data, j_down_ewo_var1, xx,yy,zz
    yy *= 1e-3
    store_data, j_down_ewo_var1, xx,yy,zz
    options, j_down_ewo_var1, 'ztitle', 'Downward current (kA)'
    options, j_down_ewo_var1, 'zrange', [0,50]
    options, j_down_ewo_var1, 'zcharsize', 0.8
    options, j_down_ewo_var1, 'zposition', cbpos
    options, j_down_ewo_var1, 'xticklen', xticklen
    options, j_down_ewo_var1, 'yticklen', yticklen
    options, j_down_ewo_var1, 'color_table', j_down_ewo_color_table
    options, j_down_ewo_var1, 'no_color_scale', 1
    options, j_down_ewo_var1, 'yrange', yrange
    options, j_down_ewo_var1, 'yticks', yticks
    options, j_down_ewo_var1, 'ytickv', ytickv
    options, j_down_ewo_var1, 'ytitle', ytitle
    tplot, j_down_ewo_var1, trange=time_range, position=tpos, /noerase, /nouttick

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

    ; Add msg.
    foreach panel_id, [2,3] do begin
        tpos = poss[*,panel_id]
        plot, xrange, yrange, /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos

        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*0.8
        msg = 'EWOgram avg in!CMLat ['+strjoin(string(ewo_mlat_range,format='(I0)'),',')+'] deg'
        xyouts, tx,ty,/normal, msg, charsize=label_size
    endforeach


;---Add start location and slopes.
    start_psym = 8
    tmp = smkarthm(0,2*!dpi,20,'n')
    txs = cos(tmp)
    tys = sin(tmp)
    thick = 2
    usersym, txs, tys, thick=thick
    start_symsize = 1

    foreach panel_id, [1,2,3] do begin
        tpos = poss[*,panel_id]
        plot, xrange, yrange, /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos

        omega1 = the_start.omega1
        omega2 = the_start.omega2

        ; Add start location.
        txs = the_start.times
        tys = the_start.mlts
        plots, txs,tys, psym=start_psym, symsize=start_symsize, color=the_line_color
        nstart = n_elements(txs)
        for ii=0,nstart-1 do begin
            tmp1 = [txs[ii],time_range[1]]
            the_omega = omega1[ii]
            if finite(the_omega) then begin
                slope = the_omega/60./15
                tmp2 = tys[ii]-[0,slope*(tmp1[1]-tmp1[0])]
                oplot, tmp1, tmp2, linestyle=the_line_style, color=the_line_color

                if panel_id eq 1 then begin
                    ; Find the position to print omega.
                    ty = tpos[3]-ychsz*0.9
                    tmp = convert_coord(tpos[0],ty, /normal, /to_data)
                    tx = -(tmp[1]-tys[ii])/slope+txs[ii]
                    tmp = convert_coord(tx,tmp[1], /data, /to_normal)
                    tx = tmp[0]-xchsz*1
                    ty = tmp[1]
                    xyouts, tx,ty,/normal, alignment=1, charsize=label_size, $
                        strtrim(string(the_omega,format='(F5.1)'))+' deg/min'
                endif
            endif
            the_omega = omega2[ii]
            if finite(the_omega) then begin
                slope = the_omega/60./15
                tmp2 = tys[ii]-[0,slope*(tmp1[1]-tmp1[0])]
                oplot, tmp1, tmp2, linestyle=the_line_style, color=the_line_color

                if panel_id eq 1 then begin
                    ; Find the position to print omega.
                    ty = tpos[1]+ychsz*0.3
                    tmp = convert_coord(tpos[0],ty, /normal, /to_data)
                    tx = -(tmp[1]-tys[ii])/slope+txs[ii]
                    tmp = convert_coord(tx,tmp[1], /data, /to_normal)
                    tx = tmp[0]-xchsz*1
                    ty = tmp[1]
                    xyouts, tx,ty,/normal, alignment=1, charsize=label_size, $
                        strtrim(string(the_omega,format='(F5.1)'))+' deg/min'
                endif
            endif
        endfor
    endforeach



;---Panel e. KEOgram.
    tpos = poss[*,4]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    cbpos[3] = poss[3,2]
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    ytickn = string(ytickv,format='(I0)')
    ytickn[1:*:2] = ' '
    yminor = ystep
    ytitle = 'MLat!C(deg)'

    j_ver_keo_var = 'thg_j_ver_keo'
    if check_if_update(j_ver_keo_var,time_range) then $
        themis_read_j_ver_keo, time_range, mlt_range=mlt_range

    j_ver_keo_var1 = 'thg_j_ver_keo1'
    copy_data, j_ver_keo_var, j_ver_keo_var1
    get_data, j_ver_keo_var1, xx,yy,zz
    yy *= 1e-3
    store_data, j_ver_keo_var1, xx,yy,zz
    options, j_ver_keo_var1, 'ztitle', 'Vertical current (kA), negative value for upward current'
    options, j_ver_keo_var1, 'zrange', [-1,1]*50
    options, j_ver_keo_var1, 'zcharsize', 0.8
    options, j_ver_keo_var1, 'zposition', cbpos
    options, j_ver_keo_var1, 'xticklen', xticklen
    options, j_ver_keo_var1, 'yticklen', yticklen
    options, j_ver_keo_var1, 'color_table', j_ver_keo_color_table
    options, j_ver_keo_var1, 'zposition', cbpos
    options, j_ver_keo_var1, 'yrange', yrange
    options, j_ver_keo_var1, 'yticks', yticks
    options, j_ver_keo_var1, 'ytickv', ytickv
    options, j_ver_keo_var1, 'ytickname', ytickn
    options, j_ver_keo_var1, 'ytitle', ytitle
    options, j_ver_keo_var1, 'yminor', yminor
    tplot, j_ver_keo_var1, trange=time_range, position=tpos, /noerase

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

    ; Add msg.
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*0.8
    msg = 'KEOgram avg in!CMLT ['+strjoin(string(mlt_range,format='(I0)'),',')+'] hour'
    xyouts, tx,ty,/normal, msg, charsize=label_size



;---Done.
    if keyword_set(test) then stop
    sgclose
end
