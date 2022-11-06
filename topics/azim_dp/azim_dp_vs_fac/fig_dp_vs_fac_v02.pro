;+
; Generate the plot for dp vs fac
; fig_dp_vs_fac_v02.pdf
; This version is for HGIO 2022 proposal.
;-


;---Settings.
test = 0

    plot_file = join_path([srootdir(),'fig_dp_vs_fac_v02.pdf'])
    if keyword_set(test) then plot_file = 0

    event_list = list()
    ; 2016-10-13 event.
    event_list.add, dictionary($
        'time_range', time_double(['2016-10-13/12:00','2016-10-13/13:30']), $
        'probes', ['rbspa','rbspb','g13','g14','g15','thd'], $
        'mlt_range', [-1,11], $
        'rxy_range', [4.,30], $
        'mlat_range', [60,80], $
        'ewo_mlat_range', [60,75], $
        'dp_info', dictionary($
            'times', time_double(['2016-10-13/12:20']), $
            'mlts', [1], $
            'omega_east', [-4.6], $
            'omega_west', [!values.f_nan] ) )
    ; 2013-06-07 event.
    event_list.add, dictionary($
        'time_range', time_double(['2013-06-07/04:25','2013-06-07/05:45']), $
        'probes', ['g13','rbspa','rbspb','g15','tha','thd','the'], $
        'snapshot_times', time_double('2013-06-07/'+['04:47','04:59','05:21','05:42']), $
        'snapshot_label', 'f-', $
        'mlt_range', [-8,8], $
        'rxy_range', [4.,30], $
        'mlat_range', [52,80], $
        'ewo_mlat_range', [55,75], $
        'dp_info', dictionary($
            'times', time_double(['2013-06-07/04:40']), $
            'mlts', [0], $
            'omega_east', [-3.5], $
            'omega_west', [2.0] ) )
    ; 2014-08-28 event.
    event_list.add, dictionary($
        'time_range', time_double(['2014-08-28/09:50','2014-08-28/11:10']), $
        'probes', ['tha','thd','g15','the','rbspb','g13'], $
        'mlt_range', [-2,7], $
        'rxy_range', [4.,30], $
        'mlat_range', [60,80], $
        'ewo_mlat_range', [60,75], $
        'dp_info', dictionary($
        'times', time_double(['2014-08-28/10:08']), $
        'mlts', [0], $
        'omega_east', [-2.1], $
        'omega_west', [3] ) )

;---Constant.
    secofday = 86400d


;---Calculate fig size.
    ; margins.
    margins = [9,4,6,1]
    ; panel sizes.
    ypans = [0.6,1,1,1]*0.9
    nypanel = n_elements(ypans)
    ypads = 0.4+dblarr(nypanel-1)
    fig_labels = $
        ['AE','J up','J down','Scaled '+tex2str('theta')]
    panel_vars = ['ae','j_up','j_down','scaled_theta']
    fig_letters = letters(nypanel)

    nxpanel = n_elements(event_list)
    xpans = dblarr(nxpanel)
    panel_x_ratio = 1.2 ; inch per hour.
    foreach event, event_list, event_id do begin
        time_range = event.time_range
        xpans[event_id] = abs(total(time_range*[-1,1]))/3600*panel_x_ratio
    endforeach
    xpads = 5+dblarr(nxpanel-1)
    ; figure sizes.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    fig_xsize = total(xpans)+$
        total(xpads)*abs_xchsz+$
        total(margins[[0,2]])*abs_xchsz
    fig_ysize = total(ypans)+$
        total(ypads)*abs_ychsz+$
        total(margins[[1,3]])*abs_ychsz

    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz
    all_poss = sgcalcpos(nypanel,nxpanel, $
        margins=margins, xpans=xpans, ypans=ypans, xpad=xpads, ypad=ypads)


;---J color bar.
    tpos = all_poss[*,nxpanel-1,2]
    j_cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*cb_width
    j_cbpos[3] = all_poss[3,1,1]
    j_ver_keo_color_table = 70
    j_up_ewo_color_table = 62
    j_down_ewo_color_table = 57



;---Adopted from fig_ewogram_of_dp_and_up_down_current.pro.
    for panel_xid=0, nxpanel-1 do begin
        poss = reform(all_poss[*,panel_xid,*])

    ;---Read theta, mlt, r_sm, j_up_ewo_var.
        event = event_list[panel_xid]
        time_range = event.time_range
        probes = event.probes
        nprobe = n_elements(probes)

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
        supermag_read_sme, time_range
        ae = get_var_data('ae', at=common_times)
        sme = get_var_data('sm_sme', at=common_times)
        store_data, 'ae_combo', common_times, [[ae],[sme]]
        add_setting, 'ae_combo', smart=1, dictionary($
            'display_type', 'stack', $
            'unit', 'nT', $
            'labels', ['AE','SME'], $
            'colors', sgcolor(['black','red']))



        ewo_mlat_range = event.ewo_mlat_range
        j_up_ewo_var = 'thg_j_up_ewo'
        if check_if_update(j_up_ewo_var, time_range) then begin
            themis_read_upward_current_ewo, time_range, mlt_range=mlt_range, mlat_range=ewo_mlat_range
        endif

        j_down_ewo_var = 'thg_j_down_ewo'
        if check_if_update(j_down_ewo_var, time_range) then begin
            themis_read_downward_current_ewo, time_range, mlt_range=mlt_range, mlat_range=ewo_mlat_range
        endif

    ;---Prepare to plot.
        full_ychsz = 0.7
        half_ychsz = 0.35
        label_size = 0.8
        constant_linestyle = 1

    ;---Common x-axis.
        xrange = time_range
        xticklen_chsz = -0.2
        yticklen_chsz = -0.4
        xminor = 3 ; hr.
        xstep = 30*60.
        xtickv = make_bins(xrange, xstep, /inner) ; make time line up.
        xticks = n_elements(xtickv)-1
        xtickn = strarr(xticks+1)
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


    ;---Settings for DP info.
        dp_info = event.dp_info
        dp_info_psym = 8
        dp_info_symsize = 1
        omega_east = dp_info.omega_east
        omega_west = dp_info.omega_west
        dp_times = dp_info.times
        dp_mlts = dp_info.mlts
        dp_info_linestyle = 2
        dp_info_color = sgcolor('dim_grey')

    ;---Panels in order.
        foreach panel_var, panel_vars, panel_yid do begin
            tpos = poss[*,panel_yid]
            xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
            yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
            xtickformat = '(A1)'
            if panel_yid eq nypanel-1 then xtickformat = ''
            ytickformat = ''

            if panel_var eq 'ae' then begin
                the_var = 'ae'
                yys = float(get_var_data(the_var, in=time_range, times=xxs, limits=lim))
                index = where(abs(yys) ge 1e5, count)
                if count ne 0 then yys[index] = !values.f_nan

                ystep = 100.
                yminor = 5
                yrange = minmax(make_bins(minmax(yys),ystep))
                ytickv = make_bins(yrange, ystep*yminor, /inner)
                yticks = n_elements(ytickv)-1
                ytitle = '(nT)'
                if panel_xid ne 0 then ytitle = ' '

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
            endif else if panel_var eq 'scaled_theta' then begin
                yrange = event.mlt_range
                ystep = 3
                yminor = ystep
                ytickv = make_bins(yrange,yminor, /inner)
                yticks = n_elements(ytickv)-1
                ytitle = 'MLT!C(hour)'
                if panel_xid ne 0 then ytitle = ' '
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
                rxy_range = event.rxy_range
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

                ; Add DP info.
                tmp = smkarthm(0,2*!dpi,20,'n')
                txs = cos(tmp)
                tys = sin(tmp)
                thick = 2
                usersym, txs, tys, thick=thick

                plots, dp_times,dp_mlts, psym=dp_info_psym, symsize=dp_info_symsize, color=dp_info_color
                foreach dp_time, dp_times, dp_id do begin
                    tmp1 = [dp_time,time_range[1]]
                    foreach the_omega, [omega_west[dp_id],omega_east[dp_id]] do begin
                        if ~finite(the_omega) then continue
                        slope = the_omega/60./15
                        tmp2 = dp_mlts[dp_id]-[0,slope*(tmp1[1]-tmp1[0])]
                        oplot, tmp1, tmp2, linestyle=dp_info_linestyle, color=dp_info_color

                        ; Find the position to print omega.
                        msg = strtrim(string(the_omega,format='(F5.1)'))+' deg/min'
                        ty = tpos[3]-ychsz*0.8
                        if the_omega gt 0 then ty = tpos[1]+ychsz*0.2
                        tmp = convert_coord(tpos[0],ty, normal=1, to_data=1)
                        ty = tmp[1]
                        tx = interp(tmp1, tmp2, ty)
                        tmp = convert_coord(tx,ty, data=1, to_normal=1)
                        tx = tmp[0]-xchsz*0.0
                        ty = tmp[1]
                        xyouts, tx,ty,/normal, alignment=1, charsize=label_size, msg
                    endforeach
                endforeach

                ; Add box.
                plot, xrange, yrange, position=tpos, $
                    xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
                    ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
                    /nodata, /noerase

                ; Draw color bar.
                if panel_xid eq nxpanel-1 then begin
                    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*cb_width
                    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks
                endif
            endif else if panel_var eq 'j_vert' then begin
                yrange = event.mlat_range
                ystep = 5
                ytickv = make_bins(yrange, ystep, /inner)
                yticks = n_elements(ytickv)-1
                ytickn = string(ytickv,format='(I0)')
                ytickn[1:*:2] = ' '
                yminor = ystep
                ytitle = 'MLat!C(deg)'
                if panel_xid ne 0 then ytitle = ' '
                nouttick = 1

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
                options, j_ver_keo_var1, 'xticklen', xticklen
                options, j_ver_keo_var1, 'yticklen', yticklen
                options, j_ver_keo_var1, 'color_table', j_ver_keo_color_table
                options, j_ver_keo_var1, 'zposition', j_cbpos
                options, j_ver_keo_var1, 'yrange', yrange
                options, j_ver_keo_var1, 'yticks', yticks
                options, j_ver_keo_var1, 'ytickv', ytickv
                options, j_ver_keo_var1, 'ytickname', ytickn
                options, j_ver_keo_var1, 'ytitle', ytitle
                options, j_ver_keo_var1, 'yminor', yminor
                tplot, j_ver_keo_var1, trange=time_range, position=tpos, noerase=1, nouttick=nouttick

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
                ty = tpos[3]-ychsz*1.8
                msg = 'avg ['+strjoin(string(mlt_range,format='(I0)'),',')+'] h'
                xyouts, tx,ty,/normal, msg, charsize=label_size

                ; Add DP info.
                foreach dp_time, dp_times, dp_id do begin
                    oplot, dp_time+[0,0],yrange, linestyle=dp_info_linestyle, color=dp_info_color
                endforeach
            endif else if panel_var eq 'j_up' then begin
                yrange = event.mlt_range
                ystep = 3
                ytickv = make_bins(yrange, ystep, /inner)
                yticks = n_elements(ytickv)-1
                ytickn = string(ytickv,format='(I0)')
                yminor = ystep
                ytitle = 'MLT (h)'
                if panel_xid ne 0 then ytitle = ' '
                nouttick = 1

                j_up_ewo_var = 'thg_j_up_ewo'
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
                tplot, j_up_ewo_var1, trange=time_range, position=tpos, noerase=1, nouttick=nouttick

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

                ; Add label.
                tx = tpos[0]+xchsz*0.5
                ty = tpos[3]-ychsz*1.8
                msg = 'avg ['+strjoin(string(ewo_mlat_range,format='(I0)'),',')+'] deg'
                xyouts, tx,ty,/normal, msg, charsize=label_size

                ; Add DP info.
                tmp = smkarthm(0,2*!dpi,20,'n')
                txs = cos(tmp)
                tys = sin(tmp)
                thick = 2
                usersym, txs, tys, thick=thick

                plots, dp_times,dp_mlts, psym=dp_info_psym, symsize=dp_info_symsize, color=dp_info_color
                foreach dp_time, dp_times, dp_id do begin
                    tmp1 = [dp_time,time_range[1]]
                    foreach the_omega, [omega_west[dp_id],omega_east[dp_id]] do begin
                        if ~finite(the_omega) then continue
                        slope = the_omega/60./15
                        tmp2 = dp_mlts[dp_id]-[0,slope*(tmp1[1]-tmp1[0])]
                        oplot, tmp1, tmp2, linestyle=dp_info_linestyle, color=dp_info_color
                    endforeach
                endforeach
            endif else if panel_var eq 'j_down' then begin
                yrange = event.mlt_range
                ystep = 3
                ytickv = make_bins(yrange, ystep, /inner)
                yticks = n_elements(ytickv)-1
                ytickn = string(ytickv,format='(I0)')
                yminor = ystep
                ytitle = 'MLT (h)'
                if panel_xid ne 0 then ytitle = ' '
                nouttick = 1

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
                tplot, j_down_ewo_var1, trange=time_range, position=tpos, noerase=1, nouttick=nouttick

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

                ; Add label.
                tx = tpos[0]+xchsz*0.5
                ty = tpos[3]-ychsz*1.8
                msg = 'avg ['+strjoin(string(ewo_mlat_range,format='(I0)'),',')+'] deg'
                xyouts, tx,ty,/normal, msg, charsize=label_size

                ; Add DP info.
                tmp = smkarthm(0,2*!dpi,20,'n')
                txs = cos(tmp)
                tys = sin(tmp)
                thick = 2
                usersym, txs, tys, thick=thick

                plots, dp_times,dp_mlts, psym=dp_info_psym, symsize=dp_info_symsize, color=dp_info_color
                foreach dp_time, dp_times, dp_id do begin
                    tmp1 = [dp_time,time_range[1]]
                    foreach the_omega, [omega_west[dp_id],omega_east[dp_id]] do begin
                        if ~finite(the_omega) then continue
                        slope = the_omega/60./15
                        tmp2 = dp_mlts[dp_id]-[0,slope*(tmp1[1]-tmp1[0])]
                        oplot, tmp1, tmp2, linestyle=dp_info_linestyle, color=dp_info_color
                    endforeach
                endforeach
            endif else if panel_var eq 'sme2d' then begin
                yrange = [-1,1]*12
                ystep = 6
                ytickv = make_bins(yrange, ystep, /inner)
                yticks = n_elements(ytickv)-1
                ytickn = string(ytickv,format='(I0)')
                yminor = ystep
                ytitle = 'MLT (h)'
                if panel_xid ne 0 then ytitle = ' '
                nouttick = 1
                
                sme_var = 'sm_sme2d'
                zlim, sme_var, 100,1000, 1
                options, sme_var, 'ztitle', 'SME (nT)'
                options, sme_var, 'zcharsize', 0.8
                options, sme_var, 'zposition', j_cbpos2
                options, sme_var, 'xticklen', xticklen
                options, sme_var, 'yticklen', yticklen
                options, sme_var, 'yrange', yrange
                options, sme_var, 'yticks', yticks
                options, sme_var, 'ytickv', ytickv
                options, sme_var, 'ytitle', ytitle
                options, sme_var, 'yminor', yminor
                tplot, sme_var, trange=time_range, position=tpos, noerase=1, nouttick=nouttick

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

                ; Add label.
                tx = tpos[0]+xchsz*0.5
                ty = tpos[3]-ychsz*1.8
                msg = 'avg ['+strjoin(string(ewo_mlat_range,format='(I0)'),',')+'] deg'
                xyouts, tx,ty,/normal, msg, charsize=label_size

                ; Add DP info.
                tmp = smkarthm(0,2*!dpi,20,'n')
                txs = cos(tmp)
                tys = sin(tmp)
                thick = 2
                usersym, txs, tys, thick=thick

                plots, dp_times,dp_mlts, psym=dp_info_psym, symsize=dp_info_symsize, color=dp_info_color
                foreach dp_time, dp_times, dp_id do begin
                    tmp1 = [dp_time,time_range[1]]
                    foreach the_omega, [omega_west[dp_id],omega_east[dp_id]] do begin
                        if ~finite(the_omega) then continue
                        slope = the_omega/60./15
                        tmp2 = dp_mlts[dp_id]-[0,slope*(tmp1[1]-tmp1[0])]
                        oplot, tmp1, tmp2, linestyle=dp_info_linestyle, color=dp_info_color
                    endforeach
                endforeach
                stop
            endif
        endforeach
        

    ;---Panel label.
        for panel_yid=0,nypanel-1 do begin
            tpos = poss[*,panel_yid]
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            msg = fig_letters[panel_yid]+'-'+string(panel_xid+1,format='(I0)')+'.'
            xyouts, tx,ty,normal=1, msg
            if panel_xid eq 0 then begin
                tx = xchsz*0.5
                ty = tpos[3]-ychsz*0.9
                msg = fig_letters[panel_yid]+'. '+fig_labels[panel_yid]
                xyouts, tx,ty,normal=1, msg
            endif
        endfor

    ;---Add xticks.
        tpos = poss[*,nypanel-1]
        xtickformat = ''
        ; Add axes.
        plot, xrange, yrange, position=tpos, $
            xstyle=9, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=5, $
            /nodata, /noerase
    endfor
    
;---Add color table for ewograms.
    zrange = [-50,50]
    ztitle = 'Vertical current (kA), red for upward'
    zstep = 20
    ztickv = make_bins(zrange, zstep, inner=1)
    sgcolorbar, 256-findgen(256), ct=j_ver_keo_color_table, $
        zrange=zrange, position=j_cbpos, ztitle=ztitle, ztickv=ztickv;, ztickn=ztickn, zticks=zticks


;---Done.
    if keyword_set(test) then stop
    sgclose
end
