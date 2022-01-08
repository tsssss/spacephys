;+
; Plot SC position, AE/Dst, and theta.
; 
; events.
; project=.
; dirname=.
;-

pro azim_df_subgroup_gen_diagnostic_plot, project=project, events, dirname=dirname, $
    skip_load_data=skip_load_data

test = 0
    secofhour = 3600.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(events) eq 0 then stop
    if n_elements(dirname) eq 0 then stop

    foreach event, events do begin
        event_time_range = event.time_range
        plot_time_range = event_time_range+[-1,1]*secofhour*0.5
        ;plot_time_range = minmax(make_bins(plot_time_range,600))
        time_range = plot_time_range
        all_probes = event.probes
        region_name = (strsplit(event.region,'%',/extract))[0]
        search_name = event.search_name

        case region_name of
            'post_midn': mlt_range = [0,1]*9.
            'pre_midn': mlt_range = [-1,0]*9.
            'around_midn': mlt_range = [-1,1]*3.
        endcase
        case search_name of
            'beyond_15Re': rxy_range = [4,30.]
            'within_15Re': rxy_range = [4,15.]
        endcase
        the_region = dictionary($
            'mlt_range', mlt_range, $
            'rxy_range', rxy_range)


        probe_infos = project.probe_infos
        xticklen_chsz = -0.15   ; in ychsz.
        yticklen_chsz = -0.30   ; in xchsz.
        full_ychsz = constant('full_ychsz')
        half_ychsz = full_ychsz*0.5
        lineskip = constant('lineskip')
        label_size = 0.7
        label_xshift = 10
        label_yshift = full_ychsz
        secofday = constant('secofday')
        bar_thick = keyword_set(test)? 0.5: 4

        rad = constant('rad')
        deg = constant('deg')


    ;---Load data.
        if ~keyword_set(skip_load_data) then azim_df_load_basic_data, project=project, scale_width=scale_width


    ;---Size of panel and figure.
        file_suffix = 'azim_df_event_'+strjoin(time_string(event_time_range,tformat='YYYY_MMDD_hhmm_ss'),'_')+'_v01.pdf'
        plot_file = join_path([project.plot_dir,'diagnostic_plot',dirname,file_suffix])
        sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
        sgclose, /wdelete
        ; 3 line panels: Dst/AE, MLT-UT, X-UT.
        line_panel_ysize = 1.5
        ypans = [0.5,1,1]*line_panel_ysize
        nypanel = n_elements(ypans)
        ypads = dblarr(nypanel-1)+0.5
        line_panel_aspect_ratio = 1     ; inch per hour.
        duration = total(time_range*[-1,1])/constant('secofhour')
        line_pos_xsize = line_panel_ysize*line_panel_aspect_ratio*duration
        line_pos_ysize = total(ypans)+total(ypads)*abs_ychsz
        ; 1 orbit panel: X-Y plane.
        orbit_pos_xsize = line_pos_ysize
        orbit_pos_ysize = line_pos_ysize
        ; Space b/w line and orbit panels.
        xpads = [15]
        ; Figure size.
        margins = [12,5,10,2]
        nxpanel = 2
        pos_xsize = line_pos_xsize+total(xpads)*abs_xchsz+orbit_pos_xsize
        pos_ysize = line_pos_ysize
        fig_xsize = pos_xsize+total(margins[[0,2]])*abs_xchsz
        fig_ysize = pos_ysize+total(margins[[1,3]])*abs_ychsz
        ; Panel positions.
        pos = [margins[0]*abs_xchsz/fig_xsize,margins[1]*abs_ychsz/fig_ysize, $
            1-margins[2]*abs_xchsz/fig_xsize,1-margins[3]*abs_ychsz/fig_ysize]
        line_pos = [pos[2]-line_pos_xsize/fig_xsize,pos[1],pos[2],pos[3]]
        orbit_pos = [pos[0],pos[1],pos[0]+orbit_pos_xsize/fig_xsize,pos[3]]
        pos_list = list()
        pos_list.add, orbit_pos
        pos_tops = dblarr(nypanel)+line_pos[3]
        for ii=1, nypanel-1 do pos_tops[ii] = pos_tops[ii-1]-(ypans[ii-1]+ypads[ii-1]*abs_ychsz)/fig_ysize
        for ii=0, nypanel-1 do begin
            pos_top = pos_tops[ii]
            tpos = [line_pos[0],pos_top-ypans[ii]/fig_ysize,line_pos[2],pos_top]
            pos_list.add, tpos
        endfor

        if n_elements(plot_file) eq 0 then plot_file = 0
        if keyword_set(test) then plot_file = test
        sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz
    ;    if keyword_set(test) then foreach tpos, pos_list do plot, [0,1], [0,1], /nodata, /noerase, position=tpos


    ;---Panel 1: SC orbit in XY plane.
        tpos = pos_list[0]

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        xrange = [20.,-40]
        xstep = 20
        xtickv = make_bins(xrange, xstep)
        xticks = n_elements(xtickv)-1
        xminor = 5
        xtitle = 'SM X (Re)'
        xticklen = xticklen

        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        yrange = (region_name eq 'pre_midn')? [40.,-20]: [20.,-40]
        ystep = 20
        ytickv = make_bins(yrange, ystep)
        yticks = n_elements(ytickv)-1
        yminor = 5
        ytitle = 'SM Y (Re)'
        yticklen = yticklen

        ; Set up coord.
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, /nodata, /noerase


        ; Add earth.
        nangle = 50
        angles = smkarthm(0,2*!dpi,nangle,'n')
        circle_x = cos(angles)
        circle_y = sin(angles)
        polyfill, circle_x<0, circle_y, color=sgcolor('silver')
        plots, circle_x, circle_y

        ; Add magnetopause.
        magn_test = fltarr(nangle,3)
        magn_test[*,0] = circle_x*30
        magn_test[*,1] = circle_y*30
        tmp = check_if_in_magn(magn_test, magn_pos=magn_pos)
        magn_pos = magn_pos[*,[0,2,1]]
        magn_pos = [magn_pos,magn_pos]
        magn_pos[nangle:nangle*2-1,1] *= -1
        oplot, magn_pos[*,0], magn_pos[*,1]

        ; Add ROI.
        mlt_range = the_region.mlt_range
        rxy_range = the_region.rxy_range
        foreach tmp, [5,10,20,40] do oplot, circle_x*tmp, circle_y*tmp, linestyle=1
        foreach tmp, rxy_range do oplot, circle_x*tmp, circle_y*tmp
        foreach tmp, mlt_range do oplot, -rxy_range*cos(tmp*15*rad), -rxy_range*sin(tmp*15*rad)

        ; Draw probe position for the plot_time_range.
        foreach df, event.df_list do begin
            probe = df.probe
            color = probe_infos[probe].color
            rsm = df.obs_r_sm
            plots, rsm[0], rsm[1], color=color, psym=1
        endforeach

        ; Add arrow if triad info exists.
        nvertex = 3
        dis_scale = 1.5
        vel_scale = 20
        hsize = keyword_set(test)? 6: 100
        arrow_solid = 0

        timing_keys = ['xcor','obs_time']
        timing_key = 'obs_time'
        if event.haskey('triad_list') then begin
            foreach triad, event.triad_list do begin
                rsm_center = triad.center_r_sm
                v_mag = triad['vmag_'+timing_key]
                v_hat = triad['vhat_'+timing_key]

                x0 = rsm_center[0]
                y0 = rsm_center[1]
                scale = dis_scale/vel_scale*v_mag
                x1 = x0+v_hat[0]*scale
                y1 = y0+v_hat[1]*scale
                arrow, x0,y0,x1,y1,/data, solid=arrow_solid, hsize=hsize
            endforeach
        endif




        ; Add notations.
        ty0 = (region_name eq 'post_midn')? tpos[3]: tpos[1]+ychsz*full_ychsz*(3+lineskip)
        tx = tpos[0]+xchsz*1
        step = 4

        ty = ty0-ychsz*full_ychsz*1
        xyouts, tx,ty,/normal, 'Probes of the candidate: ', charsize=label_size
        ttx = tx+15*xchsz*label_size
        foreach probe, all_probes do begin
            probe_info = probe_infos[probe]
            xyouts, ttx,ty,/normal, strupcase(probe_info.short_name), color=probe_info.color, charsize=label_size
            ttx += step*xchsz*label_size
        endforeach

        ty = ty0-ychsz*full_ychsz*2
        xyouts, tx,ty,/normal, 'Probes of the search: ', charsize=label_size
        ttx = tx+15*xchsz*label_size
        foreach probe, all_probes do begin
            probe_info = probe_infos[probe]
            xyouts, ttx,ty,/normal, strupcase(probe_info.short_name), color=probe_info.color, charsize=label_size
            ttx += step*xchsz*label_size
        endforeach

        ty = ty0-ychsz*full_ychsz*3
        duration = total(time_range*[-1,1])/constant('secofhour')
        msg = 'Orbit from '+strjoin(time_string(time_range,tformat='YYYY-MM-DD/hh:mm'), ' to ')+', duration is '+sgnum2str(duration,ndec=1)+' hour(s)'
        xyouts, tx,ty,/normal, msg, charsize=label_size

        ; Add axes.
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtitle=xtitle, xticklen=xticklen, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, yticklen=yticklen, $
            position=tpos, /noerase, /nodata

        ; Add labels.
        fig_label = 'a. XY!C    plane'
        tx = tpos[0]-label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, fig_label


    ;---The common x-axis setting.
        xrange = time_range
        xstep = constant('secofhour')
        xticks = floor(total(xrange*[-1,1])/xstep)
        xtickv = smkarthm(xrange[0], xstep, xticks+1, 'x0')
        xminor = 6
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


    ;---Panel b. Dst/AE.
        tpos = pos_list[1]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        the_var = 'dst'
        ystep = 50.
        constants = [-50,0]
        yys = get_var_data(the_var, in=time_range, times=xxs)
        yrange = minmax([yys,constants])
        ytickv = make_bins(yrange, ystep)
        yticks = n_elements(ytickv)-1
        yrange = minmax(ytickv)
        ytitle = '(nT)'

        ; Set up coord.
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, /nodata, /noerase

        ; Add data.
        oplot, xxs, yys
        foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

        ; Add axes.
        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=9, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase

        the_var = 'ae'
        ystep = 500.
        constants = [0,500]
        yys = get_var_data(the_var, in=time_range, times=xxs)
        yrange = minmax([yys,constants])
        ytickv = make_bins(yrange, ystep)
        yticks = n_elements(ytickv)-1
        yrange = minmax(ytickv)
        ytitle = '(nT)'

        ; Set up coord.
        ae_color = sgcolor('red')
        axis, yaxis=1, /save, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            color=ae_color

        ; Add data.
        oplot, xxs, yys, color=ae_color
        foreach constant, constants do oplot, xrange, constant+[0,0], color=ae_color, linestyle=1

        tx = tpos[0]-label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, 'b. Dst/'
        xyouts, tx,ty,/normal, '           AE', color=ae_color


    ;---Common z-range, colorbar, symbol.
        theta_var = 'theta'
        ztitle = 'Log!D10!N[Detrended tilt angle (deg)]'
        theta_range = [-1,1]*64
        zrange = theta_range
        ztickv = [-64,-16,-4,0,4,16,64]
        ztickn = string(ztickv,format='(I0)')
        ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
        zrange = alog10(abs(zrange))*[-1,1]
        zticks = n_elements(ztickv)-1
        spec_psym = 8
        spec_symsize = 0.3
        usersym, [1,1,-1,-1,1]*0.5, [-1,1,1,-1,-1]*1, /fill
        spec_ct = 70


    ;---Common x data.
        time_step = 30.
        xxs = make_bins(time_range, time_step)
        nxx = n_elements(xxs)

    ;---MLT-UT.
        tpos = pos_list[2]
        pos_var = 'pseudo_mlt'
        xtickformat = '(A1)'
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        ytitle = 'MLT (hr)'
        yminor = 3
        yrange = (region_name eq 'pre_midn')? mlt_range+[0,1]*yminor: mlt_range+[-1,0]*yminor
        ytickv = make_bins(yrange,yminor, /inner)
        yticks = n_elements(ytickv)-1
        constants = 0


        ; Set up coord.
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, /nodata, /noerase

        ; Plot data.
        foreach probe, all_probes do begin
            prefix = probe_infos[probe].prefix
            yys = get_var_data(prefix+pos_var, at=xxs)
            zzs = get_var_data(prefix+theta_var, at=xxs)
            mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
            zzs = azim_df_scale_theta(zzs, mlt) ; apply MLT scaling.
            zzs = azim_df_normalize_theta(zzs, zrange=theta_range, ct=spec_ct, /reverse_ct)

            ; Remove data outside ROI.
            index = lazy_where(mlt, '][', mlt_range, count=count)
            if count ne 0 then yys[index] = !values.f_nan
            rsm = get_var_data(prefix+'r_sm', at=xxs)
            rxy = snorm(rsm[*,0:1])
            index = lazy_where(rxy, '][', rxy_range, count=count)
            if count ne 0 then yys[index] = !values.f_nan
            index = where(finite(zzs,/nan), count)
            if count ne 0 then yys[index] = !values.f_nan

            ; Plot data.
            for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=spec_psym, symsize=spec_symsize
        endforeach

        ; Add DF.
        foreach df, event.df_list do begin
            ty = azim_df_calc_pseudo_mlt(df.obs_r_sm)
            tx = df.obs_time
            plots, tx,ty,/data, psym=1, symsize=0.5
        endforeach

        ; Add notations.
        foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

        ; Add axes.
        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase

        ; Add label.
        tx = tpos[0]-label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, 'c. UT-MLT'

        cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
        sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


    ;---R-UT.
        tpos = pos_list[3]
        xtickformat = ''
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        ytitle = 'R!Dxy!N (Re)'
        yrange = (search_name eq 'beyond_15Re')? [0,35]: [0,15]
        ystep = (search_name eq 'beyond_15Re')? 10: 5
        ytickv = make_bins(yrange, ystep, /inner)
        yticks = n_elements(ytickv)-1
        yminor = 5


        ; Set up coord.
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, /nodata, /noerase

        ; Plot data.
        foreach probe, all_probes do begin
            prefix = probe_infos[probe].prefix
            rsm = get_var_data(prefix+'r_sm', at=xxs)
            yys = snorm(rsm[*,0:1])
            zzs = get_var_data(prefix+theta_var, at=xxs)
            mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
            zzs = azim_df_normalize_theta(zzs, mlt, zrange=theta_range, ct=spec_ct, /reverse_ct)

            ; Remove data outside ROI.
            index = lazy_where(mlt, '][', mlt_range, count=count)
            if count ne 0 then yys[index] = !values.f_nan
            rsm = get_var_data(prefix+'r_sm', at=xxs)
            rxy = snorm(rsm[*,0:1])
            index = lazy_where(rxy, '][', rxy_range, count=count)
            if count ne 0 then yys[index] = !values.f_nan
            index = where(finite(zzs,/nan), count)
            if count ne 0 then yys[index] = !values.f_nan

            ; Plot data.
            for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=spec_psym, symsize=spec_symsize
        endforeach

        ; Add DF.
        foreach df, event.df_list do begin
            ty = df.obs_rxy
            tx = df.obs_time
            plots, tx,ty,/data, psym=1, symsize=0.5
        endforeach

        ; Add axes.
        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase

        ; Add label.
        tx = tpos[0]-label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, 'd. UT-R!Dxy!N'

        cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
        sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;    ;---X-UT.
;        tpos = pos_list[3]
;        pos_var = 'x_sm'
;        xtickformat = ''
;        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
;        ytitle = 'SM X (Re)'
;        yrange = (search_name eq 'beyond_15Re')? [-40,10]: [-20,10]
;        ystep = (search_name eq 'beyond_15Re')? 20: 10
;        ytickv = make_bins(yrange, ystep, /inner)
;        yticks = n_elements(ytickv)-1
;        yminor = 5
;
;
;        ; Set up coord.
;        plot, xrange, yrange, $
;            xstyle=5, xrange=xrange, $
;            ystyle=5, yrange=yrange, $
;            position=tpos, /nodata, /noerase
;
;        ; Plot data.
;        foreach probe, all_probes do begin
;            prefix = probe_infos[probe].prefix
;            rsm = get_var_data(prefix+'r_sm', at=xxs)
;            yys = rsm[*,0]
;            zzs = get_var_data(prefix+theta_var, at=xxs)
;            mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
;            zzs = azim_df_normalize_theta(zzs, mlt, zrange=theta_range, ct=spec_ct, /reverse_ct)
;
;            ; Remove data outside ROI.
;            index = lazy_where(mlt, '][', mlt_range, count=count)
;            if count ne 0 then yys[index] = !values.f_nan
;            rsm = get_var_data(prefix+'r_sm', at=xxs)
;            rxy = snorm(rsm[*,0:1])
;            index = lazy_where(rxy, '][', rxy_range, count=count)
;            if count ne 0 then yys[index] = !values.f_nan
;            index = where(finite(zzs,/nan), count)
;            if count ne 0 then yys[index] = !values.f_nan
;
;            ; Plot data.
;            for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=spec_psym, symsize=spec_symsize
;        endforeach
;
;        ; Add DF.
;        foreach df, event.df_list do begin
;            ty = df.obs_r_sm[0]
;            tx = df.obs_time
;            plots, tx,ty,/data, psym=1, symsize=0.5
;        endforeach
;
;        ; Add axes.
;        plot, xrange, yrange, position=tpos, $
;            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
;            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
;            /nodata, /noerase
;
;        ; Add label.
;        tx = tpos[0]-label_xshift*xchsz
;        ty = tpos[3]-label_yshift*ychsz
;        xyouts, tx,ty,/normal, 'd. UT-X!USM'
;
;        cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
;        sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


        if keyword_set(test) then stop
        sgclose
    endforeach

end




if n_elements(project) eq 0 then project = azim_df_load_project()

azim_df_subgroup_gen_diagnostic_plot, events, dirname=dirname, project=project

end