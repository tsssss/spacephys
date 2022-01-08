pro azim_df_gen_diagnostic_plot_search_event_plot, event, project=project, filename=plot_file

test = 0

;---event properties, general settings.
    time_range = event.time_range
    if event.haskey('df_probes') then all_probes = event.df_probes else all_probes = event.probes
    region_name = event.region
    the_region = project.search_candidate.search_roi.regions[region_name]
    search_name = event.search_type
    search_settings = project.search_settings
    search_setting = (search_settings[0].name eq search_name)? search_settings[0]: search_settings[1]
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


;---Size of panel and figure.
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
    rxy_range = project.search_candidate.search_roi.rxy_range
    foreach tmp, [5,10,20,40] do oplot, circle_x*tmp, circle_y*tmp, linestyle=1
    foreach tmp, rxy_range do oplot, circle_x*tmp, circle_y*tmp
    foreach tmp, mlt_range do oplot, -rxy_range*cos(tmp*15*rad), -rxy_range*sin(tmp*15*rad)

    ; Draw probe position for the full_time_range.
    foreach probe, all_probes do begin
        prefix = probe_infos[probe].prefix
        color = probe_infos[probe].color
        rsm = get_var_data(prefix+'r_sm', in=time_range)
        xxs = rsm[*,0]
        yys = rsm[*,1]
        oplot, xxs, yys, color=color
    endforeach


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
    foreach probe, search_setting.probes do begin
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
    foreach probe, all_probes do begin
        prefix = probe_infos[probe].prefix
        foreach df_info, event.df_list do begin
            if df_info.probe ne probe then continue
            tx = df_info.arrival_time
            ty = get_var_data(prefix+pos_var, at=tx)
            plots, tx,ty,/data, psym=1, symsize=0.5
        endforeach
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



;---X-UT.
    tpos = pos_list[3]
    pos_var = 'x_sm'
    xtickformat = ''
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    ytitle = 'SM X (Re)'
    yrange = (search_name eq 'beyond_15Re')? [-40,10]: [-20,10]
    ystep = (search_name eq 'beyond_15Re')? 20: 10
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
        yys = rsm[*,0]
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
    foreach probe, all_probes do begin
        prefix = probe_infos[probe].prefix
        foreach df_info, event.df_list do begin
            if df_info.probe ne probe then continue
            tx = df_info.arrival_time
            ty = get_var_data(prefix+pos_var, at=tx)
            plots, tx,ty,/data, psym=1, symsize=0.5
        endforeach
    endforeach

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Add label.
    tx = tpos[0]-label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'd. UT-X!USM'

    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks

    if keyword_set(test) then stop
    sgclose


end




pro azim_df_gen_diagnostic_plot_search_event_search_df_group, project=project, $
    test_time=test_time, region_names=region_names, search_names=search_names

;---Handle which events are plotted.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_search_event(project=project)
    if keyword_set(test_time) then begin
        event_ids = list()
        foreach event, events do begin
            if product(event.time_range-test_time) gt 0 then continue else event_ids.add, event.id
        endforeach
        event_ids = event_ids.toarray()-1   ; id starts from 1.
        events = events[event_ids]
    endif

    ; Filter by region_names and search_names.
    regions = project.search_candidate.search_roi.regions
    if n_elements(region_names) eq 0 then begin
        region_names = list()
        foreach region, regions do region_names.add, region.name
    endif
    if ~isa(region_names,'list') then region_names = list(region_names, /extract)

    search_settings = project.search_settings
    if n_elements(search_names) eq 0 then begin
        search_names = list()
        foreach search_setting, search_settings do search_names.add, search_setting.name
    endif
    if ~isa(search_names,'list') then search_names = list(search_names, /extract)

    event_ids = list()
    foreach event, events do begin
        if region_names.where(event.region) eq !null then continue
        if search_names.where(event.search_type) eq !null then continue
        event_ids.add, event.id
    endforeach
    event_ids = event_ids.toarray()-1
    events = events[event_ids]

;---Loop through event to generate plots.
    probe_infos = project.probe_infos
    foreach search_setting, search_settings do begin
        if search_names.where(search_setting.name) eq !null then continue
        ; Load data for one time for all events in this search.
        time_range = search_setting.time_range
        time_step = project.time_step
        common_times = make_bins(time_range, time_step)
        pdyn = project.search_candidate.search_roi.pdyn
        rxy_range = project.search_candidate.search_roi.rxy_range
        mlt_range = project.search_candidate.search_roi.mlt_range
        foreach probe, search_setting.probes do begin
            lprmsg, 'Processing '+strupcase(probe)+' ...'
            prefix = probe_infos[probe].prefix
            r_sm_var = prefix+'r_sm'
            if check_if_update(r_sm_var, time_range) then azim_df_read_data, 'r_sm', time_range=time_range, probe=probe, project=project
            mlt_var = prefix+'pseudo_mlt'
            if check_if_update(mlt_var, time_range) then begin
                get_data, r_sm_var, times, rsm
                mlt = azim_df_calc_pseudo_mlt(rsm)
                store_data, mlt_var, times, mlt
            endif
            x_sm_var = prefix+'x_sm'
            if check_if_update(x_sm_var, time_range) then begin
                get_data, r_sm_var, times, rsm
                store_data, x_sm_var, times, rsm[*,0]
            endif

            theta_var = prefix+'theta'
            if check_if_update(theta_var, time_range) then azim_df_read_data, 'theta', time_range=time_range, probe=probe, project=project

            ; Remove data outside the general ROI.
            roi_flags = bytarr(n_elements(common_times))+1

            ; Magnetopause.
            r_gsm_var = prefix+'r_gsm'
            if check_if_update(r_gsm_var, time_range) then azim_df_read_data, 'r_gsm', time_range=time_range, probe=probe, project=project
            r_gsm = get_var_data(r_gsm_var, at=common_times)
            index = where(check_if_in_magn(r_gsm, dynamic_pressure=pdyn) eq 0, count)
            if count ne 0 then roi_flags[index] = 0
            ; rxy.
            r_sm = get_var_data(r_sm_var, at=common_times)
            rxy = snorm(r_sm[*,0:1])
            index = lazy_where(rxy, '][', rxy_range, count=count)
            if count ne 0 then roi_flags[index] = 0
            ; mlt.
            mlt = azim_df_calc_pseudo_mlt(r_sm)
            index = lazy_where(mlt, '][', mlt_range, count=count)
            if count ne 0 then roi_flags[index] = 0

            index = where(roi_flags eq 0, count)
            if count ne 0 then begin
                get_data, prefix+'theta', common_times, theta
                theta[index] = !values.f_nan
                store_data, prefix+'theta', common_times, theta
            endif
        endforeach
        dst_var = 'dst'
        if check_if_update(dst_var, time_range) then omni_read_index, time_range

        search_name = search_setting.name
        foreach region, regions do begin
            region_name = region.name
            if region_names.where(region_name) eq !null then continue
            foreach event, events do begin
                if event.search_type ne search_name then continue
                if event.region ne region_name then continue
                ; Extend time range to at least one hour.
                event.time_range += [-1,1]*0.5*constant('secofhour')
                plot_file = join_path([project.plot_dir,'diagnostic_plot','search_event','search_df_group',project.name+'_event_'+strjoin(time_string(event.time_range,tformat='YYYY_MMDD_hhmm'),'_')+'_v01.pdf'])
                azim_df_gen_diagnostic_plot_search_event_plot, event, project=project, filename=plot_file
            endforeach
        endforeach
    endforeach

end




pro azim_df_gen_diagnostic_plot_search_event_search_large_df, project=project, $
    test_time=test_time, region_names=region_names, search_names=search_names

;---Handle which events are plotted.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_search_event(project=project, /get_large_df_events)
    if keyword_set(test_time) then begin
        event_ids = list()
        foreach event, events do begin
            if product(event.time_range-test_time) gt 0 then continue else event_ids.add, event.id
        endforeach
        event_ids = event_ids.toarray()-1   ; id starts from 1.
        events = events[event_ids]
    endif

    ; Filter by region_names and search_names.
    regions = project.search_candidate.search_roi.regions
    if n_elements(region_names) eq 0 then begin
        region_names = list()
        foreach region, regions do region_names.add, region.name
    endif
    if ~isa(region_names,'list') then region_names = list(region_names, /extract)

    search_settings = project.search_settings
    if n_elements(search_names) eq 0 then begin
        search_names = list()
        foreach search_setting, search_settings do search_names.add, search_setting.name
    endif
    if ~isa(search_names,'list') then search_names = list(search_names, /extract)

    event_ids = list()
    foreach event, events do begin
        if region_names.where(event.region) eq !null then continue
        if search_names.where(event.search_type) eq !null then continue
        event_ids.add, event.id
    endforeach
    event_ids = event_ids.toarray()-1
    events = events[event_ids]

    ; Further filter by scale down the # plots.
    nevent = n_elements(events)
    expected_plot_number = 50
    if nevent gt expected_plot_number then begin
        step = ceil(nevent/expected_plot_number)
        event_ids = make_bins([0,nevent], step, /inner)
        events = events[event_ids]
    endif


;---Loop through event to generate plots.
    probe_infos = project.probe_infos
    foreach search_setting, search_settings do begin
        if search_names.where(search_setting.name) eq !null then continue
        ; Load data for one time for all events in this search.
        time_range = search_setting.time_range
        time_step = project.time_step
        common_times = make_bins(time_range, time_step)
        pdyn = project.search_candidate.search_roi.pdyn
        rxy_range = project.search_candidate.search_roi.rxy_range
        mlt_range = project.search_candidate.search_roi.mlt_range
        foreach probe, search_setting.probes do begin
            lprmsg, 'Processing '+strupcase(probe)+' ...'
            prefix = probe_infos[probe].prefix
            r_sm_var = prefix+'r_sm'
            if check_if_update(r_sm_var, time_range) then azim_df_read_data, 'r_sm', time_range=time_range, probe=probe, project=project
            mlt_var = prefix+'pseudo_mlt'
            if check_if_update(mlt_var, time_range) then begin
                get_data, r_sm_var, times, rsm
                mlt = azim_df_calc_pseudo_mlt(rsm)
                store_data, mlt_var, times, mlt
            endif
            x_sm_var = prefix+'x_sm'
            if check_if_update(x_sm_var, time_range) then begin
                get_data, r_sm_var, times, rsm
                store_data, x_sm_var, times, rsm[*,0]
            endif

            theta_var = prefix+'theta'
            if check_if_update(theta_var, time_range) then azim_df_read_data, 'theta', time_range=time_range, probe=probe, project=project

            ; Remove data outside the general ROI.
            roi_flags = bytarr(n_elements(common_times))+1

            ; Magnetopause.
            r_gsm_var = prefix+'r_gsm'
            if check_if_update(r_gsm_var, time_range) then azim_df_read_data, 'r_gsm', time_range=time_range, probe=probe, project=project
            r_gsm = get_var_data(r_gsm_var, at=common_times)
            index = where(check_if_in_magn(r_gsm, dynamic_pressure=pdyn) eq 0, count)
            if count ne 0 then roi_flags[index] = 0
            ; rxy.
            r_sm = get_var_data(r_sm_var, at=common_times)
            rxy = snorm(r_sm[*,0:1])
            index = lazy_where(rxy, '][', rxy_range, count=count)
            if count ne 0 then roi_flags[index] = 0
            ; mlt.
            mlt = azim_df_calc_pseudo_mlt(r_sm)
            index = lazy_where(mlt, '][', mlt_range, count=count)
            if count ne 0 then roi_flags[index] = 0

            index = where(roi_flags eq 0, count)
            if count ne 0 then begin
                get_data, prefix+'theta', common_times, theta
                theta[index] = !values.f_nan
                store_data, prefix+'theta', common_times, theta
            endif
        endforeach
        dst_var = 'dst'
        if check_if_update(dst_var, time_range) then omni_read_index, time_range

        search_name = search_setting.name
        foreach region, regions do begin
            region_name = region.name
            if region_names.where(region_name) eq !null then continue
            foreach event, events do begin
                if event.search_type ne search_name then continue
                if event.region ne region_name then continue
                plot_file = join_path([project.plot_dir,'diagnostic_plot','search_event','search_large_df',project.name+'_event_'+strjoin(time_string(event.time_range,tformat='YYYY_MMDD_hhmm'),'_')+'_v01.pdf'])
                azim_df_gen_diagnostic_plot_search_event_plot, event, project=project, filename=plot_file
            endforeach
        endforeach
    endforeach

end

pro azim_df_gen_diagnostic_plot_search_candidate_plot, candidate, project=project, filename=plot_file

test = 0

;---candidate properties, general settings.
    time_range = candidate.time_range
    all_probes = candidate.all_probes
    region_name = candidate.region
    nsection = candidate.nsection
    the_region = project.search_candidate.search_roi.regions[region_name]
    search_name = candidate.search_type
    search_settings = project.search_settings
    search_setting = (search_settings[0].name eq search_name)? search_settings[0]: search_settings[1]
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


;---Size of panel and figure.
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
    rxy_range = project.search_candidate.search_roi.rxy_range
    foreach tmp, [5,10,20,40] do oplot, circle_x*tmp, circle_y*tmp, linestyle=1
    foreach tmp, rxy_range do oplot, circle_x*tmp, circle_y*tmp
    foreach tmp, mlt_range do oplot, -rxy_range*cos(tmp*15*rad), -rxy_range*sin(tmp*15*rad)

    ; Draw probe position for the full_time_range.
    foreach probe, all_probes do begin
        prefix = probe_infos[probe].prefix
        color = probe_infos[probe].color
        rsm = get_var_data(prefix+'r_sm', in=time_range)
        xxs = rsm[*,0]
        yys = rsm[*,1]
        oplot, xxs, yys, color=color
    endforeach

    ; Highlight the sections and probes within ROI.
    for ii=0, nsection-1 do begin
        section_time_range = time_range_list[ii]
        section_probes = probe_list[ii]
        foreach probe, section_probes do begin
            prefix = probe_infos[probe].prefix
            color = probe_infos[probe].color
            rsm = get_var_data(prefix+'r_sm', in=time_range)
            xxs = rsm[*,0]
            yys = rsm[*,1]
            index = where(section_probes eq probe, count)
            oplot, xxs, yys, color=color, thick=bar_thick
        endforeach
    endfor

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
    foreach probe, search_setting.probes do begin
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



;---X-UT.
    tpos = pos_list[3]
    pos_var = 'x_sm'
    xtickformat = ''
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    ytitle = 'SM X (Re)'
    yrange = (search_name eq 'beyond_15Re')? [-40,10]: [-20,10]
    ystep = (search_name eq 'beyond_15Re')? 20: 10
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
        yys = rsm[*,0]
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

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Add label.
    tx = tpos[0]-label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'd. UT-X!USM'

    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks

    if keyword_set(test) then stop
    sgclose


end

pro azim_df_gen_diagnostic_plot_search_candidate, project=project, $
    test_time=test_time, region_names=region_names, search_names=search_names
    ; Internal, plot sc location, ROI; plot Dst/AE; plot MLT-UT and X-UT.


;---Handle which candidates are plotted.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    candidates = azim_df_search_candidate(project=project)
    if keyword_set(test_time) then begin
        candidate_ids = list()
        foreach candidate, candidates do begin
            if product(candidate.time_range-test_time) gt 0 then continue else candidate_ids.add, candidate.id
        endforeach
        candidate_ids = candidate_ids.toarray()-1   ; id starts from 1.
        candidates = candidates[candidate_ids]
    endif

    ; Filter by region_names and search_names.
    regions = project.search_candidate.search_roi.regions
    if n_elements(region_names) eq 0 then begin
        region_names = list()
        foreach region, regions do region_names.add, region.name
    endif
    if ~isa(region_names,'list') then region_names = list(region_names, /extract)

    search_settings = project.search_settings
    if n_elements(search_names) eq 0 then begin
        search_names = list()
        foreach search_setting, search_settings do search_names.add, search_setting.name
    endif
    if ~isa(search_names,'list') then search_names = list(search_names, /extract)

    candidate_ids = list()
    foreach candidate, candidates do begin
        if region_names.where(candidate.region) eq !null then continue
        if search_names.where(candidate.search_type) eq !null then continue
        candidate_ids.add, candidate.id
    endforeach
    candidate_ids = candidate_ids.toarray()-1
    candidates = candidates[candidate_ids]

    ; Further filter by scale down the # of plots.
    ncandidate = n_elements(candidates)
    expected_plot_number = 50
    if ncandidate gt expected_plot_number then begin
        step = ceil(ncandidate/expected_plot_number)
        candidate_ids = make_bins([0,ncandidate], step, /inner)
        candidates = candidates[candidate_ids]
    endif


;---Settings for detrend tilt angle.
    time_step = project.time_step
    smooth_window = constant('secofhour')
    smooth_width = smooth_window/time_step

;---Loop through candidate to generate plots.
    probe_infos = project.probe_infos
    foreach search_setting, search_settings do begin
        if search_names.where(search_setting.name) eq !null then continue
        ; Load data for one time for all candidates in this search.
        time_range = search_setting.time_range
        foreach probe, search_setting.probes do begin
            prefix = probe_infos[probe].prefix
            r_sm_var = prefix+'r_sm'
            if check_if_update(r_sm_var, time_range) then azim_df_read_data, 'r_sm', time_range=time_range, probe=probe, project=project
            mlt_var = prefix+'pseudo_mlt'
            if check_if_update(mlt_var, time_range) then begin
                get_data, r_sm_var, times, rsm
                mlt = azim_df_calc_pseudo_mlt(rsm)
                store_data, mlt_var, times, mlt
            endif
            x_sm_var = prefix+'x_sm'
            if check_if_update(x_sm_var, time_range) then begin
                get_data, r_sm_var, times, rsm
                store_data, x_sm_var, times, rsm[*,0]
            endif
            theta_var = prefix+'theta'
            if check_if_update(theta_var, time_range) then azim_df_read_data, 'theta', time_range=time_range, probe=probe, project=project
        endforeach
        dst_var = 'dst'
        if check_if_update(dst_var, time_range) then omni_read_index, time_range

        search_name = search_setting.name
        foreach region, regions do begin
            region_name = region.name
            if region_names.where(region_name) eq !null then continue
            foreach candidate, candidates do begin
                if candidate.search_type ne search_name then continue
                if candidate.region ne region_name then continue
                plot_file = join_path([project.plot_dir,'diagnostic_plot','search_candidate','summary_plot',project.name+'_candidate_'+strjoin(time_string(candidate.time_range,tformat='YYYY_MMDD_hhmm'),'_')+'_v01.pdf'])
                azim_df_gen_diagnostic_plot_search_candidate_plot, candidate, project=project, filename=plot_file
            endforeach
        endforeach
    endforeach

end




pro azim_df_gen_diagnostic_plot_search_candidate_search_triad, project=project
    ; Internal, plot sc location, ROI.

;test = 1
;test_time = time_double('2007-12-08/12:15')

    if n_elements(project) eq 0 then project = azim_df_load_project()
    triad_candidates = azim_df_search_candidate(project=project, /get_triad_candidates)
    ntriad_candidate = n_elements(triad_candidates)
    step = ceil(ntriad_candidate/50)
    index = make_bins([0,ntriad_candidate], step, /inner)
    candidates = triad_candidates[index]


    ncircle = 50
    tmp = smkarthm(0,1,ncircle,'n')*2*!dpi
    circle_x = cos(tmp)
    circle_y = sin(tmp)

    tmp = smkarthm(0,1,ncircle,'n')*!dpi
    magn_test = fltarr(ncircle,3)
    magn_test[*,0] = circle_x*30
    magn_test[*,1] = circle_y*30
    tmp = check_if_in_magn(magn_test, magn_pos=magn_gsm)
    magn_gsm = magn_gsm[*,[0,2,1]]
    magn_gsm = [magn_gsm,magn_gsm]
    magn_gsm[ncircle:ncircle*2-1,1] *= -1


;---Settings for figure.
    panel_xsize = 4.
    panel_aspect_ratio = 0.20
    panel_ysize = panel_xsize
    ypans = [panel_ysize]
    nypanel = n_elements(ypans)
    nxpanel = 1
    ypads = 0
    margins = [10,5,8,2]
    sgopen, 0, xsize=1,ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    fig_xsize = panel_xsize*nxpanel+total(margins[[0,2]])*abs_xchsz
    fig_ysize = total(ypans)+total(margins[[1,3]])*abs_ychsz+total(ypads)*abs_ychsz
    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    sgopen, 0, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz
    poss = sgcalcpos(nypanel, position=pos, ypans=ypans, ypad=ypads)
    sgclose, /wdelete
    var_labels = string(findgen(nypanel)+1,format='(I0)')+'.'

    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    full_ychsz = constant('full_ychsz')
    half_ychsz = full_ychsz*0.5
    lineskip = constant('lineskip')
    label_size = 0.7
    label_xshift = 1
    label_yshift = full_ychsz
    secofday = constant('secofday')
    bar_thick = keyword_set(test)? 0.5: 4


    mlt_limits = list()
    foreach region, project.search_candidate.search_roi.regions do mlt_limits.add, region.mlt_range
    mlt_limits = sort_uniq(mlt_limits.toarray())
    probe_infos = project.probe_infos

    ; Draw triad at the middle time.
    triad_color = sgcolor('red')
    nvertex = 3
    ndim = 3
    triad_angle_range = project.search_candidate.search_triad.triad_angle_range


;---Loop through candidate for each search.
    search_settings = project.search_settings
    foreach search_setting, search_settings do begin
        full_time_range = search_setting.time_range
        search_name = search_setting.name
        triad_min_count = search_setting.triad_min_count
;if search_name eq 'beyond_15Re' then continue

        ; Load data for 1 time.
        probes = search_setting.probes
        foreach probe, probes do azim_df_read_data, 'r_sm', time_range=full_time_range, probe=probe, project=project

        foreach candidate, candidates do begin
            search_type = candidate.search_type
            if search_type ne search_name then continue

            if keyword_set(test_time) then begin
                lprmsg, time_string(candidate.time_range)
                if product(candidate.time_range-test_time) gt 0 then continue else stop
            endif

            for section_id=0, candidate.nsection-1 do begin
                section_time_range = candidate.time_range_list[section_id]
                section_probes = candidate.probe_list[section_id]
                duration = total(section_time_range*[-1,1])
                if duration lt 3600 then continue
                region = candidate.region

                base_name = 'triad_'+strjoin(time_string(section_time_range,tformat='YYYY_MMDD_hhmm'),'_')+'.pdf'
                plot_file = join_path([project.plot_dir,'diagnostic_plot','search_candidate','search_triad',base_name])
                if keyword_set(test) then plot_file = test
                if keyword_set(test_time) then plot_file = test
                sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz

            ;---Panel a. SC location.
                tpos = poss[*,0]
                yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
                xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])

                xrange = [20.,-40]
                xstep = 20
                xtickv = make_bins(xrange, xstep)
                xticks = n_elements(xtickv)-1
                xminor = 5
                xtitle = 'SM X (Re)'
                xticklen = xticklen

                yrange = (region eq 'pre_midn')? [40.,-20]: [20.,-40]
                ystep = 20
                ytickv = make_bins(yrange, ystep)
                yticks = n_elements(ytickv)-1
                yminor = 5
                ytitle = 'SM Y (Re)'
                yticklen = yticklen

                ; Set up coord.
                plot, xrange, yrange, $
                    xstyle=5, xrange=xrange, xtickformat='(A1)', $
                    ystyle=5, yrange=yrange, ytickformat='(A1)', $
                    position=tpos, /noerase, /nodata, /iso

                ; Add earth.
                polyfill, circle_x<0, circle_y, color=sgcolor('silver')
                plots, circle_x, circle_y

                ; Add magnteopause.
                oplot, magn_gsm[*,0], magn_gsm[*,1]

                ; Add more circles.
                foreach tmp, [5,10,20,40] do oplot, circle_x*tmp, circle_y*tmp, linestyle=1
                foreach tmp, [4,40] do oplot, circle_x*tmp, circle_y*tmp

                ; Add mlt_range.
                rad = constant('rad')
                foreach tmp, mlt_limits do oplot, -[4,40]*cos(tmp*15*rad), -[4,40]*sin(tmp*15*rad)

                foreach probe, probes do begin
                    prefix = probe_infos[probe].prefix
                    color = probe_infos[probe].color
                    rsm = get_var_data(prefix+'r_sm', in=section_time_range)
                    xxs = rsm[*,0]
                    yys = rsm[*,1]
                    index = where(section_probes eq probe, count)
                    thick = (count eq 0)? 0.5*bar_thick: 2*bar_thick
                    oplot, xxs, yys, color=color, thick=thick
                endforeach


                ; Add triads.
                ;triad_time = section_time_range[0]
                ;triad_time = section_time_range[1]
                triad_time = mean(section_time_range)
                ;triad_time = triad_time-(triad_time mod project.search_candidate_time_step)
                nprobe = n_elements(section_probes)
                rsms = fltarr(ndim,nprobe)
                foreach probe, section_probes, ii do rsms[*,ii] = get_var_data(probe+'_r_sm', at=triad_time)
                rsms[2,*] = 0

                if search_name eq 'beyond_15Re' then begin
                    combos = list()
                    foreach probe, ['a','d','e'] do combos.add, 'th'+['b','c',probe]
                endif else combos = choose_from(section_probes, nvertex)
                foreach combo, combos do combo = combo[sort(combo)]
                ncombo = n_elements(combos)
                combo_flags = bytarr(ncombo)

                foreach combo, combos do begin
                    combo_index = (combos.where(combo))[0]
                    probe_index = intarr(nvertex)
                    foreach probe, combo, ii do probe_index[ii] = where(section_probes eq probe)
                    triad_angles = fltarr(nvertex)
                    foreach probe, combo, ii do begin
                        probe_index = shift(probe_index, 1)
                        triad_angles[ii] = sang($
                            rsms[*,probe_index[1]]-rsms[*,probe_index[0]], $
                            rsms[*,probe_index[2]]-rsms[*,probe_index[0]], /deg)
                    endforeach
                    index = lazy_where(triad_angles, '[]', triad_angle_range, count=count)
                    combo_flags[combo_index] = count eq nvertex
                endforeach
                index = where(combo_flags eq 1, triad_count)
                if triad_count lt triad_min_count then message, 'Inconsistency, stop here ...'
                combos = combos[index]
                foreach combo, combos do begin
                    probe_index = intarr(nvertex)
                    foreach probe, combo, ii do probe_index[ii] = where(section_probes eq probe)
                    foreach probe, combo, ii do begin
                        probe_index = shift(probe_index, 1)
                        oplot, rsms[0,probe_index[0:1]], rsms[1,probe_index[0:1]], color=triad_color
                    endforeach
                endforeach


                ; Add notations.
                ty0 = (region eq 'post_midn')? tpos[3]: tpos[1]+ychsz*full_ychsz*(3+lineskip)
                tx = tpos[0]+xchsz*1
                step = 4

                ty = ty0-ychsz*full_ychsz*1
                xyouts, tx,ty,/normal, 'Probes of the event: ', charsize=label_size
                ttx = tx+15*xchsz*label_size
                foreach probe, section_probes do begin
                    probe_info = probe_infos[probe]
                    xyouts, ttx,ty,/normal, strupcase(probe_info.short_name), color=probe_info.color, charsize=label_size
                    ttx += step*xchsz*label_size
                endforeach
                xyouts, tpos[2]-xchsz*1,ty,/normal, alignment=1, '# of good triad: '+string(triad_count,format='(I0)'), charsize=label_size


                ty = ty0-ychsz*full_ychsz*2
                xyouts, tx,ty,/normal, 'All available probes: ', charsize=label_size
                ttx = tx+15*xchsz*label_size
                foreach probe, probes do begin
                    probe_info = probe_infos[probe]
                    xyouts, ttx,ty,/normal, strupcase(probe_info.short_name), color=probe_info.color, charsize=label_size
                    ttx += step*xchsz*label_size
                endforeach

                ty = ty0-ychsz*full_ychsz*3
                duration = total(section_time_range*[-1,1])/constant('secofhour')
                msg = 'Orbit from '+strjoin(time_string(section_time_range,tformat='YYYY-MM-DD/hh:mm'), ' to ')+', duration is '+sgnum2str(duration,ndec=1)+' hour(s)'
                xyouts, tx,ty,/normal, msg, charsize=label_size


                ; Add axes.
                plot, xrange, yrange, $
                    xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtitle=xtitle, xticklen=xticklen, xtickformat=xtickformat, $
                    ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, yticklen=yticklen, $
                    position=tpos, /noerase, /nodata

                ; Add labels.
                fig_label = 'a. XY!C    plane'
                tx = label_xshift*xchsz
                ty = tpos[3]-ychsz*full_ychsz
                xyouts, tx,ty,/normal, fig_label


;                if keyword_set(test) then stop
                sgclose

            endfor
        endforeach
    endforeach

end



pro azim_df_gen_diagnostic_plot_search_candidate_search_roi, project=project
    ; Internal, plot sc location, ROI.

test = 0
    if n_elements(project) eq 0 then project = azim_df_load_project()
    roi_candidates = azim_df_search_candidate(project=project, /get_roi_candidates)
    nroi_candidate = n_elements(roi_candidates)
    step = ceil(nroi_candidate/50)
    index = make_bins([0,nroi_candidate], step, /inner)
    candidates = roi_candidates[index]


    ncircle = 50
    tmp = smkarthm(0,1,ncircle,'n')*2*!dpi
    circle_x = cos(tmp)
    circle_y = sin(tmp)

    tmp = smkarthm(0,1,ncircle,'n')*!dpi
    magn_test = fltarr(ncircle,3)
    magn_test[*,0] = circle_x*30
    magn_test[*,1] = circle_y*30
    tmp = check_if_in_magn(magn_test, magn_pos=magn_gsm)
    magn_gsm = magn_gsm[*,[0,2,1]]
    magn_gsm = [magn_gsm,magn_gsm]
    magn_gsm[ncircle:ncircle*2-1,1] *= -1


;---Settings for figure.
    panel_xsize = 4.
    panel_aspect_ratio = 0.20
    panel_ysize = panel_xsize
    ypans = [panel_ysize]
    nypanel = n_elements(ypans)
    nxpanel = 1
    ypads = 0
    margins = [10,5,8,2]
    sgopen, 0, xsize=1,ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    fig_xsize = panel_xsize*nxpanel+total(margins[[0,2]])*abs_xchsz
    fig_ysize = total(ypans)+total(margins[[1,3]])*abs_ychsz+total(ypads)*abs_ychsz
    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    sgopen, 0, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz
    poss = sgcalcpos(nypanel, position=pos, ypans=ypans, ypad=ypads)
    sgclose, /wdelete
    var_labels = string(findgen(nypanel)+1,format='(I0)')+'.'

    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    full_ychsz = constant('full_ychsz')
    half_ychsz = full_ychsz*0.5
    lineskip = constant('lineskip')
    label_size = 0.7
    label_xshift = 1
    label_yshift = full_ychsz
    secofday = constant('secofday')
    bar_thick = keyword_set(test)? 0.5: 4


    mlt_limits = list()
    foreach region, project.search_candidate.search_roi.regions do mlt_limits.add, region.mlt_range
    mlt_limits = sort_uniq(mlt_limits.toarray())
    probe_infos = project.probe_infos


;---Loop through candidate for each search.
    search_settings = project.search_settings
    foreach search_setting, search_settings do begin
        full_time_range = search_setting.time_range
        search_name = search_setting.name

        ; Load data for 1 time.
        probes = search_setting.probes
        foreach probe, probes do azim_df_read_data, 'r_sm', time_range=full_time_range, probe=probe, project=project

        foreach candidate, candidates do begin
            search_type = candidate.search_type
            if search_type ne search_name then continue

            for section_id=0, candidate.nsection-1 do begin
                section_time_range = candidate.time_range_list[section_id]
                section_probes = candidate.probe_list[section_id]
                duration = total(section_time_range*[-1,1])
                if duration lt 3600 then continue
                region = candidate.region

                if keyword_set(test_time) then begin
                    index = lazy_where(test_time, '[]', section_time_range, count=count)
                    if count ne 0 then stop
                endif

                base_name = 'roi_'+strjoin(time_string(section_time_range,tformat='YYYY_MMDD_hhmm'),'_')+'.pdf'
                plot_file = join_path([project.plot_dir,'diagnostic_plot','search_candidate','search_roi',base_name])
                if keyword_set(test) then plot_file = test
                if keyword_set(test_time) then plot_file = test
                sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz


            ;---Panel a. SC location.
                tpos = poss[*,0]
                yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
                xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])

                xrange = [20.,-40]
                xstep = 20
                xtickv = make_bins(xrange, xstep)
                xticks = n_elements(xtickv)-1
                xminor = 5
                xtitle = 'SM X (Re)'
                xticklen = xticklen

                yrange = (region eq 'pre_midn')? [40.,-20]: [20.,-40]
                ystep = 20
                ytickv = make_bins(yrange, ystep)
                yticks = n_elements(ytickv)-1
                yminor = 5
                ytitle = 'SM Y (Re)'
                yticklen = yticklen

                ; Set up coord.
                plot, xrange, yrange, $
                    xstyle=5, xrange=xrange, xtickformat='(A1)', $
                    ystyle=5, yrange=yrange, ytickformat='(A1)', $
                    position=tpos, /noerase, /nodata, /iso

                ; Add earth.
                polyfill, circle_x<0, circle_y, color=sgcolor('silver')
                plots, circle_x, circle_y

                ; Add magnteopause.
                oplot, magn_gsm[*,0], magn_gsm[*,1]

                ; Add more circles.
                foreach tmp, [5,10,20,40] do oplot, circle_x*tmp, circle_y*tmp, linestyle=1
                foreach tmp, [4,40] do oplot, circle_x*tmp, circle_y*tmp

                ; Add mlt_range.
                rad = constant('rad')
                foreach tmp, mlt_limits do oplot, -[4,40]*cos(tmp*15*rad), -[4,40]*sin(tmp*15*rad)

                foreach probe, probes do begin
                    prefix = probe_infos[probe].prefix
                    color = probe_infos[probe].color
                    rsm = get_var_data(prefix+'r_sm', in=section_time_range)
                    xxs = rsm[*,0]
                    yys = rsm[*,1]
                    index = where(section_probes eq probe, count)
                    thick = (count eq 0)? 0.5*bar_thick: 2*bar_thick
                    oplot, xxs, yys, color=color, thick=thick
                endforeach


                ; Add notations.
                ty0 = (region eq 'post_midn')? tpos[3]: tpos[1]+ychsz*full_ychsz*(3+lineskip)
                tx = tpos[0]+xchsz*1
                step = 4

                ty = ty0-ychsz*full_ychsz*1
                xyouts, tx,ty,/normal, 'Probes of the event: ', charsize=label_size
                ttx = tx+15*xchsz*label_size
                foreach probe, section_probes do begin
                    probe_info = probe_infos[probe]
                    xyouts, ttx,ty,/normal, strupcase(probe_info.short_name), color=probe_info.color, charsize=label_size
                    ttx += step*xchsz*label_size
                endforeach

                ty = ty0-ychsz*full_ychsz*2
                xyouts, tx,ty,/normal, 'All available probes: ', charsize=label_size
                ttx = tx+15*xchsz*label_size
                foreach probe, probes do begin
                    probe_info = probe_infos[probe]
                    xyouts, ttx,ty,/normal, strupcase(probe_info.short_name), color=probe_info.color, charsize=label_size
                    ttx += step*xchsz*label_size
                endforeach

                ty = ty0-ychsz*full_ychsz*3
                duration = total(section_time_range*[-1,1])/constant('secofhour')
                msg = 'Orbit from '+strjoin(time_string(section_time_range,tformat='YYYY-MM-DD/hh:mm'), ' to ')+', duration is '+sgnum2str(duration,ndec=1)+' hour(s)'
                xyouts, tx,ty,/normal, msg, charsize=label_size

                ; Add axes.
                plot, xrange, yrange, $
                    xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtitle=xtitle, xticklen=xticklen, xtickformat=xtickformat, $
                    ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, yticklen=yticklen, $
                    position=tpos, /noerase, /nodata

                ; Add labels.
                fig_label = 'a. XY!C    plane'
                tx = label_xshift*xchsz
                ty = tpos[3]-ychsz*full_ychsz
                xyouts, tx,ty,/normal, fig_label


                if keyword_set(test) then stop
                sgclose
            endfor
        endforeach
    endforeach

end





pro azim_df_gen_diagnostic_plot_primitive_data, project=project
    ; Internal, plot dis and bmag on month-base.

    if n_elements(project) eq 0 then project = azim_df_load_project()
    search_settings = project.search_settings
    foreach search_setting, search_settings do begin
        time_range = search_setting.time_range
        probes = search_setting.probes
test = 0
        ; See if data file exists.
        data_file = join_path([project.data_dir,search_setting.data_file_suffix])
        if file_test(data_file) eq 0 then begin
            foreach probe, probes do azim_df_load_primitive_data, time_range=time_range, probe=probe, project=project, data_file=data_file
        endif

        ; Setting.
        dis_range = [4.,40]
        bmag_range = [1,1000]
        plot_vars = []
        foreach probe, probes do plot_vars = [plot_vars,probe+'_'+['dis','bmag']]


        margins = [10,5,8,2]
        nypanel = n_elements(plot_vars)
        nxpanel = 1
        ypans = intarr(nypanel)+1
        ypads = intarr(nypanel)+0.5
        ypads[1:*:2] = 1.5    ; increase separation between probes.
        panel_ysize = 0.6
        panel_xsize = 8.
        sgopen, 0, xsize=1,ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
        sgclose, /wdelete
        fig_xsize = panel_xsize*nxpanel+total(margins[[0,2]])*abs_xchsz
        fig_ysize = panel_ysize*nypanel+total(margins[[1,3]])*abs_ychsz+total(ypads)*abs_ychsz
        pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
        pos = [pos[0:1],1-pos[2:3]]
        sgopen, 0, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz
        poss = sgcalcpos(nypanel, position=pos, ypans=ypans, ypad=ypads)
        sgclose, /wdelete
        var_labels = string(findgen(nypanel)+1,format='(I0)')+'.'


        ; Load data.
        probe_infos = dictionary()
        foreach probe, probes do probe_infos[probe] = resolve_probe(probe)
        foreach probe, probes do begin
            probe_label = strupcase(probe_infos[probe].short_name)
            prefix = probe_infos[probe].prefix

            the_var = prefix+'dis'
            if check_if_update(the_var, time_range) then begin
                azim_df_read_data, 'r_gsm', time_range=time_range, probe=probe, project=project
                get_data, prefix+'r_gsm', times, r_gsm
                store_data, the_var, times, snorm(r_gsm)
                add_setting, the_var, /smart, {$
                    display_type: 'scalar', $
                    short_name: 'R', $
                    unit: 'Re', $
                    yrange: dis_range, $
                    ylog: 1}
                options, the_var, 'labels', probe_label+' |R|'
            endif

            the_var = prefix+'bmag'
            if check_if_update(the_var, time_range) then begin
                azim_df_read_data, 'b_gsm', time_range=time_range, probe=probe, project=project
                get_data, prefix+'b_gsm', times, b_gsm
                store_data, the_var, times, snorm(b_gsm)
                add_setting, the_var, /smart, {$
                    display_type: 'scalar', $
                    short_name: '|B|', $
                    unit: 'nT', $
                    yrange: bmag_range, $
                    ylog: 1}
                options, the_var, 'labels', probe_label+' |B|'
            endif
        endforeach

        ; Plot on month-base.
        time_ranges = sort_uniq([time_range,break_down_times(time_range,'month')])
        ntime_range = n_elements(time_ranges)-1
        for ii=0, ntime_range-1 do begin
            the_time_range = time_ranges[ii:ii+1]
            base_name = search_setting.name+'_'+strjoin(time_string(the_time_range,tformat='YYYY_MM'),'_')+'.pdf'
            plot_file = join_path([project.plot_dir,'diagnostic_plot','primitive_data',base_name])
            if keyword_set(test) then plot_file = test
            sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz
            tplot, plot_vars, position=poss, trange=the_time_range
            foreach var, plot_vars, jj do begin
                tpos = poss[*,jj]
                tx = 2*xchsz
                ty = tpos[3]-ychsz*constant('full_ychsz')
                xyouts, tx,ty,/normal, var_labels[jj]
            endforeach
            if keyword_set(test) then stop
            sgclose
        endfor
    endforeach

end

pro azim_df_gen_diagnostic_plot, type, project=project, _extra=ex

    if n_elements(type) eq 0 then message, 'No type, stop here ...'
    call_procedure, 'azim_df_gen_diagnostic_plot_'+type, project=project, _extra=ex

end

; azim_df_gen_diagnostic_plot, 'primitive_data', project=project
; azim_df_gen_diagnostic_plot, 'search_candidate_search_roi', project=project
;region_names = ['pre','post']+'_midn'
;foreach region_name, region_names do azim_df_gen_diagnostic_plot, 'search_candidate', project=project, search_name='within_15Re', region_name=region_name;, test_time=time_double('2014-08-28/10:00')

;; This is for searching large DF events.
;region_names = ['pre','post']+'_midn'
;search_names = ['within_15Re','beyond_15Re']
;foreach search_name, search_names do foreach region_name, region_names do $
;    azim_df_gen_diagnostic_plot, 'search_event_search_large_df', project=project, $
;    region_name=region_name, search_name=search_name

; This is for searching DF groups.
azim_df_gen_diagnostic_plot, 'search_event_search_df_group', project=project

end
