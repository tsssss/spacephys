;+
; DP after time lag removed, and with keV electron.
;-

    test = 0
    short_time_range = time_double(['2016-10-13/12:00','2016-10-13/13:30'])

;---Check inputs.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(short_time_range) ne 2 then begin
        errmsg = handle_error('No input short_time_range ...')
        return
    endif

;---Find the event.
    events = azim_df_search_all_events(project=project)
    foreach event, events do if product(event.time_range-short_time_range) lt 0 then break
    azim_df_load_basic_data, project=project

;---Sort based on obs_time.
    df_list = event.df_list
    ndf = df_list.length
    obs_times = dblarr(ndf)
    foreach df, df_list, ii do obs_times[ii] = df.obs_time
    index = sort(obs_times)
    obs_times = obs_times[index]
    df_list = df_list[index]
    probes = strarr(ndf)
    foreach df, df_list, ii do probes[ii] = df.probe
    sorted_probes = probes
    nprobe = n_elements(probes)
    time_lags = obs_times-obs_times[0]
    time_lags[*] = 0

;---Read keV e- flux.
    time_range = short_time_range+[-1,1]*30*60
    if check_if_update(probes[0]+'_kev_e_flux') then begin
        foreach probe, ['a','b'] do rbsp_read_kev_electron, time_range, probe=probe
        foreach probe, ['d'] do themis_read_kev_electron, time_range, probe=probe
        foreach probe, ['13','14','15'] do goes_read_kev_electron, time_range, probe=probe
    endif

    ; select a few energy bands.
    nenergy = 4
    foreach probe, probes, probe_id do begin
        var = probe+'_kev_e_flux'
        get_data, var, times, eflux, ens, limits=lim
        nen = n_elements(eflux[0,*])
        the_index = smkarthm(0,nen-1,nenergy+1,'n')
        index = round(the_index[0:-2])
        index = findgen(nenergy)
        store_data, var+'_simple', times, $
            eflux[*,index], ens[index], limits=lim
        options, var+'_simple', 'colors', lim.colors[index]
        options, var+'_simple', 'labels', lim.labels[index]
        options, var+'_simple', 'ytitle', ' '

        the_time_lag = time_lags[probe_id]
        data = get_var_data(var+'_simple', in=short_time_range+the_time_lag)
        if n_elements(data) ne 0 then begin
            log_yrange = alog10(minmax(data))>1
            log_yrange = [floor(log_yrange[0]),ceil(log_yrange[1])]
            yrange = 10d^log_yrange
            ylim, var+'_simple', yrange
            log_ytickv = make_bins(log_yrange, 1)
            ytickv = 10d^log_ytickv
            ytickn = '10!U'+string(log_ytickv,format='(I0)')
            ytickn[0:*:2] = ' '
            yminor = 10
            yticks = n_elements(ytickv)-1
            options, var+'_simple', 'ytickv', ytickv
            options, var+'_simple', 'yticks', yticks
            options, var+'_simple', 'yminor', yminor
            options, var+'_simple', 'ytickname', ytickn
        endif else begin
            xxs = short_time_range+the_time_lag
            yys = fltarr(2,nenergy)+1
            store_data, var+'_simple', xxs, yys
            options, var+'_simple', 'ytickv', [1,10]
            options, var+'_simple', 'yticks', 2
            options, var+'_simple', 'yminor', 10
            options, var+'_simple', 'ytickname', ['1','10']
        endelse

;        tplot, probe+'_'+['kev_e_flux'+['','_simple'],'theta'], trange=short_time_range
;        stop
    endforeach

;---Prepare colors.
    probe_colors = smkarthm(90,240,nprobe, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)


;---Get the size of the figure, and position of the panels.
    plot_file = join_path([homedir(),'Dropbox','mypapers','dp_vs_fac','plot',$
        'fig_dp_vs_injection_2016_1013_v01.pdf'])
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=6, ysize=5, /inch, xchsz=xchsz, ychsz=ychsz

    xpad = [5]
    ypad = 0.4
    margins = [10,3,5,2]
    all_poss = sgcalcpos(nprobe,2, xpad=xpad, ypad=ypad, margins=margins)
    label_size = constant('label_size')
    half_ychsz = constant('half_ychsz')

;---Common x-axis.
    xrange = short_time_range
    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    xminor = 10 ; min.
    xstep = 30*60
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    secofday = 86400d
    for ii=0, xticks do begin
        the_time = xtickv[ii]
        xtickn[ii] = time_string(the_time,tformat='hh:mm')
        date = time_string(the_time,tformat='YYYY-MM-DD')+'        '
        if ii eq 0 then begin
            xtickn[ii] = date
            continue
        endif
    endfor


    line_colors = sgcolor(['black','gray'])
    foreach probe, sorted_probes, probe_id do begin
        prefix = probe+'_'
        theta = get_var_data(prefix+'theta', times=times)
        index = where_pro(times, '[]', short_time_range)
        yrange = minmax(theta[index])
        yrange = [-1,1]*max(abs(make_bins(yrange,2)))
        yticks = 2
        if max(yrange) le 5 then begin
            yminor = max(yrange)
            ytickv = [-1,0,1]*yminor
        endif else begin
            yminor = 5
            ytickv = [-1,0,1]*floor(max(yrange)/5)*5
        endelse
        xtickformat = (probe_id eq nprobe-1)? '': '(A1)'


        tpos = reform(all_poss[*,0,probe_id])

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        ytickformat = ''
        ytitle = ' '
        the_time_range = short_time_range+time_lags[probe_id]

        index = where_pro(times, '[]', the_time_range)
        xxs = times[index]
        xxs -= time_lags[probe_id]
        yys = theta[index]

        the_xtickn = xtickn

        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, /nodata, /noerase

        oplot, xxs, yys, color=line_colors[0]
        plots, xrange, [0,0], linestyle=1, color=line_colors[0]
        ref_time = obs_times[probe_id]
        ref_time -= time_lags[probe_id]
        plots, ref_time, yrange, linestyle=1, color=line_colors[0]


        ; Draw box.
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=the_xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytickformat=ytickformat, ytitle=ytitle, $
            position=tpos, /nodata, /noerase


        ; Add title.
        if probe_id eq 0 then begin
            tx = (tpos[0]+tpos[2])*0.5
            ty = tpos[3]+ychsz*0.5
            title = tex2str('theta')+' (deg), time lag removed'
            xyouts, tx,ty,/normal, title, color=line_colors[0], alignment=0.5
        endif

        ; Add label.
        label = 'a-'
        label += string(probe_id+1,format='(I0)')+'.'
        ty = tpos[3]-ychsz*1
        tx = tpos[0]+xchsz*0.5;: tpos[2]-xchsz*0.5
        alignment = 0;: 1
        xyouts, tx,ty,/normal, label, alignment=alignment


        ; Add probe info.
        tx = tpos[0]-xchsz*5
        ty = (tpos[1]+tpos[3])*0.5+ychsz*half_ychsz
        probe_info = resolve_probe(probe)
        xyouts, tx,ty,/normal, alignment=1, strupcase(probe_info.short_name), color=probe_colors[probe_id]
        df = df_list[probe_id]
        ref_mlt = df.obs_mlt
        ref_rxy = df.obs_rxy
        ty = (tpos[1]+tpos[3])*0.5-ychsz*half_ychsz
        xyouts, tx,ty,/normal, alignment=1, sgnum2str(ref_mlt,ndec=1)+' MLT!C'+$
            sgnum2str(ref_rxy,ndec=1)+' Re', charsize=label_size
    endforeach


;---Draw eflux.
    vars = probes+'_kev_e_flux_simple'
    poss = reform(all_poss[*,1,*])
    foreach probe, sorted_probes, probe_id do begin
        prefix = probe+'_'
        the_time_lag = time_lags[probe_id]


        var = prefix+'kev_e_flux_simple'
        yys = get_var_data(var, in=short_time_range+the_time_lag, times=xxs, limits=lim)
        yrange = lim.yrange
        yticks = lim.yticks
        yminor = lim.yminor
        ytickv = lim.ytickv
        ytickname = lim.ytickname
        xtickformat = (probe_id eq nprobe-1)? '': '(A1)'


        tpos = reform(all_poss[*,1,probe_id])

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        ytickformat = '';: '(A1)'
        ytitle = lim.ytitle

        plot, xrange, yrange, $
            xstyle=5, xlog=0, xrange=xrange, $
            ystyle=5, ylog=1, yrange=yrange, $
            position=tpos, /nodata, /noerase

        foreach color, lim.colors, jj do oplot, xxs-the_time_lag, yys[*,jj], color=color
        ref_time = obs_times[probe_id]-the_time_lag
        plots, ref_time, yrange, linestyle=1, color=line_colors[0]
;        if jj eq 0 then begin
;            tmp = convert_coord(ref_time,yrange[0], /data, /to_normal)
;            tx = tmp[0]+xchsz*0.5
;            ty = tmp[1]+ychsz*0.25
;            xyouts, tx,ty,/normal, time_string(ref_time,tformat='hh:mm:ss'), charsize=label_size
;        endif

        labels = lim.labels
        ndy = n_elements(labels)+1
        ys = smkarthm(tpos[3],tpos[1],ndy,'n')
        for jj=0,ndy-2 do begin
            tx = tpos[2]+xchsz*0.8
            ty = ys[jj]-ychsz*label_size
            xyouts, tx,ty,normal=1, labels[jj], color=lim.colors[jj], charsize=label_size
        endfor

        ; Draw box.
        the_xtickn = xtickn
        the_xtickn[0] = ' '

        plot, xrange, yrange, $
            xstyle=1, xlog=0, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=the_xtickn, xtickformat=xtickformat, $
            ystyle=1, ylog=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, ytickname=ytickname, $
            position=tpos, /nodata, /noerase


        ; Add title.
        if probe_id eq 0 then begin
            tx = (tpos[0]+tpos[2])*0.5
            ty = tpos[3]+ychsz*0.5
            title = 'e- flux (#/cm!U2!N-s-sr-keV), time lag removed'
            xyouts, tx,ty,/normal, title, color=line_colors[0], alignment=0.5
        endif

        ; Add label.
        label = 'b-'
        label += string(probe_id+1,format='(I0)')+'.'
        ty = tpos[3]-ychsz*1
        tx = tpos[0]+xchsz*0.5
        alignment = 0
        xyouts, tx,ty,/normal, label, alignment=alignment
    endforeach
;    tplot, vars, trange=short_time_range, position=poss, noerase=1


    if keyword_set(test) then stop
    sgclose
end
