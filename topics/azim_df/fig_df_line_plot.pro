;+
; DF in real time and after time lag removed.
;-

    test = 0
    short_time_range = time_double(['2014-08-28/09:30','2014-08-28/11:30'])

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

;---Prepare colors.
    probe_colors = sgcolor(['firebrick','orange','sea_green','deep_sky_blue','medium_blue','dark_violet'])
    probe_colors = smkarthm(90,240,nprobe, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)


;---Get the size of the figure, and position of the panels.
    plot_file = join_path([project.plot_dir,'fig_df_line_plot.pdf'])
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=6, ysize=5, /inch, xchsz=xchsz, ychsz=ychsz

    xpad = 2
    ypad = 0.4
    margins = [12.,3,2,2]
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

    time_lags = obs_times-obs_times[0]
    line_colors = sgcolor(['black','gray'])
    foreach probe, sorted_probes, ii do begin
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
        xtickformat = (ii eq nprobe-1)? '': '(A1)'


        for jj=0,1 do begin
            tpos = reform(all_poss[*,jj,ii])

            xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
            yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
            ytickformat = (jj eq 0)? '': '(A1)'
            ytitle = (jj eq 0)? '(deg)': ' '
            the_time_range = (jj eq 0)? short_time_range: short_time_range+time_lags[ii]

            index = where_pro(times, '[]', the_time_range)
            xxs = times[index]
            if jj eq 1 then xxs -= time_lags[ii]
            yys = theta[index]

            the_xtickn = xtickn
            if jj eq 1 then the_xtickn[0] = ' '

            plot, xrange, yrange, $
                xstyle=5, xrange=xrange, $
                ystyle=5, yrange=yrange, $
                position=tpos, /nodata, /noerase

            oplot, xxs, yys, color=line_colors[jj]
            plots, xrange, [0,0], linestyle=1, color=line_colors[jj]
            ref_time = obs_times[ii]
            if jj eq 1 then ref_time -= time_lags[ii]
            plots, ref_time, yrange, linestyle=1, color=line_colors[jj]
            if jj eq 0 then begin
                tmp = convert_coord(ref_time,yrange[0], /data, /to_normal)
                tx = tmp[0]+xchsz*0.5
                ty = tmp[1]+ychsz*0.25
                xyouts, tx,ty,/normal, time_string(ref_time,tformat='hh:mm:ss'), charsize=label_size
            endif

;            ; Add width and height.
;            if ii eq 0 and jj eq 1 then begin
;                df = df_list[ii]
;
;                width = df.width
;                txs = ref_time+[-1,1]*0.5*width
;                foreach tx, txs, kk do begin
;                    tmp = convert_coord(tx,yrange[0], /data, /to_normal)
;                    txs[kk] = tmp[0]
;                endforeach
;                tmp = convert_coord(xrange[0],mean(df.value_range), /data, /to_normal)
;                ty = tmp[1]
;                tys = ty;+(tpos[3]-tpos[1])*0.1
;                plots, txs,tys, /normal
;                for kk=0,1 do begin
;                    plots, txs[kk]+[0,0], tys[0]+[-1,1]*0.5*ychsz*0.2, /normal
;                endfor
;
;                tys = df.value_range
;                foreach ty, tys, kk do begin
;                    tmp = convert_coord(xrange[0],ty, /data, /to_normal)
;                    tys[kk] = tmp[1]
;                endforeach
;                tmp = convert_coord(ref_time,yrange[0], /data, /to_normal)
;                tx = tmp[0]
;                txs = tx;-(tpos[2]-tpos[0])*0.1
;                plots, txs,tys, /normal
;                for kk=0,1 do begin
;                    plots, txs[0]+[-1,1]*0.5*xchsz*0.5, tys[kk]+[0,0], /normal
;                endfor
;            endif

            plot, xrange, yrange, $
                xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=the_xtickn, xtickformat=xtickformat, $
                ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytickformat=ytickformat, ytitle=ytitle, $
                position=tpos, /nodata, /noerase

            ; Add probe info.
            if jj eq 0 then begin
                tx = tpos[0]-xchsz*7
                ty = (tpos[1]+tpos[3])*0.5+ychsz*half_ychsz
                probe_info = resolve_probe(probe)
                xyouts, tx,ty,/normal, alignment=1, strupcase(probe_info.short_name), color=probe_colors[ii]
                df = df_list[ii]
                ref_mlt = df.obs_mlt
                ref_rxy = df.obs_rxy
                ty = (tpos[1]+tpos[3])*0.5-ychsz*half_ychsz
                xyouts, tx,ty,/normal, alignment=1, sgnum2str(ref_mlt,ndec=1)+' MLT!C'+$
                    sgnum2str(ref_rxy,ndec=1)+' Re', charsize=label_size
            endif

            ; Add title.
            if ii eq 0 then begin
                tx = (tpos[0]+tpos[2])*0.5
                ty = tpos[3]+ychsz*0.3
                title = (jj eq 0)? 'Detrended tilt angle '+tex2str('theta')+' in absolute time': $
                    'Time lag removed'
                xyouts, tx,ty,/normal, title, color=line_colors[jj], alignment=0.5
            endif

            ; Add label.
            label = (jj eq 0)? 'a-': 'b-'
            label += string(ii+1,format='(I0)')+'.'
            ty = tpos[3]-ychsz*1
            tx = (jj eq 0)? tpos[0]+xchsz*0.5: tpos[2]-xchsz*0.5
            alignment = (jj eq 0)? 0: 1
            xyouts, tx,ty,/normal, label, alignment=alignment
        endfor
    endforeach



    if keyword_set(test) then stop
    sgclose
end
