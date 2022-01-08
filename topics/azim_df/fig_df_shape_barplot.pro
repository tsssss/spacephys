;+
; Plot aspect ratio.
;-


test = 0
    magnify = (keyword_set(test))? 2: 1

    if n_elements(project) eq 0 then project = azim_df_load_project()
    events = azim_df_find_dfgroup(project=project)


;test_time = time_double('2014-08-28/10:30')
;if keyword_set(test_time) then begin
;    foreach event, events, event_id do begin
;        if product(event.time_range-test_time) lt 0 then break
;    endforeach
;    events = list(events[event_id])
;endif

    event_info = list()
    ndim = 3
    foreach event, events do begin
        current_info = dictionary()

        ; Only check the azim DFs.
        direction = (strsplit(event.region,'%',/extract))[1]
        if direction eq 'earthward' then continue
        if direction eq 'outward' then continue

        ; event_duration.
        df_list = event.df_list
        ndf = df_list.length
        df_times = dblarr(ndf)
        df_mlts = fltarr(ndf)
        df_rxys = fltarr(ndf)
        foreach df, df_list, ii do begin
            df_times[ii] = df.obs_time
            df_mlts[ii] = df.obs_mlt
            df_rxys[ii] = df.obs_rxy
        endforeach
        time_range = minmax(df_times)
        current_info['time_range'] = time_range
        mlt_range = minmax(df_mlts)
        current_info['mlt_range'] = mlt_range
        rxy_range = minmax(df_rxys)
        current_info['rxy_range'] = rxy_range
;        zsm_range = minmax(df_rsms[*,2])
;        current_info['zsm_range'] = zsm_range
        ;print, time_string(time_range)

        keys = current_info.keys()
        foreach key, keys do begin
            pos = strpos(key,'range')
            if pos[0] eq -1 then continue
            new_key = strmid(key,0,pos)+'extent'
            current_info[new_key] = total(current_info[key]*[-1,1])
;            if key eq 'mlt_range' then begin
;                current_info[new_key] = max(abs(current_info[key]))
;            endif
        endforeach


    ;---Velocity.
        triad_list = event.triad_list
        ntriad = triad_list.length
        v2d = fltarr(ntriad)
        foreach triad, triad_list, ii do begin
            v2d[ii] = triad.vmag_obs_time
        endforeach
        current_info['v2d_mag'] = mean(v2d)
        current_info['v2d_mag_stddev'] = stddev(v2d)

    ;---Omega.
        omega = fltarr(ntriad)
        foreach triad, triad_list, ii do begin
            omega[ii] = abs(triad.omega_obs_time)
        endforeach
        current_info['omega'] = mean(omega)
        current_info['omega_stddev'] = stddev(omega)

        event_info.add, current_info

    ;---Width.
        widths = fltarr(ndf)
        foreach df, df_list, ii do widths[ii] = df.width
        ;current_info['width_re'] = mean(widths)*current_info['v2d_mag']/constant('re')
        ;current_info['width_deg'] = mean(widths)*current_info['omega']/60
        current_info['width_deg'] = abs(current_info['omega'])*min(widths)/60
        current_info['width_re'] = current_info.width_deg*constant('rad')*10
        
    ;---Aspect ratio.
        rxy = current_info['rxy_range']
        r_center = mean(rxy)
        width_deg = current_info['width_deg']
        w_center = r_center*width_deg*constant('rad')
        asp = total(rxy*[-1,1])/w_center
        current_info['aspect_ratio'] = asp
    endforeach




    keys = ['time_extent','mlt_extent','rxy_extent','v2d_mag','v2d_mag_stddev','omega','width_re','width_deg','aspect_ratio']
    overall_info = dictionary()
    foreach key, keys do begin
        data = list()
        foreach info, event_info do data.add, info[key]
        overall_info[key] = data.toarray()
    endforeach

    keys = ['rxy_range']
    foreach key, keys do begin
        data = list()
        foreach info, event_info do data.add, info[key]
        overall_info[key] = data.toarray()
    endforeach



    nbin0 = 8

    data_info = list()
;    data_info.add, dictionary($
;        'label', 'Life time', $
;        'data', overall_info.time_extent/60, $
;        'xtitle', tex2str('tau')+' (min)')
;    data_info.add, dictionary($
;        'label', 'MLT extent', $
;        'data', overall_info.mlt_extent, $
;        'xtitle', tex2str('Delta')+'MLT (hr)')
;    data_info.add, dictionary($
;        'label', 'Radial extent', $
;        'data', overall_info.rxy_extent, $
;        'xtitle', tex2str('Delta')+'Rxy (Re)')
;    data_info.add, dictionary($
;        'label', 'Azimuthal speed', $
;        'data', overall_info.v2d_mag, $
;        'xtitle', '|v!D2D!N| (km/s)')
;    data_info.add, dictionary($
;        'label', 'Angular speed', $
;        'data', overall_info.omega, $
;        'xtitle', '|'+tex2str('omega')+'| (deg/min)')
;    data_info.add, dictionary($
;        'label', 'Azimuthal width', $
;        'data', overall_info.width_re, $
;        'xtitle', 'W (Re)')
    data_info.add, dictionary($
        'label', 'Length over width', $
        'data', overall_info.aspect_ratio, $
        'xtitle', 'L/W (#)')

    fig_xsize = 3
    fig_ysize = 3
    plot_file = join_path([project.plot_dir,'fig_df_shape_barplot.pdf'])
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize
    npanel = data_info.length
    nxpanel = 1
    nypanel = npanel/nxpanel
    margins = [8,5,2,2]
    poss = sgcalcpos(nypanel, nxpanel, xchsz=xchsz, ychsz=ychsz, margins=margins, ypad=6, xpad=5)

    xticklen_chsz = -0.30
    yticklen_chsz = -0.40

    fig_letters = letters(npanel)
    for xid=0,nxpanel-1 do begin
        for yid=0,nypanel-1 do begin
            tpos = poss[*,xid,yid]
            xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
            yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

            tid = yid*nxpanel+xid
            info = data_info[tid]
            data = info.data
            xtitle = info.xtitle
            label = info.label

            xrange = [floor(min(data)),ceil(max(data))]
            binsize = total(xrange*[-1,1])/nbin0
            bins = smkarthm(xrange[0], xrange[1], binsize, 'dx')
            nbin = n_elements(bins)-1

            yys = fltarr(nbin)
            for ii=0,nbin-1 do begin
                index = lazy_where(data, '[)', bins[ii:ii+1], count=count)
                yys[ii] = count
            endfor
            yys = yys/total(yys)*100


            xrange = xrange+[-1,1]*binsize*0.5
            yminor = 5
            ystep = 10.
            yrange = [0,ceil(abs(max(yys)/ystep)+0.1)*ystep]
            yticks = total(yrange*[-1,1])/ystep
            if yticks ge 5 then begin
                ystep = 20
                yrange = [0,ceil(abs(max(yys)/ystep)+0.1)*ystep]
                yticks = total(yrange*[-1,1])/ystep
            endif
            ytitle = (xid eq 0)? 'Percentage (%)': ' '

            xminor = 2
            plot, xrange, yrange, $
                xstyle=1, xrange=xrange, xticklen=xticklen, xtitle=xtitle, xminor=xminor, $
                ystyle=1, yrange=yrange, yticklen=yticklen, ytitle=ytitle, yminor=yminor, yticks=yticks, $
                /nodata, /noerase, position=tpos
            for ii=0,nbin-1 do begin
                plots, bins[ii]+[0,0], [0,yys[ii]]
                plots, bins[ii+1]+[0,0], [0,yys[ii]]
                plots, bins[ii:ii+1], yys[ii]+[0,0]
            endfor
            mean_value = mean(data)
            plots, mean_value+[0,0], yrange, linestyle=1
            tmp = convert_coord(mean_value,0, /data, /to_normal)
            tx = tmp[0]
            ty = tpos[3]-ychsz*0.8
            msg = sgnum2str(mean_value,nsgn=2)
            if tid eq 0 then msg = '    mean:'+msg else msg = '    '+msg
            unit = strmid(xtitle, strpos(xtitle,'('))
            msg += ' '+unit
            xyouts, tx,ty,/normal, msg, charsize=0.8


            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            ;xyouts, tx,ty,/normal, fig_letters[tid]+'. '+label
            xyouts, tx,ty,/normal, label
        endfor
    endfor

    if keyword_set(test) then stop
    sgclose



;    median_value = median(data)
;    plots, median_value+[0,0], yrange, linestyle=1


;    sgopen, 1, xsize=4, ysize=4
;    yys = overall_info.v2d_mag
;    yys_err = overall_info.v2d_mag_stddev
;    yys = overall_info.width_deg
;    yys_err = fltarr(n_elements(yys))
;    yrange = [floor(min(yys-yys_err)),ceil(max(yys+yys_err))]
;    xxs = overall_info.rxy_range
;    xrange = [floor(min(xxs)),ceil(max(xxs))]
;    xxs_err = 0.5*(xxs[*,1]-xxs[*,0])
;    xxs = 0.5*(xxs[*,1]+xxs[*,0])
;    xtitle = 'Rxy (Re)'
;    ytitle = 'v!D2D!N (km/s)'
;    xlog = 0
;    ylog = 0
;    plot, xrange, yrange, $
;        xstyle=1, xlog=xlog, xtitle=xtitle, xrange=xrange, xticklen=xticklen, $
;        ystyle=1, ylog=ylog, ytitle=ytitle, yrange=yrange, yticklen=yticklen, $
;        /nodata, /noerase, position=tpos
;    foreach tmp, xxs, ii do begin
;        plots, xxs[ii], yys[ii], psym=6
;        plots, xxs[ii]+[-1,1]*xxs_err[ii], yys[ii]+[0,0]
;        plots, xxs[ii]+[0,0], yys[ii]+[-1,1]*yys_err[ii]
;    endforeach

end
