;+
; Generate low-res plots to survey on the data quality of the dB field generated in pflux_grant_preprocess_bfield.
;-

pro pflux_grant_survey_on_bfield, project=project, probe=probe, time_range, $
    test=test, local_root=local_root

    if keyword_set(test) eq 0 then test = 0

    if n_elements(project) eq 0 then project = pflux_grant_load_project()

    secofday = constant('secofday')
    plot_dir = join_path([project.plot_dir,'rbsp_bfield_survey'])
    perigee_dis = 4.
    xyz = constant('xyz')
    colors = constant('rgb')
    small_range = [-1,1]*50.    ; nT.
    large_range = small_range*10
    ndim = 3


    fig_xsize = 7.
    fig_ysize = 5.
    ypans = [1.,2]
    nypanel = n_elements(ypans)
    margins = [8,3,5,2]
    sgopen, 0, xsize=fig_xsize, ysize=fig_ysize, /inch
    poss = sgcalcpos(nypanel, ypans=ypans, margins=margins, xchsz=xchsz, ychsz=ychsz)
    sgclose, /wdelete

    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    full_ychsz = constant('full_ychsz')
    half_ychsz = constant('half_ychsz')
    label_xshift = 6

    the_probe = probe
    prefix = 'rbsp'+the_probe+'_'
    rbspx = 'rbsp'+the_probe
    ndate = total(time_range*[-1,1])/secofday
    dates = time_range[0]+findgen(ndate)*secofday

    foreach date, dates do begin
        date_time_range = date+[0,secofday]
        del_data, '*'

        ; Read B1 and B0.
        pflux_grant_read_bfield, date_time_range, probe=probe, local_root=local_root
        db_var = prefix+'b1_gsm'
        uniform_time, db_var, 1     ; Down sample.
        db_gsm = get_var_data(db_var, times=times)
        dbmag = snorm(db_gsm)
        dbmag_var = prefix+'dbmag'
        store_data, dbmag_var, times, dbmag
        add_setting, dbmag_var, /smart, {$
            display_type: 'scalar', $
            unit: 'nT', $
            short_name: '|dB|'}


        ; Read orbit data.
        r_gsm_var = prefix+'r_gsm'
        r_gse_var = prefix+'r_gse'
        rbsp_read_orbit, date_time_range, probe=the_probe
        r_gse = get_var_data(r_gse_var, times=orbit_times)
        r_gsm = cotran(r_gse, orbit_times, 'gse2gsm')
        store_data, r_gsm_var, orbit_times, r_gsm
        add_setting, r_gsm_var, /smart, dictionary($
            'display_type', 'vector', $
            'short_name', 'R', $
            'unit', 'Re', $
            'coord', 'GSM', $
            'coord_labels', xyz )

        dis_var = prefix+'dis'
        dis = snorm(r_gsm)
        store_data, dis_var, orbit_times, dis
        add_setting, dis_var, /smart, dictionary($
            'display_type', 'scalar', $
            'short_name', '|R|', $
            'unit', 'Re' )


    ;---Settings for plot.
        ; Common x-axis.
        xrange = date_time_range
        xstep = constant('secofhour')*4
        xminor = 8
        xtickv = make_bins(xrange, xstep)
        xticks = n_elements(xtickv)-1
        xtickn = strarr(xticks+1)
        for ii=0, xticks do xtickn[ii] = time_string(xtickv[ii],tformat='hh:mm')
        xtickn[0] = time_string(date,tformat='YYYY-MM-DD')+'      '

        year = time_string(date,tformat='YYYY')
        base_name = rbspx+'_dbfield_survey_'+time_string(date,tformat='YYYY_MMDD')+'_v01.pdf'
        plot_file = join_path([plot_dir,rbspx,year,base_name])
        if keyword_set(test) then plot_file = test

        sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch


    ;---Perigee times.
        r_var = prefix+'r_gsm'
        dis = snorm(get_var_data(r_var, times=times))
        orbit_time_step = total(times[0:1]*[-1,1])
        index = where(dis le perigee_dis)
        perigee_time_ranges = time_to_range(times[index], time_step=orbit_time_step)
        color = sgcolor('red')
        tpos = poss[*,0]
        tpos[1] = min(poss[1,*])
        tpos[3] = max(poss[3,*])
        yrange = [0,1]
        plot, xrange, yrange, /nodata, /noerase, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos
        foreach time, perigee_time_ranges do begin
            oplot, time+[0,0], yrange, linestyle=1, color=color
        endforeach

    ;---Panel 1. |dB|.
        dbmag_pos = poss[*,0]
        the_poss = sgcalcpos(2, ypans=[1,1], position=dbmag_pos, ypad=0)
        yys = get_var_data(dbmag_var, times=xxs)

        ; Small part.
        tpos = the_poss[*,1]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = small_range>0
        ystep = 20
        ytickv = make_bins(yrange, ystep, /inner)
        yticks = n_elements(ytickv)
        yminor = 4
        ylog = 0
        ytitle = ' '
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, $
            ystyle=5, ylog=ylog, yrange=yrange, $
            position=tpos, /nodata, /noerase
        oplot, xxs, yys
        plot, xrange, yrange, $
            xstyle=9, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, ylog=ylog, yrange=yrange, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, $
            position=tpos, /nodata, /noerase
        plots, xrange, yrange[1]+[0,0], linestyle=1

        ; Large part.
        tpos = the_poss[*,0]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = abs([small_range[1],large_range[1]])
        ylog = 1
        ytickv = yrange
        yticks = n_elements(ytickv)-1
        yminor = 10
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, $
            ystyle=5, ylog=ylog, yrange=yrange, $
            position=tpos, /nodata, /noerase
        oplot, xxs, yys
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, ylog=ylog, yrange=yrange, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, $
            position=tpos, /nodata, /noerase
        axis, xaxis=1, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat

        ; Overall box.
        tpos = dbmag_pos
        ty = (tpos[1]+tpos[3])*0.5-ychsz*half_ychsz
        tx = tpos[2]+xchsz*1
        xyouts, tx,ty,/normal, '|dB|'
        ty = (tpos[1]+tpos[3])*0.5
        tx = tpos[0]-xchsz*4
        xyouts, tx,ty,/normal, '(nT)', orientation=90, alignment=0.5
        tx = tpos[0]-xchsz*label_xshift
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, 'a.'
        tx = tpos[0]
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,/normal, 'RBSP-'+strupcase(the_probe)+' dB (B1 after preprocessing) on '+time_string(date,tformat='YYYY-MM-DD')


    ;---Panel 2. dB.
        db_pos = poss[*,1]
        the_poss = sgcalcpos(3, ypans=[1,2,1], position=db_pos, ypad=0)
        yys = get_var_data(db_var, times=xxs)

        ; Small part.
        tpos = the_poss[*,1]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = small_range
        ystep = 20
        ytickv = make_bins(yrange, ystep, /inner)
        yticks = n_elements(ytickv)
        yminor = 4
        ylog = 0
        ytitle = ' '
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, $
            ystyle=5, ylog=ylog, yrange=yrange, $
            position=tpos, /nodata, /noerase
        for ii=0, ndim-1 do oplot, xxs, yys[*,ii], color=colors[ii]
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, ylog=ylog, yrange=yrange, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, $
            position=tpos, /nodata, /noerase
        foreach ty, [yrange,0] do plots, xrange, ty+[0,0], linestyle=1

        ; Positive part.
        tpos = the_poss[*,0]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = abs([small_range[1],large_range[1]])
        ylog = 1
        ytickv = yrange
        yticks = n_elements(ytickv)-1
        yminor = 10
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, $
            ystyle=5, ylog=ylog, yrange=yrange, $
            position=tpos, /nodata, /noerase
        for ii=0, ndim-1 do oplot, xxs, yys[*,ii], color=colors[ii]
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, ylog=ylog, yrange=yrange, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, $
            position=tpos, /nodata, /noerase
        axis, xaxis=1, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat

        ; Negative part.
        tpos = the_poss[*,2]
        xtickformat = ''
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = reverse(abs([small_range[1],large_range[1]]))
        ylog = 1
        ytickv = yrange
        yticks = n_elements(ytickv)-1
        yminor = 10
        ytickn = '-'+string(ytickv,format='(I0)')
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, $
            ystyle=5, ylog=ylog, yrange=yrange, $
            position=tpos, /nodata, /noerase
        for ii=0, ndim-1 do oplot, xxs, -yys[*,ii], color=colors[ii]
        plot, xrange, yrange, $
            xstyle=9, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, ylog=ylog, yrange=yrange, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytickname=ytickn, $
            position=tpos, /nodata, /noerase

        ; Overall box.
        tpos = db_pos
        labels = 'dB!D'+xyz
        nlabel = n_elements(labels)
        dy = (tpos[3]-tpos[1])/(nlabel+1)
        tys = tpos[1]+(findgen(nlabel)+1)*dy-ychsz*half_ychsz
        tx = tpos[2]+xchsz*1
        foreach ty, tys, ii do xyouts, tx,ty,/normal, labels[ii], color=colors[ii]
        ty = (tpos[1]+tpos[3])*0.5
        tx = tpos[0]-xchsz*4
        xyouts, tx,ty,/normal, '(nT)', orientation=90, alignment=0.5
        tx = tpos[0]-xchsz*label_xshift
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, 'b.'

        if keyword_set(test) then stop
        sgclose

    endforeach
end



;probe = 'a'
;time_range = time_double(['2015-03-01','2015-10-01'])
;pflux_grant_survey_on_bfield, project=project, probe=probe, time_range


probe = 'a'
time_range = time_double(['2012-10-01','2015-10-01'])
pflux_grant_survey_on_bfield, project=project, probe=probe, time_range
end