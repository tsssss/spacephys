;+
; Generate low-res plots to survey on the data quality of the pflux generated in pflux_grant_preprocess_pflux.
;-

pro pflux_grant_survey_on_pflux, project=project, probe=probe, time_range=time_range

test = 0

    if n_elements(project) eq 0 then project = pflux_grant_load_project()

    secofday = constant('secofday')
    plot_dir = join_path([project.plot_dir,'rbsp_pflux_survey'])
    efield_time_step = 1d/16
    perigee_dis = 4.
    orbit_time_step = 1.

    pf_small_range = [-1,1]*5.      ; mW/m^2.
    pf_large_range = pf_small_range*200
    pf_lines = [50,500.]
    v_small_range = [-1,1]*5.      ; mW/m^2.
    v_large_range = v_small_range*200
    small_pf = 2.

    fig_xsize = 7.
    fig_ysize = 5.
    ypans = [1.,2]
    nypanel = n_elements(ypans)

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
    paths = [default_local_root(),'sdata','rbsp',rbspx,'ebfield']

    foreach date, dates do begin
        date_time_range = date+[0,secofday]
        del_data, '*'

        ; Read orbit data.
        r_var = prefix+'r_gse'
        rbsp_read_orbit, date_time_range, probe=the_probe

        ; Read pflux data.
        year = time_string(date,tformat='YYYY')
        pflux_grant_read_pflux, date_time_range, probe=probe
        pf_var = prefix+'pf_fac_norm'
        lowres_times = make_bins(date_time_range,orbit_time_step)
        interp_time, pf_var, lowres_times

        pflux_grant_read_efield, date_time_range, probe=probe
        e_var = prefix+'e_mgse'
        get_data, e_var, times, e_mgse
        index = where(finite(snorm(e_mgse),/nan), count)
        if count ne 0 then begin
            nan_time_ranges = times[time_to_range(index,time_step=1)]
            durations = nan_time_ranges[*,1]-nan_time_ranges[*,0]
            index = where(durations ge 600., nnan_time_range)
            if nnan_time_range eq 0 then continue
            nan_time_ranges = nan_time_ranges[index,*]
            pf_fac = get_var_data(pf_var)
            fillval = !values.f_nan
            pad_time = 1800.
            for ii=0,nnan_time_range-1 do begin
                index = lazy_where(lowres_times, '[]', nan_time_ranges[ii,*]+[-1,1]*pad_time, count=count)
                if count eq 0 then continue
                pf_fac[index,*] = fillval
            endfor
            store_data, pf_var, lowres_times, pf_fac
        endif




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
        base_name = rbspx+'_pflux_survey_'+time_string(date,tformat='YYYY_MMDD')+'_v01.pdf'
        plot_file = join_path([plot_dir,rbspx,year,base_name])
        if keyword_set(test) then plot_file = test

        sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch
        margins = [8,3,5,2]
        poss = sgcalcpos(nypanel, ypans=ypans, margins=margins, xchsz=xchsz, ychsz=ychsz)


    ;---Perigee times.
        dis = snorm(get_var_data(r_var, times=times))
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

    ;---Panel 1. |S_para|/|S|.
        vsc_pos = poss[*,0]
        the_poss = sgcalcpos(3, ypans=[1,2,1], position=vsc_pos, ypad=0)
        pfs = get_var_data(pf_var, times=xxs)
        yys = [[pfs[*,0]],[snorm(pfs)]]
        xyz = ['S!Dearth','|S|']
        colors = sgcolor(['red','black'])
        ndim = n_elements(colors)

        ; Small part.
        tpos = the_poss[*,1]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = v_small_range
        ystep = 4
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

        yrange = abs([v_small_range[1],v_large_range[1]])
        ylog = 1
        ytickv = min(yrange)*[1,10,100]
        yticks = n_elements(ytickv)-1
        yminor = 10
        ytickn = string(ytickv,format='(I0)')
        ytickn[1] = ' '
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, $
            ystyle=5, ylog=ylog, yrange=yrange, $
            position=tpos, /nodata, /noerase
        for ii=0, ndim-1 do oplot, xxs, yys[*,ii], color=colors[ii]
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, ylog=ylog, yrange=yrange, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytickn=ytickn, $
            position=tpos, /nodata, /noerase
        axis, xaxis=1, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat

        ; Negative part.
        tpos = the_poss[*,2]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = reverse(abs([v_small_range[1],v_large_range[1]]))
        ylog = 1
        ytickv = min(yrange)*[100,10,1]
        yticks = n_elements(ytickv)-1
        yminor = 10
        ytickn = '-'+string(ytickv,format='(I0)')
        ytickn[1] = ' '
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
        tpos = vsc_pos
        labels = xyz
        nlabel = n_elements(labels)
        dy = (tpos[3]-tpos[1])/(nlabel+1)
        tys = tpos[1]+(findgen(nlabel)+1)*dy-ychsz*half_ychsz
        tx = tpos[2]+xchsz*1
        foreach ty, reverse(tys), ii do xyouts, tx,ty,/normal, labels[ii], color=colors[ii]
        ty = (tpos[1]+tpos[3])*0.5
        tx = tpos[0]-xchsz*4
        xyouts, tx,ty,/normal, '(mW/m!U2!N)', orientation=90, alignment=0.5
        tx = tpos[0]-xchsz*label_xshift
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, 'a.'
        tx = tpos[0]
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,/normal, 'RBSP-'+strupcase(the_probe)+' S (preprocessed ExB, normalized to 100 km) on '+time_string(date,tformat='YYYY-MM-DD')


    ;---Panel 2. S.
        de_pos = poss[*,1]
        the_poss = sgcalcpos(3, ypans=[1,1,1], position=de_pos, ypad=0)
        yys = pfs
        xyz = ['earth','west','out']
        colors = constant('rgb')
        ndim = n_elements(colors)

        ; Small part.
        tpos = the_poss[*,1]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = pf_small_range
        ystep = 4
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
        foreach ty, [yrange,[-1,1]*small_pf,0] do plots, xrange, ty+[0,0], linestyle=1

        ; Positive part.
        tpos = the_poss[*,0]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = abs([pf_small_range[1],pf_large_range[1]])
        ylog = 1
        ytickv = min(yrange)*[1,10,100]
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
        foreach ty, pf_lines do plots, xrange, ty+[0,0], linestyle=1

        ; Negative part.
        tpos = the_poss[*,2]
        xtickformat = ''
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = reverse(abs([pf_small_range[1],pf_large_range[1]]))
        ylog = 1
        ytickv = min(yrange)*[100,10,1]
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
        foreach ty, pf_lines do plots, xrange, ty+[0,0], linestyle=1

        ; Overall box.
        tpos = de_pos
        labels = 'S!D'+xyz
        nlabel = n_elements(labels)
        dy = (tpos[3]-tpos[1])/(nlabel+1)
        tys = tpos[1]+(findgen(nlabel)+1)*dy-ychsz*half_ychsz
        tx = tpos[2]+xchsz*1
        foreach ty, reverse(tys), ii do xyouts, tx,ty,/normal, labels[ii], color=colors[ii]
        ty = (tpos[1]+tpos[3])*0.5
        tx = tpos[0]-xchsz*4
        xyouts, tx,ty,/normal, '(mW/m!U2!N)', orientation=90, alignment=0.5
        tx = tpos[0]-xchsz*label_xshift
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, 'b.'

        if keyword_set(test) then stop
        sgclose


    endforeach
end


;probe = 'a'
;time_range = time_double(['2012-10-01','2015-10-01'])
;pflux_grant_survey_on_pflux, project=project, probe=probe, time_range=time_range

probe = 'b'
time_range = time_double(['2012-10-01','2015-10-01'])
;time_range = time_double(['2013-06-07','2013-06-08'])
pflux_grant_survey_on_pflux, project=project, probe=probe, time_range=time_range
end
