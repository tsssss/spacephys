;+
; Generate low-res plots to survey on the data quality of the dE field generated in pflux_grant_preprocess_efield.
;-

pro pflux_grant_survey_on_efield, project=project, probe=probe, time_range=time_range, test=test


    if n_elements(project) eq 0 then project = pflux_grant_load_project()

    secofday = constant('secofday')
    plot_dir = join_path([project.plot_dir,'rbsp_efield_survey'])
    efield_time_step = 1d/16
    perigee_dis = 4.
    orbit_time_step = 1.

    xyz = ['y','z']
    colors = sgcolor(['blue','red'])
    e_small_range = [-1,1]*10.    ; mV/m.
    e_large_range = e_small_range*30
    e_lines = [25,50,100.]
    v_small_range = [-1,1]*10.   ; V.
    v_large_range = [-1,1]*200
    ndim = 2


    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    full_ychsz = constant('full_ychsz')
    half_ychsz = constant('half_ychsz')
    label_xshift = 6

    the_probe = probe
    prefix = 'rbsp'+the_probe+'_'
    rbspx = 'rbsp'+the_probe
    the_time_range = (n_elements(time_range) eq 2)? time_range: rbsp_info('emfisis_l3_data_range', probe=the_probe)
    ndate = total(the_time_range*[-1,1])/secofday
    dates = the_time_range[0]+findgen(ndate)*secofday
    paths = [default_local_root(),'sdata','rbsp',rbspx,'ebfield']

    foreach date, dates do begin
        date_time_range = date+[0,secofday]
        del_data, '*'

        ; Read orbit data.
        r_var = prefix+'r_gse'
        rbsp_read_orbit, date_time_range, probe=the_probe

        ; Read Vsc_median.
        rbsp_efw_read_boom_flag, date_time_range, probe=the_probe

        ; Read efield data.
        pflux_grant_read_e_mgse, date_time_range, probe=the_probe



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
        base_name = rbspx+'_efield_survey_'+time_string(date,tformat='YYYY_MMDD')+'_v01.pdf'
        plot_file = join_path([plot_dir,rbspx,year,base_name])
        if keyword_set(test) then plot_file = 0

        fig_xsize = 7.
        fig_ysize = 5.
        sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch, $
            xchsz=xchsz, ychsz=ychsz

        ypans = [1.,0.3,2]
        nypanel = n_elements(ypans)
        margins = [8,3,5,2]
        poss = sgcalcpos(nypanel, ypans=ypans, margins=margins)


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



    ;---Panel 1. Vsc.
        vsc_var = prefix+'vsc_median'
        vsc_pos = poss[*,0]
        the_poss = sgcalcpos(3, ypans=[1,2,1], position=vsc_pos, ypad=0)
        yys = get_var_data(vsc_var, times=xxs)

        ; Small part.
        tpos = the_poss[*,1]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = v_small_range
        ystep = 5
        ytickv = make_bins(yrange*0.5, ystep, /inner)
        yticks = n_elements(ytickv)
        yminor = 5
        ylog = 0
        ytitle = ' '
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, $
            ystyle=5, ylog=ylog, yrange=yrange, $
            position=tpos, /nodata, /noerase
        oplot, xxs, yys
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
        ytickv = min(yrange)*[1,10]
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

        ; Negative part.
        tpos = the_poss[*,2]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = reverse(abs([v_small_range[1],v_large_range[1]]))
        ylog = 1
        ytickv = min(yrange)*[10,1]
        yticks = n_elements(ytickv)-1
        yminor = 10
        ytickn = '-'+string(ytickv,format='(I0)')
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, $
            ystyle=5, ylog=ylog, yrange=yrange, $
            position=tpos, /nodata, /noerase
        oplot, xxs, -yys
        plot, xrange, yrange, $
            xstyle=9, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, ylog=ylog, yrange=yrange, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytickname=ytickn, $
            position=tpos, /nodata, /noerase

        ; Overall box.
        tpos = vsc_pos
        ty = (tpos[1]+tpos[3])*0.5-ychsz*half_ychsz
        tx = tpos[2]+xchsz*1
        xyouts, tx,ty,/normal, '|Vsc|'
        ty = (tpos[1]+tpos[3])*0.5
        tx = tpos[0]-xchsz*4
        xyouts, tx,ty,/normal, '(V)', orientation=90, alignment=0.5
        tx = tpos[0]-xchsz*label_xshift
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, 'a.'
        tx = tpos[0]
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,/normal, 'RBSP-'+strupcase(the_probe)+' dE (E after preprocessing) on '+time_string(date,tformat='YYYY-MM-DD')




    ;---Panel 2. flag.
        flag_pos = poss[*,1]
        tpos = flag_pos
        flag_var = prefix+'boom_flag'
        yys = total(get_var_data(flag_var, times=xxs),2) ne 4
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = [-0.2,1.2]
        ylog = 0
        ytickv = [0,1]
        yticks = n_elements(ytickv)
        yminor = 0
        ytickn = ['0','1']
        plot, xrange, yrange, $
            xstyle=5, xlog=xlog, xrange=xrange, $
            ystyle=5, ylog=ylog, yrange=yrange, $
            position=tpos, /nodata, /noerase
        oplot, xxs, yys
        plot, xrange, yrange, $
            xstyle=1, xlog=xlog, xrange=xrange, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, ylog=ylog, yrange=yrange, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytickname=ytickn, $
            position=tpos, /nodata, /noerase

        ; Overall box.
        tpos = flag_pos
        ty = (tpos[1]+tpos[3])*0.5-ychsz*half_ychsz
        tx = tpos[2]+xchsz*1
        xyouts, tx,ty,/normal, 'flag'
        ty = (tpos[1]+tpos[3])*0.5
        tx = tpos[0]-xchsz*4
        xyouts, tx,ty,/normal, '', orientation=90, alignment=0.5
        tx = tpos[0]-xchsz*label_xshift
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, 'b.'



    ;---Panel 3. E.
        de_pos = poss[*,2]
        the_poss = sgcalcpos(3, ypans=[1,1,1], position=de_pos, ypad=0)
        e_var = prefix+'e_mgse'
        yys = get_var_data(e_var, times=xxs)
        step = 16.
        xxs = xxs[0:*:step]
        yys = yys[0:*:step,1:2]

        boom_flags = total(get_var_data(prefix+'boom_flag', times=uts),2) ne 4
        index = where(boom_flags eq 1, count)
        pad_window = 600.   ; sec.
        if count ne 0 then begin
            bad_times = uts[time_to_range(index,time_step=1)]
            nbad_time = n_elements(bad_times)*0.5
            for section_id=0,nbad_time-1 do begin
                index = lazy_where(xxs, '[]', bad_times[section_id,*]+[-1,1]*pad_window, count=count)
                if count eq 0 then continue
                yys[index,*] = !values.f_nan
            endfor
        endif

        ; Small part.
        tpos = the_poss[*,1]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = e_small_range
        ystep = 5
        ytickv = make_bins(yrange*0.5, ystep, /inner)
        yticks = n_elements(ytickv)
        yminor = 5
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
        foreach ty, [yrange,yrange*0.5,0] do plots, xrange, ty+[0,0], linestyle=1

        ; Positive part.
        tpos = the_poss[*,0]
        xtickformat = '(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = abs([e_small_range[1],e_large_range[1]])
        ylog = 1
        ytickv = min(yrange)*[1,10]
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
        foreach ty, e_lines do plots, xrange, ty+[0,0], linestyle=1

        ; Negative part.
        tpos = the_poss[*,2]
        xtickformat = ''
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        yrange = reverse(abs([e_small_range[1],e_large_range[1]]))
        ylog = 1
        ytickv = min(yrange)*[10,1]
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
        foreach ty, e_lines do plots, xrange, ty+[0,0], linestyle=1

        ; Overall box.
        tpos = de_pos
        labels = 'dE!D'+xyz
        nlabel = n_elements(labels)
        dy = (tpos[3]-tpos[1])/(nlabel+1)
        tys = tpos[1]+(findgen(nlabel)+1)*dy-ychsz*half_ychsz
        tx = tpos[2]+xchsz*1
        foreach ty, tys, ii do xyouts, tx,ty,/normal, labels[ii], color=colors[ii]
        ty = (tpos[1]+tpos[3])*0.5
        tx = tpos[0]-xchsz*4
        xyouts, tx,ty,/normal, '(mV/m)', orientation=90, alignment=0.5
        tx = tpos[0]-xchsz*label_xshift
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, 'c.'

        if keyword_set(test) then stop
        sgclose
    endforeach
end


;probe = 'a'
;time_range = time_double(['2013-06-19','2015-10-01'])
;pflux_grant_survey_on_efield, project=project, probe=probe, $
;    time_range=time_range, test=1
;stop


time_range = time_double(['2012-10-01','2015-10-01'])
time_range = time_double(['2014-10-01','2015-10-01'])
foreach probe, ['b'] do $
    pflux_grant_survey_on_efield, project=project, probe=probe, time_range=time_range
end