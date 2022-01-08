;+
; Survey the whole mission for data quality of the spinfit E field.
;-

pro rbsp_survey_on_spinfit_efield, project=project, time_range=time_range

test = 1

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    probe_infos = project.probe_infos
    probes = probe_infos.keys()
    secofday = constant('secofday')
    plot_dir = join_path([project.plot_dir,'rbsp_efield_survey'])
    efield_resolution = 'survey'
    efield_nrec = secofday/11
    efield_min_nrec = 0.5*efield_nrec
    model = 't89'
    ndim = 3
    par = 2
    xyz = constant('xyz')
    colors = constant('rgb')
    small_range = [-1,1]*10.    ; mV/m.
    large_range = small_range*30

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

    log_file = join_path([project.data_dir,'rbsp_survey_on_spinfit_efield.log'])
    if file_test(log_file) eq 0 then ftouch, log_file

    foreach probe, probes do begin
        the_probe = strmid(probe,0,1,/reverse)
        prefix = 'rbsp'+the_probe+'_'

        the_time_range = (n_elements(time_range) eq 2)? time_range: rbsp_info('efw_l3_data_range', probe=the_probe)
        ndate = total(the_time_range*[-1,1])/secofday
        dates = smkarthm(the_time_range[0],secofday,ndate,'x0')

        foreach date, dates do begin
            date_time_range = date+[0,secofday]


            ; Read E data.
            e_var = prefix+'e_gsm'
            del_data, e_var
            rbsp_read_efield, date_time_range, probe=the_probe, resolution=efield_resolution
            get_data, e_var, times, e_gsm
            data_quality = 1
            if n_elements(e_gsm[*,0]) le efield_min_nrec then data_quality = 0
            index = where(finite(snorm(e_gsm)), count)
            if count le efield_min_nrec then data_quality = 0
            if data_quality eq 0 then begin
                msg = strupcase(probe)+ ': not enough good E data on '+time_string(date,tformat='YYYY/MM-DD')
                errmsg = handle_error(msg)
                lprmsg, msg, log_file
                continue
            endif

            emag = snorm(e_gsm)
            emag_var = prefix+'emag'
            store_data, emag_var, times, emag
            add_setting, emag_var, /smart, {$
                display_type: 'scalar', $
                unit: 'mV/m', $
                short_name: '|E|'}


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
            base_name = probe+'_efield_survey_'+time_string(date,tformat='YYYY_MMDD')+'_v01.pdf'
            plot_file = join_path([plot_dir,probe,year,base_name])
            if keyword_set(test) then plot_file = test

            sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch

        ;---Panel 1. |E|.
            dbmag_pos = poss[*,0]
            the_poss = sgcalcpos(2, ypans=[1,1], position=dbmag_pos, ypad=0)
            yys = get_var_data(emag_var, times=xxs)

            ; Small part.
            tpos = the_poss[*,1]
            xtickformat = '(A1)'
            xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
            yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

            yrange = small_range>0
            ystep = 4
            ytickv = make_bins(yrange, ystep, /inner)
            yticks = n_elements(ytickv)
            yminor = 2
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
            ytickv = large_range[1]*[0.1,1]
            ytickv = small_range[1]*[1,10]
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
            xyouts, tx,ty,/normal, '|E|'
            ty = (tpos[1]+tpos[3])*0.5
            tx = tpos[0]-xchsz*4
            xyouts, tx,ty,/normal, '(mV/m)', orientation=90, alignment=0.5
            tx = tpos[0]-xchsz*label_xshift
            ty = tpos[3]-ychsz*full_ychsz
            xyouts, tx,ty,/normal, 'a.'
            tx = tpos[0]
            ty = tpos[3]+ychsz*0.5
            xyouts, tx,ty,/normal, 'RBSP-'+strupcase(the_probe)+' E = E_0 in GSM on '+time_string(date,tformat='YYYY-MM-DD')


        ;---Panel 2. E.
            db_pos = poss[*,1]
            the_poss = sgcalcpos(3, ypans=[1,2,1], position=db_pos, ypad=0)
            yys = get_var_data(e_var, times=xxs)

            ; Small part.
            tpos = the_poss[*,1]
            xtickformat = '(A1)'
            xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
            yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

            yrange = small_range
            ystep = 4
            ytickv = make_bins(yrange, ystep, /inner)
            yticks = n_elements(ytickv)
            yminor = 2
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
            ytickv = small_range[1]*[1,10]
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
            ytickv = large_range[1]*[1,0.1]
            ytickv = small_range[1]*[10,1]
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
            labels = 'E!D'+xyz
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
            xyouts, tx,ty,/normal, 'b.'

            if keyword_set(test) then stop
            sgclose
        endforeach
    endforeach

end

rbsp_survey_on_spinfit_efield, project=project
end
