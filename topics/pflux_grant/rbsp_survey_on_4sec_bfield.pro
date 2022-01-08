;+
; Survey on data quality of the 4-sec B field.
;-

pro rbsp_survey_on_4sec_bfield, project=project, time_range=time_range

test = 0

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    probe_infos = project.probe_infos
    probes = probe_infos.keys()
    secofday = constant('secofday')
    plot_dir = join_path([project.plot_dir,'rbsp_bfield_survey'])
    bfield_resolution = '4sec'
    bfield_nrec = secofday/4.
    bfield_min_nrec = 0.5*bfield_nrec
    model = 't89'
    ndim = 3
    par = 2
    xyz = constant('xyz')
    colors = constant('rgb')
    small_range = [-1,1]*50.    ; nT.
    large_range = small_range*10

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
    
    log_file = join_path([project.data_dir,'rbsp_survey_on_4sec_bfield.log'])
    ;if file_test(log_file) eq 1 then file_delete, log_file
    if file_test(log_file) eq 0 then ftouch, log_file

    foreach probe, probes do begin
        the_probe = strmid(probe,0,1,/reverse)
        prefix = 'rbsp'+the_probe+'_'

        the_time_range = (n_elements(time_range) eq 2)? time_range: rbsp_info('emfisis_l3_data_range', probe=the_probe)      
        ndate = total(the_time_range*[-1,1])/secofday
        dates = smkarthm(the_time_range[0],secofday,ndate,'x0')

        foreach date, dates do begin
            date_time_range = date+[0,secofday]


            ; Read B data.
            b_var = prefix+'b_gsm'
            del_data, b_var
            rbsp_read_bfield, date_time_range, probe=the_probe, resolution=bfield_resolution
            get_data, b_var, times, b_gsm
            if n_elements(b_gsm[*,0]) le bfield_min_nrec then begin
                msg = 'Not enough good data on '+time_string(date,tformat='YYYY/MM-DD')
                errmsg = handle_error(msg)
                lprmsg, msg, log_file
                continue
            endif

            ; Read orbit data.
            r_var = prefix+'r_gsm'
            rbsp_read_orbit, date_time_range, probe=the_probe
            r_gsm = get_var_data(r_var, at=times)
            ntime = n_elements(times)
            b0_gsm = fltarr(ntime,ndim)
            for ii=0, ntime-1 do begin
                tilt = geopack_recalc(times[ii])
                ; in-situ position
                rx = r_gsm[ii,0]
                ry = r_gsm[ii,1]
                rz = r_gsm[ii,2]
                ; in-situ B field.
                geopack_igrf_gsm, rx,ry,rz, bx,by,bz
                geopack_t89, par, rx,ry,rz, dbx,dby,dbz
                b0_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
            endfor
            b0_var = prefix+'b0_gsm'
            store_data, b0_var, times, b0_gsm
            add_setting, b0_var, /smart, {$
                display_type: 'vector', $
                unit: 'nT', $
                short_name: 'B0!S!UT89!N!R', $
                coord: 'GSM', $
                coord_labels: xyz}

            db_gsm = b_gsm-b0_gsm
            db_var = prefix+'db_gsm'
            store_data, db_var, times, db_gsm
            add_setting, db_var, /smart, {$
                display_type: 'vector', $
                unit: 'nT', $
                short_name: 'dB', $
                coord: 'GSM', $
                coord_labels: xyz}

            dbmag = snorm(db_gsm)
            dbmag_var = prefix+'dbmag'
            store_data, dbmag_var, times, dbmag
            add_setting, dbmag_var, /smart, {$
                display_type: 'scalar', $
                unit: 'nT', $
                short_name: '|dB|'}


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
            base_name = probe+'_bfield_survey_'+time_string(date,tformat='YYYY_MMDD')+'_v01.pdf'
            plot_file = join_path([plot_dir,probe,year,base_name])
            if keyword_set(test) then plot_file = test

            sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch

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
            xyouts, tx,ty,/normal, 'RBSP-'+strupcase(the_probe)+' dB = B-B!DT89!N in GSM on '+time_string(date,tformat='YYYY-MM-DD')


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
    endforeach

end

rbsp_survey_on_4sec_bfield, project=project
end
