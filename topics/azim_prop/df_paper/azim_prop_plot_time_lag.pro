;+
; Based on azim_prop_calc_time_lag.pro. Here we make plots to show the before/after of the time lag.
;-

pro azim_prop_plot_time_lag, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    test = 0

;---Plot settings.
    margins = [12,3,2,2]
    fig_xsize = 6
    aspect_ratio = 0.3
    ytitle = '(deg)'
    xticklen = -0.06
    yticklen = xticklen*aspect_ratio
    half_ychsz = 0.35
    full_ychsz = 0.7
    line_skip = 0.35
    label_size = 0.7
    title_size = 1.0
    xpad = 10
    ypad = 0.2


    foreach event_id, events.keys() do begin
        event = events[event_id]
        data_file = event.file
        if file_test(data_file) eq 0 then azim_prop_gen_survey_plot, /save_data
        if file_test(data_file) eq 0 then message, 'Something is wrong when loading the data ...'
        tplot_restore, file = data_file

        cc_time_range = event.cc_time_range
        sorted_probes = event.sorted_probes


    ;---Plot settings.
        nvar = n_elements(sorted_probes)
        sgopen, 0, xsize=fig_xsize, ysize=fig_xsize
        tpos = sgcalcpos(1,2, xchsz=xchsz, ychsz=ychsz, ypad=ypad, xpad=xpad, margins=margins)
        sgclose, /wdelete
        pan_xsize = fig_xsize*(tpos[2,0]-tpos[0,0])
        pan_ysize = pan_xsize*aspect_ratio
        fig_ysize = pan_ysize*nvar+(ypad*(nvar-1)+(margins[1]+margins[3]))*ychsz*fig_xsize
        if keyword_set(test) then begin
            file = test
            magnify = 2
        endif else begin
            file = join_path([project.plot_dir,'time_lag','fig_db_tilt_'+event_id+'.pdf'])
            magnify = 1
        endelse

        thick = (size(file,/type) eq 7)? 4:2
        sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify
        poss = sgcalcpos(nvar,2, margins=margins, xchsz=xchsz,ychsz=ychsz, ypad=ypad, xpad=xpad)


    ;---Plot B tilt angle.
        str_theta = '!9'+string(113b)+'!X'
        str_delta = '!9'+string(68b)+'!X'
        plot_title = str_delta+str_theta+', tilt angle variation from T89'
        plot_title = 'Tilt angle in absolute time'
        event_time_range = events[event_id].time_range

        xrange = event_time_range
        xstep = 30*60
        xtickv = smkarthm(xrange[0],xrange[1],xstep, 'dx')
        xticks = n_elements(xtickv)-1
        xminor = 6
        while xticks ge 4 do begin
            xstep *= 2
            xtickv = smkarthm(xrange[0],xrange[1],xstep, 'dx')
            xticks = n_elements(xtickv)-1
            xminor *= 2
        endwhile
        xtickn = strarr(xticks+1)
        for ii=0, xticks do xtickn[ii] = time_string(xtickv[ii],tformat='hh:mm')
        xtickn[0] = time_string(xtickv[0],tformat='YYYY-MM-DD') & for ii=0, strlen(xtickn[0])-1 do xtickn[0] += ' '

        foreach probe, sorted_probes, ii do begin
            get_data, probe+'_db_tilt', times, dbtilt
            tpos = reform(poss[*,0,ii])
            xtickformat = (ii eq nvar-1)? '': '(A1)'

            yrange = sg_autolim(dbtilt)
            yrange = round(yrange/2)*2
            yticks = 2
            ytickv = smkarthm(yrange[0],yrange[1],yticks+1,'n')
            ytickn = strarr(yticks+1)
            yminor = 5
            for jj=0, yticks do ytickn[jj] = sgnum2str(ytickv[jj])
            plot, times, dbtilt, $
                xstyle=1, xrange=xrange, xtickformat=xtickformat, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, $
                ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytickname=ytickn, $
                position=tpos, /noerase, /nodata
            oplot, times, dbtilt, color=sgcolor('black')
            tx = event[probe].ref_time
            plots, tx+[0,0], yrange, linestyle=1
            ;times = times+(events[event_id])[probe].time_lag
            ;oplot, times, dbtilt, color=sgcolor('salmon'), linestyle=1
            ; Add label.
            ;tx = tpos[2]+xchsz*0.5
            ;ty = mean(tpos[[1,3]])-ychsz*half_ychsz
            ;xyouts, tx,ty,/normal, strupcase(project[probe].short_name)
            if ii eq 0 then begin
                tx = (tpos[0]+tpos[2])*0.5
                ty = tpos[3]+ychsz*line_skip
                xyouts, tx,ty,/normal,alignment=0.5, plot_title, charsize=title_size
                ; add the time range for cc.
                plot, event_time_range, [0,1], /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
                ty = 0.1
                txs = cc_time_range
                tys = ty+[0,0]
                plots, txs,tys
                tmp = convert_coord(mean(txs),ty, /data, /to_normal)
                xyouts, tmp[0], tmp[1]+0.2*ychsz*label_size, alignment=0.5, /normal, 'data for c.c', charsize=label_size

                for jj=0,1 do begin
                    tmp = convert_coord(txs[jj],tys[0], /data, /to_normal)
                    tx = tmp[0]
                    ty = tmp[1]
                    plots, tx+[0,0], ty+[-1,1]*ychsz*0.15, /normal
                endfor
            endif
            tx = tpos[0]-xchsz*6
            ty = tpos[3]-ychsz*half_ychsz
            xyouts, tx,ty,/normal, 'a-'+sgnum2str(ii+1)+'.'
            
            ; Add the error in c.c.
            if ii ne 0 then begin
                plot, event_time_range, [0,1], /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
                ty = 0.1
                cc_error = event[probe].time_lag_error
                txs = event[probe].ref_time+[-1,1]*cc_error
                tys = ty+[0,0]
                plots, txs,tys, thick=thick
                tmp = convert_coord(max(txs),ty, /data, /to_normal)
                if ii eq 1 then begin
                    xyouts, tmp[0]+xchsz*0.25, tmp[1]-0.2*ychsz*label_size, /normal, 'uncertainty of c.c', charsize=label_size
                endif

                for jj=0,1 do begin
                    tmp = convert_coord(txs[jj],tys[0], /data, /to_normal)
                    tx = tmp[0]
                    ty = tmp[1]
                    plots, tx+[0,0], ty+[-1,1]*ychsz*0.05, /normal
                endfor
            endif
            
            ; s/c name, MLT and R.
            tx = xchsz*5
            ty = mean(tpos[[1,3]])+ychsz*half_ychsz
            alignment = 1
            xyouts, tx,ty,/normal, alignment=alignment, strupcase(project[probe].short_name), color=project[probe].color
            xyouts, tx,ty-ychsz*label_size*1,/normal, alignment=alignment, sgnum2str(event[probe].ref_mlt,ndec=1)+' MLT', charsize=label_size
            ref_rsm = event[probe].ref_rsm
            xyouts, tx,ty-ychsz*label_size*2,/normal, alignment=alignment, sgnum2str(snorm(ref_rsm[0:1]),ndec=1)+' Re', charsize=label_size
            if ii eq 0 then begin
                tx = (tpos[0]+tpos[2])*0.5
                ty = tpos[3]+ychsz*line_skip
                xyouts, tx,ty,/normal,alignment=0.5, plot_title, charsize=title_size
            endif
        endforeach

        plot_title = 'Time lag removed'
        foreach probe, sorted_probes, ii do begin
            get_data, probe+'_db_tilt', times, dbtilt
            tpos = reform(poss[*,1,ii])
            xtickformat = (ii eq nvar-1)? '': '(A1)'

            yrange = sg_autolim(dbtilt)
            yrange = round(yrange/2)*2
            yticks = 2
            ytickv = smkarthm(yrange[0],yrange[1],yticks+1,'n')
            ytickn = strarr(yticks+1)
            yminor = 5
            for jj=0, yticks do ytickn[jj] = sgnum2str(ytickv[jj])
            plot, times, dbtilt, $
                xstyle=1, xrange=xrange, xtickformat=xtickformat, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, $
                ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytickname=ytickn, $
                position=tpos, /noerase, /nodata
            time_lag = event[probe].time_lag
            max_corr = event[probe].max_corr
            times = times+time_lag
            oplot, times, dbtilt, color=sgcolor('gray')
            tx = event.ref_time
            plots, tx+[0,0], yrange, linestyle=1, color=sgcolor('gray')
            ; Add label.
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*full_ychsz
            xyouts, tx,ty,/normal, str_delta+'t = '+sgnum2str(time_lag)+' sec', charsize=label_size
            if max_corr ne 0 then begin
                ;tx = tpos[0]+xchsz*8
                ty = tpos[3]-ychsz*full_ychsz*2
                xyouts, tx,ty,/normal, 'c.corr = '+sgnum2str(max_corr,ndec=2), charsize=label_size
            endif
            
            ; Add the error in c.c.
            if ii ne 0 then begin
                plot, event_time_range, [0,1], /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
                ty = 0.1
                cc_error = event[probe].time_lag_error
                txs = event.ref_time+[-1,1]*cc_error
                tys = ty+[0,0]
                plots, txs,tys, thick=thick
;                tmp = convert_coord(max(txs),ty, /data, /to_normal)
;                if ii eq 1 then begin
;                    xyouts, tmp[0]+xchsz*0.25, tmp[1]-0.25*ychsz*label_size, /normal, 'uncertainty of c.c', charsize=label_size
;                endif
    
                for jj=0,1 do begin
                    tmp = convert_coord(txs[jj],tys[0], /data, /to_normal)
                    tx = tmp[0]
                    ty = tmp[1]
                    plots, tx+[0,0], ty+[-1,1]*ychsz*0.05, /normal
                endfor
            endif
            
            tx = tpos[0]-xchsz*6
            ty = tpos[3]-ychsz*half_ychsz
            xyouts, tx,ty,/normal, 'b-'+sgnum2str(ii+1)+'.'
            
            if ii eq 0 then begin
                tx = (tpos[0]+tpos[2])*0.5
                ty = tpos[3]+ychsz*line_skip
                xyouts, tx,ty,/normal,alignment=0.5, plot_title, charsize=title_size
            endif
        endforeach

        if keyword_set(test) then stop
        sgclose
    endforeach

end
