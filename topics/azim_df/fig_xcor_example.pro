;+
; Based on the results of azim_df_search_candidate and azim_df_search_event. We only use the results from search_large_df because search_df_group is not good enough and its job is done here in a more complicated way.
;
; Search final events and write results to a text file that is in a different
; format from azim_df_search_candidate and azim_df_search_event
;-

test = 0

    test_time = time_double('2014-08-28/10:00')    ; should get thd,the,g15,tha,rbspb,g13

;---Check inputs.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(time_range) ne 2 then begin
        errmsg = handle_error('No input time_range ...')
        return
    endif


;---Find the event.
    events = azim_df_search_all_events(project=project)
    foreach event, events do if product(event.time_range-short_time_range) lt 0 then break
    azim_df_load_basic_data, project=project


;---Settings for DF group.
    xcor_section_ratio = [-1,1]*2.5
    common_data_rate = project.time_step

    probe_infos = project.probe_infos

    df_list = event.df_list
    ndf = n_elements(df_list)

    plot_file = join_path([project.plot_dir,'fig_xcor_example.pdf'])
    if keyword_set(test) then plot_file = 0
    blue = sgcolor('blue')
    red = sgcolor('red')
    label_size = constant('label_size')
    full_ychsz = constant('full_ychsz')
    xticklen_chsz = -0.3
    yticklen_chsz = -0.5
    sgopen, plot_file, xsize=5, ysize=7, xchsz=xchsz, ychsz=ychsz
    margins = [10,5,5,1]
    poss = sgcalcpos(ndf-1, 2, xpad=10, xpans=[2,1], ypad=5, margins=margins)
    fig_labels = letters(ndf)

    foreach the_df, df_list, jj do begin
        if jj eq 0 then continue
        pre_df = df_list[jj-1]
        pre_df_width = pre_df.width
        the_df_width = the_df.width
        the_width = min([pre_df_width,the_df_width])
        the_time_range = the_width*xcor_section_ratio
        nlag = floor(the_width*2/common_data_rate)
        section_times = make_bins(the_time_range, common_data_rate)

        pre_time = pre_df.obs_time+section_times
        pre_data = get_var_data(pre_df.probe+'_theta', at=pre_time)
        the_time = the_df.obs_time+section_times
        the_data = get_var_data(the_df.probe+'_theta', at=the_time)
        adjust_time_lag = the_df.obs_time-pre_df.obs_time

        nlag = floor(the_width*2/common_data_rate)
        lags = findgen(nlag)-round(nlag/2)
        xcors = c_correlate(pre_data, the_data, lags)
        xcor_max = round(max(xcors, index)*10)/10.
        xcor_time_lag = lags[index]*common_data_rate
        if xcor_max gt 0 and finite(xcor_max) then begin
            xx = pre_data
            del_x = xx-mean(xx)
            dx_dt = deriv(xx)/common_data_rate
            nrec = n_elements(xx)
            xcor_err = sqrt(1./(nrec-1)*(1-xcor_max)/xcor_max*2*mean(del_x^2)/mean(dx_dt^2))
        endif else begin
            xcor_max = 2.
            xcor_err = 0.
            xcor_time_lag = 0.
        endelse
        xcor_time = pre_df.obs_time+xcor_time_lag+adjust_time_lag
        time_lag = xcor_time_lag+adjust_time_lag


    ;---Panel 1: theta.
        tpos = poss[*,0,jj-1]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        xxs = pre_time-pre_df.obs_time
        xrange = minmax(xxs)
        xtitle = 'Second from '+time_string(pre_df.obs_time,tformat='YYYY-MM-DD/hh:mm:ss')
        xticks = 2
        yticks = 2
        yminor = 5
        ytitle = '(deg)'
        yrange = ceil(max(abs(pre_data))/2)*2*[-1,1]
        plot, xxs, pre_data, $
            xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xtitle=xtitle, xticklen=xticklen, $
            ystyle=9, yrange=yrange, yticks=yticks, yminor=yminor, ytitle=ytitle, yticklen=yticklen, $
            position=tpos, /noerase
        plots, xrange, [0,0], linestyle=1

        yrange = ceil(max(abs(the_data))/2)*2*[-1,1]
        plot, the_time-xcor_time, the_data, color=red, $
            xstyle=5, $
            ystyle=5, yrange=yrange, $
            position=tpos, /noerase
        axis, yaxis=1, yrange=yrange, ystyle=1, yticks=yticks, yminor=yminor, ytitle=' ', yticklen=yticklen, color=red

        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = strupcase(probe_infos[pre_df.probe].short_name)
        xyouts, tx,ty,/normal, msg
        tx = tpos[0]+xchsz*4.5
        msg = strupcase(probe_infos[the_df.probe].short_name)
        xyouts, tx,ty,/normal, msg, color=red

        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, fig_labels[jj-1]+'-1.'


    ;---Panel 2: time lag.
        tpos = poss[*,1,jj-1]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        xxs = lags*common_data_rate+adjust_time_lag
        yys = xcors
        xrange = minmax(xxs)
        xminor = 5
        xticks = 2
        xtitle = 'Time lag (sec)'
        yrange = [0.,1]
        yticks = 2
        yminor = 5
        ytitle = 'xcor'
        plot, xxs, yys, $
            xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xticklen=xticklen, xtitle=xtitle, $
            ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            position=tpos, /noerase

        the_x = xcor_time_lag+adjust_time_lag
        the_y = max(yys)
        plots, the_x+[0,0], yrange, linestyle=1
        plots, the_x,the_y, color=blue, psym=6, symsize=0.5
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.5
        msg = 't_xcor = '+string(the_x,format='(I0)')+' sec'
        xyouts, tx,ty, /normal, msg, charsize=label_size, color=blue
        ty = tpos[1]+ychsz*(0.5+label_size*2)
        msg = 'Max xcor = '+string(the_y,format='(F3.1)')
        xyouts, tx,ty, /normal, msg, charsize=label_size, color=blue
        ty = tpos[1]+ychsz*(0.5+label_size*1)
        msg = 'Error = '+string(xcor_err,format='(I0)')+' sec'
        xyouts, tx,ty, /normal, msg, charsize=label_size, color=blue

        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*full_ychsz
        xyouts, tx,ty,/normal, fig_labels[jj-1]+'-2.'
    endforeach

    if keyword_set(test) then stop
    sgclose
end
