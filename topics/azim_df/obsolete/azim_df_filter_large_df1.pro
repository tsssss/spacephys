;+
; Check if a given DF is a large DF.
;-

function azim_df_filter_large_df1, df_list, project=project

    retval = dictionary()
    log_file = -1
    tab = constant('4space')
    if n_elements(project) eq 0 then project = azim_df_load_project()
    time_step = project.time_step


;---Settings.
    clear_window_size = 60*10. ; sec.
    height_range = [0.,80]      ; deg.
    width_range = [1.,25]*60.   ; sec.
    min_scaled_height = 8.      ; deg.
    scale_width = 4.            ; MLT.
    scale_width = 4.            ; MLT.
    section_time_range = [0.,10]*60 ; sec.
    section_min_scaled_height = 2.  ; deg.
    section_min_ratio = 0.6         ; #.


    ; Derived quantities.
    section_duration = total(section_time_range*[-1,1])
    section_count = section_duration/time_step



;---Collect info.
    obs_times = list()
    foreach df, df_list do obs_times.add, df.obs_time
    event_time_range = minmax(obs_times.toarray())

    probe_list = list()
    foreach df, df_list do begin
        probe = df.probe
        if probe_list.where(probe) eq !null then probe_list.add, probe
    endforeach


;---Load data.
    data_time_range = event_time_range+section_time_range
    foreach probe, probe_list do begin
        prefix = probe+'_'
        theta_var = prefix+'scaled_theta'
        if check_if_update(theta_var, data_time_range) then begin
            azim_df_read_data, 'theta', time_range=data_time_range, probe=probe, project=project
            azim_df_read_data, 'r_sm', time_range=data_time_range, probe=probe, project=project
            theta = get_var_data(prefix+'theta', times=times)
            r_sm = get_var_data(prefix+'r_sm', at=times)
            mlt = pseudo_mlt(r_sm)
            store_data, prefix+'pseudo_mlt', times, mlt
            store_data, theta_var, times, azim_df_scale_theta(theta, mlt, width=scale_width)
        endif
    endforeach


;---Ensure scaled_height is up to date.
    foreach df, df_list do df.scaled_height = azim_df_scale_theta(df.height, df.obs_mlt)


;---Sort DF according to probe.
    df_dict = dictionary()
    foreach df, df_list do begin
        probe = df.probe
        if ~df_dict.haskey(probe) then df_dict[probe] = list()
        df_dict[probe].add, df
    endforeach


;---For each probe, filter large DFs.
    foreach probe, probe_list do begin
        the_df_list = df_dict[probe]
        flags = bytarr(the_df_list.length)+1
        foreach df, the_df_list, df_id do begin

;if probe eq 'the' and df.obs_time-time_double('2007-11-20/17:30:50') eq 0 then stop
;if probe eq 'thc' and df.obs_time-time_double('2007-11-20/17:38:20') eq 0 then stop
if probe eq 'the' and df.obs_time-time_double('2007-12-20/09:10:45') eq 0 then stop


        ;---Check clear window.
            if df_id lt the_df_list.length-1 then begin
                dtime = the_df_list[df_id+1].obs_time-df.obs_time
                if dtime lt clear_window_size then begin
                    flags[df_id] = 0
                    continue
                endif
            endif

        ;---Check width.
            index = where_pro(df.width, '][', width_range, count=count)
            if count ne 0 then begin
                flags[df_id] = 0
                continue
            endif

        ;---Check height.
            index = where_pro(df.height, '][', height_range, count=count)
            if count ne 0 then begin
                flags[df_id] = 0
                continue
            endif

        ;---Check scaled_height.
            if round(df.scaled_height) le min_scaled_height then begin
                flags[df_id] = 0
                continue
            endif

        ;---Check shape.
            prefix = probe+'_'
            theta_var = prefix+'scaled_theta'
            the_time_range = df.obs_time+section_time_range
            scaled_heights = get_var_data(theta_var, in=the_time_range)
            index = where(scaled_heights ge section_min_scaled_height, count)
            section_ratio = round(float(count)/section_count*10)/10.
            azim_df_write_df, df, filename=log_file
            print, mean(scaled_heights)
            print, stddev(scaled_heights)
            msg = 'Check scaled_height: >'+string(section_min_scaled_height,format='(I0)')+$
                ' deg for '+string(section_ratio,format='(F3.1)')+' of ['+$
                strjoin(string(section_time_range,format='(I0)'),',')+'] sec around obs_time ...'
            lprmsg, msg, log_file
            if section_ratio lt section_min_ratio then begin
                flags[df_id] = 0
                continue
            endif
        endforeach

        index = where(flags eq 1, count)
        if count eq 0 then begin
            df_dict.remove, probe
        endif else begin
            df_dict[probe] = the_df_list[index]
        endelse
    endforeach

    large_df_list = list()
    foreach the_df_list, df_dict do large_df_list.add, the_df_list, /extract
    return, large_df_list

end




;---The wanted event.
    event_list = list()

    ; post_midn, beyond_15Re.
    region_name = 'post_midn'
    search_name = 'beyond_15Re'
    ; Good ones.
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2008-01-07/04:46','2008-01-07/05:04']), $
;        'probes', ['thc','thd','thb','tha'])
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2008-01-09/11:27','2008-01-09/11:39']), $
;        'probes', ['thc','tha','the','thd'])
;    ; bad ones.
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2007-11-22/16:46','2007-11-22/17:03']), $
;        'probes', ['thd','thc','the','thb'])

    ; post_midn, within_15Re.
    region_name = 'post_midn'
    search_name = 'within_15Re'
    ; Good ones.
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2014-08-28/10:11','2014-08-28/10:47']), $
;        'probes', ['thd','the','g15','tha','rbspb','g13'])
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2016-10-13/12:22','2016-10-13/12:45']), $
;        'probes', ['rbspa','rbspb','g15','thd','g14','g13'])
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2016-11-26/16:23','2016-11-26/16:57']), $
;        'probes', ['tha','thd','the','g15'])
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2016-12-11/09:46','2016-12-11/10:03']), $
;        'probes', ['tha','thd','the','g13'])


    ; pre_midn, beyond_15Re.
    region_name = 'pre_midn'
    search_name = 'beyond_15Re'
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2009-04-07/07:05','2009-04-07/07:08']), $
;        'probes', ['thd','tha','thb','thc'])
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2008-03-27/00:25','2008-03-27/01:02']), $
;        'probes', ['tha','thd','the','thc'])
;    event_list.add, dictionary($
;        'region', region_name, $
;        'search_name', search_name, $
;        'time_range', time_double(['2008-04-19/21:33','2008-04-19/21:42']), $
;        'probes', ['thc','thd','the','tha'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2009-03-30/01:48','2009-03-30/02:23']), $
        'probes', ['tha','thb','thc','thd','the'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2009-04-11/06:26','2009-04-11/06:44']), $
        'probes', ['thd','the','tha','thc'])

    ; pre_midn, within_15Re.
    region_name = 'pre_midn'
    search_name = 'within_15Re'
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2014-12-26/01:06','2014-12-26/01:19']), $
        'probes', ['g13','tha','the','g15'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2016-06-15/04:19','2016-06-15/04:39']), $
        'probes', ['g13','mms1','g14','g15','tha'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2017-03-28/03:00','2017-03-28/03:14']), $
        'probes', ['tha','the','thd','rbspa','g15'])
    event_list.add, dictionary($
        'region', region_name, $
        'search_name', search_name, $
        'time_range', time_double(['2017-03-31/02:47','2017-03-31/03:01']), $
        'probes', ['g13','tha','the','thd','g15'])



;---Find the wanted event.
    search_step = 'df'
    events = list()
    events.add, azim_df_search_post_midn_events(search_step=search_step), /extract
    events.add, azim_df_search_pre_midn_events(search_step=search_step), /extract
    good_events = list()
    foreach event, event_list do begin
        lprmsg, 'Finding event :'+strjoin(time_string(event.time_range),' to ')+' ...'
        foreach tmp, events do begin
            if event.region ne tmp.region then continue
            if event.search_name ne tmp.search_name then continue
            if product(event.time_range-tmp.time_range) gt 0 then continue

            tmp_time_range = event.time_range+[-1,1]*120
            df_list = tmp.df_list
            flags = bytarr(df_list.length)
            foreach df, df_list, df_id do begin
                flags[df_id] = product(df.obs_time-tmp_time_range) le 0
            endforeach
            index = where(flags eq 1, count)
            if count eq 0 then stop
            tmp.df_list = df_list[index]
            tmp.probes = event.probes
            good_events.add, tmp
            lprmsg, 'Found the event ...'
            break
        endforeach
    endforeach


    scale_width = 4.
    foreach event, good_events, event_id do begin
    ;---Load data.
        probes = event.probes
        data_time_range = event.time_range+[-1,1]*1800.
        foreach probe, probes do begin
            prefix = probe+'_'
            theta_var = prefix+'scaled_theta'
            if check_if_update(theta_var, data_time_range) then begin
                azim_df_read_data, 'theta', time_range=data_time_range, probe=probe, project=project
                azim_df_read_data, 'r_sm', time_range=data_time_range, probe=probe, project=project
                theta = get_var_data(prefix+'theta', times=times)
                r_sm = get_var_data(prefix+'r_sm', at=times)
                mlt = pseudo_mlt(r_sm)
                store_data, prefix+'pseudo_mlt', times, mlt
                store_data, theta_var, times, azim_df_scale_theta(theta, mlt, width=scale_width)
            endif
        endforeach


        df_list = list()
        df_list.add, event.df_list, /extract

        large_df_list = azim_df_filter_large_df(df_list)



    ;---Plot settings.
        xtickformat = ''
        xticklen_chsz = -0.15   ; in ychsz.
        yticklen_chsz = -0.30   ; in xchsz.
        full_ychsz = constant('full_ychsz')
        half_ychsz = full_ychsz*0.5
        lineskip = constant('lineskip')
        label_size = 0.7
        label_xshift = 10
        label_yshift = full_ychsz
        secofday = constant('secofday')
        bar_thick = keyword_set(test)? 0.5: 4

        rad = constant('rad')
        deg = constant('deg')

        mlt_range = [-1,1]*9
        rxy_range = [4.,30]


        sgopen, 0, xsize=6, ysize=4, /inch
        pos_list = list()
        pos_list.add, sgcalcpos(1, xchsz=xchsz, ychsz=ychsz, rmargin=10, lmargin=12)

        tpos = pos_list[0]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ;---The common x-axis setting.
        xrange = event_list[event_id].time_range+[-1,1]*1800.
        xstep = constant('secofhour')
        if total(xrange*[-1,1]) le xstep then xstep = 600.
        xticks = floor(total(xrange*[-1,1])/xstep)
        xtickv = smkarthm(xrange[0], xstep, xticks+1, 'x0')
        xminor = 6
        xtickn = strarr(xticks+1)
        for ii=0, xticks do begin
            the_time = xtickv[ii]
            xtickn[ii] = time_string(the_time,tformat='hh:mm')
            date = time_string(the_time,tformat='YYYY-MM-DD')
            if ii eq 0 then begin
                xtickn[ii] += '!C'+date
                continue
            endif
            if the_time mod secofday ne 0 then continue
            xtickn[ii] += '!C'+date
        endfor

    ;---Common z-range, colorbar, symbol.
        theta_var = 'theta'
        ztitle = 'Log!D10!N[Detrended tilt angle (deg)]'
        theta_range = [-1,1]*64
        zrange = theta_range
        ztickv = [-64,-16,-4,0,4,16,64]
        ztickn = string(ztickv,format='(I0)')
        ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
        zrange = alog10(abs(zrange))*[-1,1]
        zticks = n_elements(ztickv)-1
        spec_psym = 8
        spec_symsize = 0.3
        usersym, [1,1,-1,-1,1]*0.5, [-1,1,1,-1,-1]*1, /fill
        spec_ct = 70


    ;---Common x data.
        time_step = project.time_step
        xxs = make_bins(xrange, time_step)
        nxx = n_elements(xxs)

    ;---MLT-UT.
        pos_var = 'pseudo_mlt'
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        ytitle = 'MLT (hr)'
        yminor = 3
        yrange = []
        foreach probe, event.probes do yrange = [yrange,minmax(get_var_data(probe+'_pseudo_mlt'))]
        yrange = [floor(min(yrange)),ceil(max(yrange))]
        ytickv = make_bins(yrange,yminor, /inner)
        yticks = n_elements(ytickv)-1
        constants = 0


        ; Set up coord.
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, /nodata, /noerase

        ; Plot data.
        foreach probe, event.probes do begin
            prefix = probe+'_'
            yys = get_var_data(prefix+pos_var, at=xxs)
            zzs = get_var_data(prefix+theta_var, at=xxs)
            mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
            zzs = azim_df_scale_theta(zzs, mlt) ; apply MLT scaling.
            zzs = azim_df_normalize_theta(zzs, zrange=theta_range, ct=spec_ct, /reverse_ct)

            ; Remove data outside ROI.
            index = where_pro(mlt, '][', mlt_range, count=count)
            if count ne 0 then yys[index] = !values.f_nan
            rsm = get_var_data(prefix+'r_sm', at=xxs)
            rxy = snorm(rsm[*,0:1])
            index = where_pro(rxy, '][', rxy_range, count=count)
            if count ne 0 then yys[index] = !values.f_nan
            index = where(finite(zzs,/nan), count)
            if count ne 0 then yys[index] = !values.f_nan

            ; Plot data.
            for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=spec_psym, symsize=spec_symsize
        endforeach

        ; Add DF.
        foreach df, df_list do begin
            ty = df.obs_mlt
            tx = df.obs_time
            plots, tx,ty,/data, psym=1, symsize=0.5, color=sgcolor('gray')
        endforeach

        ; Add large DF.
        foreach df, large_df_list do begin
            ty = df.obs_mlt
            tx = df.obs_time
            plots, tx,ty,/data, psym=6, symsize=0.5, color=sgcolor('red')
        endforeach

        ; Add notations.
        foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

        ; Add axes.
        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase

        ; Add label.
        tx = tpos[0]-label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, 'c. UT-MLT'

        cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
        sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


        ;foreach df, df_list do azim_df_write_df, df

        stop
    endforeach





end
