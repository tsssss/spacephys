;+
; Generate survey plot for event (not candidate).
; Plot theta of all available_probes.
;-

pro azim_df_event_gen_theta_survey, event_time_range, event_id=event_id, project=project, $
    data_file=data_file, plot_file=plot_file, $
    ; To put a thick bar for each of time_ranges in [n,2].
    highlight_time_ranges=highlight_time_ranges, xsm_range=xsm_range

    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(event_time_range) eq 2 then begin
        event_id = time_string(mean(event_time_range),tformat='YYYY_MMDD_hh')
    endif else if n_elements(event_id) eq 0 then begin
        errmsg = handle_error('No input time_range or event_id ...')
        return
    endif

test = 0

;---Load basic data and event_info.
    azim_df_load_basic_data, event_time_range, event_id=event_id, project=project, reset=reset
    if n_elements(event_time_range) ne 2 then begin
        vars = tnames('*_r_sm')
        get_data, vars[0], times
        event_time_range = minmax(times)
    endif

;---Prepare to plot.
    if n_elements(plot_file) eq 0 then begin
        if keyword_set(test) then plot_file = 0 else plot_file = join_path([project.plot_dir,'fig_'+event_id+'_survey_plot.pdf'])
    endif
    ; Cut off data outside the range.
    if n_elements(mlt_range) ne 2 then mlt_range = [-1,1]*9
    if n_elements(xsm_range) ne 2 then xsm_range = [-20,4]
    if n_elements(dis_range) ne 2 then dis_range = minmax(abs(xsm_range))

;---Theta.
    constant = 0
    flags = list()
    available_prefixes = tnames('*_theta')
    foreach prefix, available_prefixes, ii do available_prefixes[ii] = strmid(prefix,0,strpos(prefix,'_')+1)
    foreach prefix, available_prefixes, ii do begin
        the_var = prefix+'theta'
        get_data, the_var, xxs, zzs, limits=lim
        if n_elements(zzs) ne n_elements(xxs) then begin
            errmsg = handle_error(data_var+': inconsistent data and time ...')
            stop
        endif

        ; Remove data out of range.
        get_data, prefix+'r_sm', times, r_sm
        mlt = get_var_data(prefix+'mlt')
        xsm = r_sm[*,0]
        dis = snorm(r_sm)
        ; Exclude data outside the MLT range.
        index = where_pro(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = where_pro(dis, '][', dis_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(cotran(r_sm, times, 'sm2gsm'))
        index = where(magn_flags eq 0, count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside time range.
        index = where_pro(xxs, '][', event_time_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan

        index = where(finite(xxs), count)
        the_flag = (count gt 0)
        flags.add, the_flag
        if count eq 0 then continue

        nxx = count
        xxs = xxs[index]
        zzs = zzs[index]

        store_data, the_var+'_plot', xxs, zzs, limits=lim
    endforeach

    flags = flags.toarray()
    index = where(flags eq 1, count)
    available_prefixes = available_prefixes[index]


;---Figure out the size of the plot.
    secofday = 86400d
    duration = total(event_time_range*[-1,1])
    nday = duration/secofday
    ysize = 6.
    aspect_ratio = 1.8
    xsize = ysize*aspect_ratio*nday ; normalize by duration.
    npanel = n_elements(available_prefixes)
    ysize = ysize*npanel/3
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=xsize, ysize=ysize, /inch

    margins = [12,5,8,2]
    poss = sgcalcpos(npanel, margins=margins, xchsz=xchsz, ychsz=ychsz)
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 2
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1
    time_step = 10.


;---The common xrange.
    xrange = event_time_range
    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    xminor = 8  ; hour.
    xtickv = make_bins(xrange, 3600*xminor, /inner) ; make time line up at hours.
    xticks = n_elements(xtickv)-1
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


;---Theta.
    prefixes = sort_uniq(available_prefixes)
    ytitle = tex2str('theta')+' (deg)'
    fig_labels = letters(npanel)
    foreach prefix, prefixes, ii do begin
        mission_probe = strmid(prefix,0,strpos(prefix,'_'))
        mission_info = resolve_probe(mission_probe)
        short_name = mission_info.short_name

        tpos = poss[*,ii]
        tx = label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, fig_labels[ii]+'. Tilt '+strupcase(short_name)

        mission_probe = strmid(prefix,0,strpos(prefix,'_'))
        add_setting, prefix+'theta_plot', {$
            ytitle: ytitle, $
            labels: strupcase(short_name)}
    endforeach


;---Print available_probes.
    msg = ['Available spacecraft: ']
    foreach prefix, prefixes do begin
        mission_probe = strmid(prefix,0,strpos(prefix,'_'))
        mission_info = resolve_probe(mission_probe)
        msg = [msg,strupcase(mission_info.short_name)]
    endforeach
    if n_elements(msg) eq 1 then msg = [msg,'none']
    msg = msg[0]+strjoin(msg[1:*],',')
    tpos = poss[*,0]
    tx = tpos[0]
    ty = tpos[3]+ychsz*half_ychsz
    xyouts, tx,ty,/normal, msg


;---Plot all data.
    tplot, prefixes+'theta_plot', position=poss, trange=event_time_range, /noerase

    if keyword_set(test) then stop
    sgclose

end


; For the 10 storms.
project = azim_df_load_project()
event_infos = project.candidate
foreach event_info, event_infos do begin
    event_id = event_info.id
    data_file = join_path([project.data_dir,event_info.data_file_suffix])
    time_range = event_info.time_range
    plot_file = join_path([project.plot_dir,'10_storm_survey',event_id+'_survey_plot3.pdf'])
    azim_df_event_gen_theta_survey, time_range, data_file=data_file, plot_file=plot_file, xsm_range=[-15,4]
endforeach


end
