;+
; Generate survey plot for event (not candidate). The difference is that here
; We use azim_df_load_basic_data to load more data, especially the triad_flag,
; which is plotted here.
;-

pro azim_df_event_gen_survey_plot, event_time_range, event_id=event_id, project=project, $
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
    ; Cut off data outside the range.
    if n_elements(mlt_range) ne 2 then mlt_range = [-1,1]*9
    if n_elements(xsm_range) ne 2 then xsm_range = [-20,4]
    if n_elements(dis_range) ne 2 then dis_range = minmax(abs(xsm_range))


;---Figure out the size of the plot.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz,ychsz=abs_ychsz
    sgclose, /wdelete
    ; Get the sizes of the panels.
    panel_ysize = 3   ; inch.
    ypans = [1,2,2,1]
    nypanel = n_elements(ypans)
    panel_ypad = 0.4
    panel_aspect_ratio = 2
    secofday = 86400d
    duration = total(event_time_range*[-1,1])
    nday = duration/secofday
    panel_xsize = panel_ysize*panel_aspect_ratio*nday ; normalize by duration.

    margins = [12.,5,8,2]
    fig_xsize = panel_xsize+total(margins[[0,2]])*abs_xchsz
    fig_ysize = panel_ysize+total(margins[[1,3]])*abs_ychsz+panel_ypad*(nypanel-1)
    if n_elements(plot_file) eq 0 then plot_file = join_path([project.plot_dir,'fig_'+event_id+'_survey_plot.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize

    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    poss = sgcalcpos(nypanel, position=pos, ypans=ypans, ypad=panel_ypad, xchsz=xchsz, ychsz=ychsz)

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


;---Dst and AE.
    tpos = poss[*,0]
    xtickformat='(A1)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    vars = ['dst','ae']
    foreach var, vars do if check_if_update(var,event_time_range) then omni_read_index, event_time_range
    options, vars, 'yticklen', yticklen
    options, vars, 'ystyle', 1
    options, vars, 'yminor', 5
    ystep = [50,500]
    constant = [0,500]
    ytitles = ['Dst','AE']+' (nT)'
    foreach var, vars, ii do begin
        get_data, var, times, data
        ytickv = make_bins(data, ystep[ii])
        yticks = n_elements(ytickv)-1
        options, var, 'ytitle', ytitles[ii]
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'yrange', minmax(ytickv)
        options, var, 'constant', constant[ii]
    endforeach

    get_data, 'dst', xxs, yys, limits=lim
    yrange = lim.yrange
    ytickv = lim.ytickv
    yticks = lim.yticks
    yminor = lim.yminor
    ytitle = lim.ytitle
    constant = lim.constant
    dst_color = sgcolor('black')
    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase
    oplot, xxs, yys, color=dst_color
    oplot, xrange, constant+[0,0], linestyle=constant_linestyle
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=9, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    get_data, 'ae', xxs, yys, limits=lim
    yrange = lim.yrange
    ytickv = lim.ytickv
    yticks = lim.yticks
    yminor = lim.yminor
    ytitle = lim.ytitle
    constant = lim.constant
    ae_color = sgcolor('red')
    axis, yaxis=1, /save, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        color=ae_color
    oplot, xxs, yys, color=ae_color
    oplot, xrange, constant+[0,0], linestyle=constant_linestyle, color=color

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'a. Dst/'
    ;ty = tpos[3]-label_yshift*ychsz-full_ychsz*ychsz
    xyouts, tx,ty,/normal, '           AE', color=ae_color


;---Common z-range, colorbar, symbol.
    data_var_suffix = 'db_tilt'
    ztitle = 'Log!D10!N[Tilt angle - T89 (deg)]'
    tilt_range = [-1,1]*64
    zrange = tilt_range
    ztickv = [-64,-16,-4,0,4,16,64]
    ztickn = string(ztickv,format='(I0)')
    ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
    zrange = alog10(abs(zrange))*[-1,1]
    zticks = n_elements(ztickv)-1

    ; Symbol.
    psym = 8
    symsize = 0.5
    usersym, [1,1,-1,-1,1]*0.1, [-1,1,1,-1,-1]*1
    ct = 70
    reverse_ct = 1

    ; About transform data.
    smooth_window = 80*60   ; sec.
    prefixes = tnames('*_mlt')
    foreach prefix, prefixes, ii do prefixes[ii] = get_prefix(prefixes[ii])


;---EWOgram.
    tpos = poss[*,1]
    xtickformat = '(A1)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    ytitle = 'MLT (hr)'
    yrange = mlt_range
    yminor = 3
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase

    pos_var_suffix = 'mlt'
    available_prefixes = prefixes
    flags = list()
    foreach prefix, prefixes do begin
        pos_var = prefix+pos_var_suffix
        data_var = prefix+data_var_suffix
        get_data, data_var, xxs, zzs
        get_data, pos_var, xxs, yys
        if n_elements(zzs) ne n_elements(xxs) then begin
            errmsg = handle_error(data_var+': inconsistent data and time ...')
            stop
        endif
        mlt = get_var_data(prefix+'mlt')
        smooth_width = smooth_window/time_step
        zzs = azim_df_normalize_theta(zzs, mlt, smooth_width=smooth_width, ct=ct, reverse_ct=reverse_ct, zrange=tilt_range)

        ; Remove data out of range.
        get_data, prefix+'r_sm', times, r_sm
        xsm = r_sm[*,0]
        dis = snorm(r_sm)
        ; Exclude data outside the MLT range.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = lazy_where(dis, '][', dis_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(cotran(r_sm, times, 'sm2gsm'))
        index = where(magn_flags eq 0, count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside time range.
        index = lazy_where(xxs, '][', event_time_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside yrange.
        index = lazy_where(yys, '][', yrange, count=count)
        if count ne 0 then xxs[index] = !values.f_nan

        index = where(finite(xxs), count)
        the_flag = (count gt 0)
        flags.add, the_flag
        if count eq 0 then continue

        nxx = count
        xxs = xxs[index]
        yys = yys[index]
        zzs = zzs[index]
        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize

    endforeach
    oplot, xrange, constant+[0,0], linestyle=constant_linestyle

    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'b. UT-'+strmid(ytitle,0,strpos(ytitle,' '))

    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks

    flags = flags.toarray()
    index = where(flags eq 1, count)
    available_prefixes = available_prefixes[index]

;---Keogram.
    tpos = poss[*,2]
    xtickformat = '(A1)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    ytitle = 'X!USM!N (Re)'
    yrange = xsm_range
    yminor = 5
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase

    pos_var_suffix = 'r_sm'
    foreach prefix, prefixes do begin
        pos_var = prefix+pos_var_suffix
        data_var = prefix+data_var_suffix
        get_data, data_var, xxs, zzs
        get_data, pos_var, xxs, yys & yys = yys[*,0]
        if n_elements(zzs) ne n_elements(xxs) then begin
            errmsg = handle_error(data_var+': inconsistent data and time ...')
            stop
        endif
        mlt = get_var_data(prefix+'mlt')
        smooth_width = smooth_window/time_step
        zzs = azim_df_normalize_theta(zzs, mlt, smooth_width=smooth_width, ct=ct, reverse_ct=reverse_ct, zrange=tilt_range)

        ; Remove data out of range.
        get_data, prefix+'r_sm', times, r_sm
        xsm = r_sm[*,0]
        dis = snorm(r_sm)
        ; Exclude data outside the MLT range.
        index = lazy_where(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = lazy_where(dis, '][', dis_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(cotran(r_sm, times, 'sm2gsm'))
        index = where(magn_flags eq 0, count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside time range.
        index = lazy_where(xxs, '][', event_time_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside yrange.
        index = lazy_where(yys, '][', yrange, count=count)
        if count ne 0 then xxs[index] = !values.f_nan

        index = where(finite(xxs), count)
        if count eq 0 then continue

        nxx = count
        xxs = xxs[index]
        yys = yys[index]
        zzs = zzs[index]
        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach

    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'c. UT-'+strmid(ytitle,0,strpos(ytitle,' '))

    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;---Print available_probes.
    msg = ['Available spacecraft: ']
    prefixes = sort_uniq(available_prefixes)
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


;---Available triad.
    tpos = poss[*,3]
    xtickformat = ''
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])

    labels = ['post','pre']
    vars = labels+'_midn_triad_flag'
    colors = sgcolor(['blue','red'])

    ytitle = '(#)'
    yrange = []
    foreach var, vars do yrange = [yrange,minmax(get_var_data(var,in=event_time_range))]
    yrange = minmax(yrange)
    yminor = 5
    ytickv = make_bins(yrange, yminor)
    yrange = minmax(ytickv)
    yticks = n_elements(ytickv)
    if yticks gt 4 then begin
        yminor *= 2
        ytickv = make_bins(yrange, yminor)
        yrange = minmax(ytickv)
        yticks = n_elements(ytickv)
    endif
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase

    dy = (tpos[3]-tpos[1])/(n_elements(vars)+1)
    foreach var, vars, ii do begin
        get_data, var, xxs, yys
        oplot, xxs,yys, color=colors[ii]
        tx = tpos[2]+xchsz*0.5
        ty = tpos[3]-dy*(ii+1)-ychsz*0.35
        xyouts, tx,ty,/normal, strupcase(labels[ii]), color=colors[ii]
    endforeach


    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'd. Triad'
    tx = tpos[0]+xchsz*0.5*label_size
    ty = tpos[3]-ychsz*1*label_size
    small_angle_limit = 15.
    angle_range = [small_angle_limit,180-small_angle_limit]
    xyouts, tx,ty,/normal, 'Good triad: angles in ['+sgnum2str(angle_range[0])+','+sgnum2str(angle_range[1])+']', charsize=label_size

;---Add highlighted time ranges.
    tpos = poss[*,1]
    bar_color = sgcolor('black')
    bar_thick = (size(plot_file,/type) eq 7)? 8: 2
    nhighlight_time_range = n_elements(highlight_time_ranges)/2
    if nhighlight_time_range gt 0 then begin
        for ii=0, nhighlight_time_range-1 do plots, highlight_time_ranges[ii,*], yrange[0]+[0,0], color=bar_color, thick=bar_thick
    endif

    min_num_triad = 5.
    min_triad_duration = 60*60.     ; sec.
    triad_section_pad = 0
    colors = sgcolor(['blue','red'])
    foreach region, ['post','pre'], ii do begin
        the_var = region+'_midn_triad_flag'
        get_data, the_var, uts, flags
        index = where(flags ge min_num_triad, count)
        if count eq 0 then continue

        bar_color = colors[ii]
        bar_shift = ychsz*0.1

        plot, xrange, [0,1], xstyle=5, ystyle=5, /nodata, /noerase, position=tpos
        bar_ypos = (region eq 'post')? tpos[3]-bar_shift: tpos[1]+bar_shift
        tmp = convert_coord([0,bar_ypos], /normal, /to_data)
        ty = tmp[1]

        msg_ypos = (region eq 'post')? tpos[3]-ychsz*full_ychsz*label_size-bar_shift*2: tpos[1]+bar_shift*4
        msg_xpos = tpos[0]+xchsz*0.5
        msg = region+'-midn'
        msg_shift = 0
        xyouts, msg_xpos,msg_ypos-msg_shift,/normal, alignment=0, msg, charsize=label_size, color=bar_color

        uts = uts[index]
        time_ranges = time_to_range(uts, time_step=time_step, pad_times=triad_section_pad)
        time_ranges >= event_time_range[0]
        time_ranges <= event_time_range[1]
        ntime_range = n_elements(time_ranges)/2

        for jj=0, ntime_range-1 do begin
            time_range = time_ranges[jj,*]
            section_duration = total(time_range*[-1,1])
            if section_duration lt min_triad_duration then continue
            plots, time_range, ty+[0,0], color=bar_color, thick=bar_thick
        endfor
    endforeach

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
    plot_file = join_path([project.plot_dir,'good_storm_survey_with_triad_flag',event_id+'_survey_plot.pdf'])
    azim_df_event_gen_survey_plot, time_range, data_file=data_file, plot_file=plot_file, xsm_range=[-15,4]
endforeach


end
