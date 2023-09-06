;+
;-

function normalize_data, zzs, mlt, smooth_width=smooth_width, zrange=range, ct=ct, reverse_ct=reverse_ct

    zzs -= smooth(zzs,smooth_width, /edge_truncate, /nan)
    mean = mean(zzs,/nan)
    zzs = (zzs-mean)*sqrt(abs(mlt)+1)+mean

    linear_range = 1
    signs = sign(zzs)
    data = alog10(abs(zzs)>linear_range)*signs
    the_range = [-1,1]*alog10(max(abs(range)))
    index_color = bytscl(data, min=the_range[0], max=the_range[1])
    if keyword_set(reverse_ct) then index_color = 256-index_color

    if n_elements(ct) eq 0 then ct = 70
    loadct, ct, rgb_table=rgb0
    rgb0 = float(rgb0)
    true_colors = 256L*(256L*rgb0[index_color,2]+rgb0[index_color,1])+rgb0[index_color,0]

    return, true_colors
end

pro azim_df_gen_survey_plot, time_range, data_file=data_file, $
    plot_file=plot_file, reload_data=reload, $
    ; To put a thick bar for each of time_ranges in [n,2].
    highlight_time_ranges=highlight_time_ranges, $
    xsm_range=xsm_range

;---Prepare data for the event.
    load_data = 0
    del_data, '*'
    if keyword_set(reload) then file_delete, data_file, /allow_nonexistent
    if n_elements(data_file) eq 0 then load_data = 1 else if file_test(data_file) eq 0 then load_data = 1
    if load_data then azim_df_load_data_per_event, time_range=time_range, filename=data_file else tplot_restore, filename=data_file
    if file_test(data_file) eq 0 then return


;---Prepare to plot.
    if n_elements(plot_file) eq 0 then begin
        if keyword_set(test) then plot_file = 0 else plot_file = join_path([homedir(),'fig_azim_prop_storm_keox_ewo_'+time_string(time_range[0],tformat='YYYY_MMDD')+'.pdf'])
    endif
    ; Cut off data outside the range.
    if n_elements(mlt_range) ne 2 then mlt_range = [-1,1]*9
    if n_elements(xsm_range) ne 2 then xsm_range = [-20,4]
    if n_elements(dis_range) ne 2 then dis_range = minmax(abs(xsm_range))

;---Figure out the size of the plot.
    secofday = 86400d
    duration = total(time_range*[-1,1])
    nday = duration/secofday
    ysize = 6.
    aspect_ratio = 1.8
    xsize = ysize*aspect_ratio*nday ; normalize by duration.
    sgopen, plot_file, xsize=xsize, ysize=ysize

    npanel = 3
    poss = sgcalcpos(npanel, ypans=[1,2,2], tmargin=2, rmargin=8, lmargin=12, xchsz=xchsz, ychsz=ychsz)
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 2
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1
    time_step = 10.


;---The common xrange.
    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    xrange = time_range
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
    if tnames(vars[0]) eq '' then omni_read_index, time_range
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

;---Add highlighted time ranges.
    bar_color = sgcolor('blue')
    bar_thick = (size(plot_file,/type) eq 7)? 8: 2
    nhighlight_time_range = n_elements(highlight_time_ranges)/2
    if nhighlight_time_range gt 0 then begin
        for ii=0, nhighlight_time_range-1 do plots, highlight_time_ranges[ii,*], yrange[0]+[0,0], color=bar_color, thick=bar_thick
    endif


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
        zzs = normalize_data(zzs, mlt, smooth_width=smooth_width, ct=ct, reverse_ct=reverse_ct, zrange=tilt_range)

        ; Remove data out of range.
        get_data, prefix+'r_sm', times, r_sm
        xsm = r_sm[*,0]
        dis = snorm(r_sm)
        ; Exclude data outside the MLT range.
        index = where_pro(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside the x-range.
        index = where_pro(xsm, '][', xsm_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = where_pro(dis, '][', dis_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(cotran(r_sm, times, 'sm2gsm'))
        index = where(magn_flags eq 0, count)
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

;---Keoram.
    tpos = poss[*,2]
    xtickformat = ''
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
        zzs = normalize_data(zzs, mlt, smooth_width=smooth_width, ct=ct, reverse_ct=reverse_ct, zrange=tilt_range)

        ; Remove data out of range.
        get_data, prefix+'r_sm', times, r_sm
        xsm = r_sm[*,0]
        dis = snorm(r_sm)
        ; Exclude data outside the MLT range.
        index = where_pro(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside the x-range.
        index = where_pro(xsm, '][', xsm_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = where_pro(dis, '][', dis_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(cotran(r_sm, times, 'sm2gsm'))
        index = where(magn_flags eq 0, count)
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

    sgclose

end

; For all storms.
;project = azim_df_load_project()
;storm_times = project.storm_times
;storm_ids = project.storm_id
;foreach time_range, storm_times, ii do begin
;    data_file = join_path([project.data_dir,storm_ids[ii]+'_data.tplot'])
;    plot_file = join_path([project.plot_dir,'storm_survey',storm_ids[ii]+'_storm_survey.pdf'])
;    azim_df_gen_survey_plot, time_range, data_file=data_file, plot_file=plot_file
;endforeach

;; For the 10 storms.
;project = azim_df_load_project()
;event_infos = project.candidate
;foreach event_info, event_infos do begin
;    event_id = event_info.id
;    data_file = join_path([project.data_dir,event_info.data_file_suffix])
;    time_range = event_info.time_range
;    plot_file = join_path([project.plot_dir,'10_storm_survey',event_id+'_survey_plot.pdf'])
;    azim_df_gen_survey_plot, time_range, data_file=data_file, plot_file=plot_file, xsm_range=[-15,4]
;endforeach

; For Runov+2009.
time_range = time_double('2009-02-27')+[0,86400]
azim_df_gen_survey_plot, time_range, plot_file=join_path([homedir(),'runov_2009_paper.pdf']), xsm_range=[-25,0]


;time_range = storm_times[50]
;time_range = time_double('2009-02-07')+[0,86400]
;time_range = time_double('2009-02-13')+[0,86400]
;time_range = time_double('2009-02-04')+[0,86400]
;time_range = time_double('2009-02-27')+[0,86400]
;time_range = time_double('2016-10-28')+[0,86400]
;time_range = time_double('2007-12-17')+[0,86400]
;azim_df_gen_survey_plot, time_range

end
