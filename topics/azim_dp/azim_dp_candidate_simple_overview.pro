;+
; Plot SC position, AE/Dst, and theta.
; Do not need to be an event, i.e., have arrows.
;
; Adopted from azim_df_event_survey_plot2.
;-

pro azim_dp_candidate_simple_overview, time_range, probes=all_probes, $
    mlt_range=mlt_range, rxy_range=rxy_range, test=test, $
    panel_x_ratio=panel_x_ratio, no_dp_plotted=no_dp_plotted

;test = 0


;---Check input.
    if n_elements(time_range) ne 2 then message, 'No time range ...'
    if n_elements(all_probes) eq 0 then message, 'No probes ...'
    if n_elements(mlt_range) ne 2 then mlt_range = [-1,1]*9
    if n_elements(rxy_range) ne 2 then rxy_range = [5.,20]

;---Settings.
    secofhour = 3600.
    probe_color_ct = 40
    color_bottom = 50
    color_top = 250
    time_type = '_ramp'

;---Load data.
    plot_time_range = time_range
    omni_read_index, plot_time_range

    ; Search DP in ROI.
    roi_dp_list = azim_dp_search_dp_in_roi($
        azim_dp_search_roi(time_range, probes=all_probes, mlt_range=mlt_range, rxy_range=rxy_range))


    ; Get the probes have DP.
    dp_list = list()
    foreach event, roi_dp_list do dp_list.add, event.dp_list, /extract
    dp_probes = list()
    foreach dp, dp_list do begin
        if dp_probes.where(dp.probe) ne !null then continue
        dp_probes.add, dp.probe
    endforeach

    ; Load data for certain probes.
    plot_probes = list()
    foreach probe, all_probes do begin
        prefix = probe+'_'
        vars = prefix+['theta','mlt','scaled_theta','r_sm']
        del_data, vars

        azim_dp_read_theta, plot_time_range, probe=probe
        azim_dp_read_mlt, plot_time_range, probe=probe
        azim_dp_read_orbit, plot_time_range, probe=probe

        get_data, prefix+'theta', times, theta
        if n_elements(times) le 2 then continue
        get_data, prefix+'mlt', times
        if n_elements(times) le 2 then continue
        get_data, prefix+'r_sm', times
        if n_elements(times) le 2 then continue

        get_data, prefix+'theta', times
        interp_time, prefix+'mlt', times
        interp_time, prefix+'r_sm', times

        mlt = get_var_data(prefix+'mlt')
        store_data, prefix+'scaled_theta', times, azim_dp_scale_theta(theta, mlt)
        plot_probes.add, probe
    endforeach


    ; Collect info for DPs.
    ndp = dp_list.length
    ndim = 3
    probes = dp_probes.toarray()
    if ndp eq 0 then begin
        ramp_times = []
        ramp_mlts = []
        ramp_r_sms = fltarr(1,ndim)+!values.f_nan
        ramp_rxys = []
    endif else begin
        ramp_times = dblarr(ndp)
        ramp_mlts = fltarr(ndp)
        ramp_r_sms = fltarr(ndp,ndim)
        foreach dp, dp_list, ii do begin
            ramp_times[ii] = dp.time
            ramp_mlts[ii] = dp.mlt
            ramp_r_sms[ii] = dp.r_sm
        endforeach
        ramp_rxys = snorm(ramp_r_sms[*,0:1])
    endelse



    ; Figure out the colors.
    nplot_probe = n_elements(plot_probes)
    center_mlts = fltarr(nplot_probe)
    center_time = mean(plot_time_range)
    foreach probe, plot_probes, probe_id do begin
        prefix = probe+'_'
        center_mlts[probe_id] = get_var_data(prefix+'mlt', at=center_time)
    endforeach
    index = sort(center_mlts)
    sorted_probes = plot_probes[index]
    nprobe = n_elements(sorted_probes)
    probe_colors = smkarthm(color_bottom,color_top,nprobe, 'n')
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)


;---Plot settings.
    pos_xy_range = minmax(make_bins(minmax([-1,1,ramp_rxys]),5))
    pos_xrange = minmax(make_bins(minmax([-1,1,ramp_r_sms[*,0]]),5))
    pos_yrange = minmax(make_bins(minmax([-1,1,ramp_r_sms[*,1]]),5))
    pos_xrange = [-30,10]
    pos_yrange = [-1,1]*20
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 2
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=xchsz, ychsz=ychsz
    sgclose, /wdelete
    abs_xchsz = xchsz
    abs_ychsz = ychsz
    ; Get the sizes of the left panels.
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1,1]
    nypanel = n_elements(ypans)
    panel_ypads = [0.4,0.4,3]
    if n_elements(panel_x_ratio) eq 0 then panel_x_ratio = 2.0     ; inch per hour.
    panel_xsize = abs(total(plot_time_range*[-1,1]))/3600*panel_x_ratio
    ; Get the size of the right panel.
    right_pan_ysize = total(panel_ypads)*abs_ychsz+total(ypans)*panel_ysize
    right_pan_xsize = right_pan_ysize


    margins = [12.,5,4,2]
    panel_xpad = 14

    fig_ysize = right_pan_ysize+total(margins[[1,3]])*abs_ychsz
    fig_xsize = panel_xsize+panel_xpad*abs_xchsz+right_pan_xsize+total(margins[[0,2]])*abs_xchsz
    event_id = time_string(mean(time_range),tformat='YYYY_MMDD')
    file = join_path([homedir(),'azim_dp_candidate','fig_azim_dp_candidate_'+event_id+'_simple_overview.pdf'])
    if keyword_set(test) then file = test
    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch

    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    lower_pos = pos
    poss = sgcalcpos(1,2, position=lower_pos, xpans=[panel_xsize,right_pan_xsize], xpad=panel_xpad)
    left_pos = poss[*,0]
    left_poss = sgcalcpos(nypanel, position=left_pos, ypans=ypans, ypad=panel_ypads, xchsz=xchsz, ychsz=ychsz)
    right_pos = poss[*,1]
    tx = total(right_pos[[0,2]])*0.5
    ty = total(right_pos[[1,3]])*0.5
    dx = total(right_pos[[0,2]]*[-1,1])*0.4
    dy = total(right_pos[[1,3]]*[-1,1])*0.4
    right_pos = [tx-dx,ty-dy,tx+dx,ty+dy]


;---Common x-axis.
    xrange = plot_time_range
    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    xstep = 30*60   ; min.
    xminor = 6
    xtickv = make_bins(xrange, xstep, /inner) ; make time line up.
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    secofday = 86400d
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



;---Panel a. Dst/AE.
    tpos = left_poss[*,0]
    xtickformat = '(A1)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    the_var = 'dst'
    ystep = 50.
    constants = [-50,0]
    yys = get_var_data(the_var, in=plot_time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    ytitle = '(nT)'

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Add data.
    oplot, xxs, yys
    foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=9, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    the_var = 'ae'
    ystep = 500.
    constants = [0,500]
    yys = get_var_data(the_var, in=plot_time_range, times=xxs)
    index = where(yys ge 0, count)
    if count eq 0 then yys = !null
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    ytitle = '(nT)'

    ; Set up coord.
    ae_color = sgcolor('red')
    axis, yaxis=1, /save, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        color=ae_color

    ; Add data.
    oplot, xxs, yys, color=ae_color
    foreach constant, constants do oplot, xrange, constant+[0,0], color=ae_color, linestyle=1

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'a. Dst/'
    xyouts, tx,ty,/normal, '           AE', color=ae_color


;---Color bar.
    tpos = left_poss[*,1]
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    tpos = left_poss[*,2]
    cbpos[1] = tpos[1]


;---Panel b. UT-MLT short.
    tpos = left_poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    ytitle = 'MLT (hr)'
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0

    ztitle = 'Scaled detr.tilt '+tex2str('theta')+' (deg)'
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

    time_step = 10.
    times = make_bins(plot_time_range+[0,-1]*time_step, time_step)
    ntime = n_elements(times)
    thetas = fltarr(ntime,nprobe)
    mlts = fltarr(ntime,nprobe)
    rxys = fltarr(ntime,nprobe)
    rsms = fltarr(ntime,3,nprobe)
    foreach probe, plot_probes, probe_id do begin
        prefix = probe+'_'
        the_var = prefix+'scaled_theta'
        if tnames(the_var) eq '' then continue
        thetas[*,probe_id] = get_var_data(the_var, at=times)
        mlts[*,probe_id] = get_var_data(prefix+'mlt', at=times)
        rsms[*,*,probe_id] = get_var_data(prefix+'r_sm', at=times)
        rxys[*,probe_id] = snorm(reform(rsms[*,0:1,probe_id]))
    endforeach

    ystep = 3
    yrange = minmax(make_bins(minmax(mlts),ystep))
    yminor = 3
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1


    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase

    foreach time, times, time_id do begin
        zzs = thetas[time_id,*]
        yys = mlts[time_id,*]
        rxy = rxys[time_id,*]
        rsm = rsms[time_id,*,*]

        index = where(finite(zzs),count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = azim_df_normalize_theta(zzs[index], zrange=tilt_range, ct=ct, /reverse_ct)
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Filter spatially.
        index = where_pro(yys, '[]', yrange, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = where_pro(rxy, '[]', rxy_range, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = where(check_if_in_magn(rsm, dynamic_pressure=pdyn) eq 1, count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Sort by tilt.
        index = reverse(sort(zzs))
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        for ii=0, count-1 do plots, time, yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach

    ; Add DP.
    square_symbol = 6
    square_size = label_size*0.8
    if ~keyword_set(no_dp_plotted) then begin
        foreach dp, dp_list, ii do begin
            ;the_color = probe_colors[where(plot_probes eq dp.probe)]
            the_color = sgcolor('black')
            plots, dp.time,dp.mlt,/data, color=the_color, psym=square_symbol, symsize=square_size
        endforeach
    endif


    ; Draw box.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


    ; Color bar.
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks

    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'b. UT-'+strmid(ytitle,0,strpos(ytitle,' '))




;---Panel b. UT-X short.
    tpos = left_poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = ''
    ytitle = 'Rx (Re)'
    yminor = 5
    yrange = minmax(make_bins(rsms[*,0,*],yminor))
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0


    ; Symbol.
    psym = 8
    symsize = 0.5
    usersym, [1,1,-1,-1,1]*0.1, [-1,1,1,-1,-1]*1
    ct = 70
    reverse_ct = 1

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase


    foreach time, times, time_id do begin
        zzs = thetas[time_id,*]
        yys = rsms[time_id,0,*]
        rxy = rxys[time_id,*]
        rsm = rsms[time_id,*,*]

        index = where(finite(zzs),count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = azim_df_normalize_theta(zzs[index], zrange=tilt_range, ct=ct, /reverse_ct)
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Filter spatially.
        index = where_pro(yys, '[]', yrange, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = where_pro(rxy, '[]', rxy_range, count=count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        index = where(check_if_in_magn(rsm, dynamic_pressure=pdyn) eq 1, count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        ; Sort by tilt.
        index = reverse(sort(zzs))
        yys = yys[index]
        zzs = zzs[index]
        rxy = rxy[index]
        rsm = rsm[index,*]

        for ii=0, count-1 do plots, time, yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach

    square_symbol = 6
    square_size = label_size*0.8
    if ~keyword_set(no_dp_plotted) then begin
        foreach dp, dp_list, ii do begin
            ;the_color = probe_colors[where(plot_probes eq dp.probe)]
            the_color = sgcolor('black')
            plots, dp.time,dp.r_sm[0],/data, color=the_color, psym=square_symbol, symsize=square_size
        endforeach
    endif


    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


    tx = label_xshift*xchsz
    ty = tpos[3]-label_yshift*ychsz
    xyouts, tx,ty,/normal, 'c. UT-'+strmid(ytitle,0,strpos(ytitle,' '))







;---Right panel.
    tpos = right_pos
    panel_xrange = reverse(pos_xrange)
    panel_yrange = reverse(pos_yrange)

    hsize = keyword_set(test)? 8: 160
    xtitle = 'SM X (Re)'
    ytitle = 'SM Y (Re)'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    plot, panel_xrange, panel_yrange, xstyle=1, ystyle=1, position=tpos, $
        xrange=panel_xrange, yrange=panel_yrange, /nodata, /noerase, /isotropic, $
        xtitle=xtitle, ytitle=ytitle, xticklen=xticklen, yticklen=yticklen


    ; Add figure label.
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1.2
    label = 'd. SM x-y plane'
    xyouts, tx,ty,/normal,alignment=0, label

    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    plots, circle_x, circle_y
    foreach r, make_bins(minmax(abs(pos_xrange)),10,/inner) do oplot, circle_x*r, circle_y*r, linestyle=1


    ; Add SC as vertex.
    foreach probe, plot_probes do begin
        prefix = probe+'_'
        the_color = probe_colors[where(plot_probes eq probe)]
        the_color = the_color[0]

        rsm = get_var_data(prefix+'r_sm', in=time_range)
        xxs = rsm[*,0]
        yys = rsm[*,1]
        plots, xxs, yys, /data, color=the_color

        rsm = get_var_data(prefix+'r_sm', at=center_time)
        plots, rsm[0],rsm[1],/data, psym=square_symbol, symsize=square_size, color=the_color
;        short_name = probe
;        tmp = convert_coord(rsm[0:1], /data, /to_normal)
;        tx = tmp[0]+xchsz*0.5
;        ty = tmp[1]
;        xyouts, tx,ty,/normal, strupcase(short_name), color=the_color
    endforeach

    ty = tpos[1]+ychsz*0.2
    foreach probe, plot_probes, probe_id do begin
        the_color = probe_colors[where(plot_probes eq probe)]
        tx = tpos[0]+xchsz*(1+probe_id*4.5)
        probe_info = resolve_probe(probe)
        msg = strupcase(probe_info.short_name)
        xyouts, tx,ty,/normal, msg, color=the_color
    endforeach

    if keyword_set(test) then stop
    sgclose


end

time_range = time_double(['2017-03-21/07:00','2017-03-21/09:00'])
probes = ['rbsp'+letters('b'),'th'+['a','d','e'],'g'+['13','14','15','16','17'],'mms1']
probes = ['th'+['a','d','e'],'g'+['13','15']]
mlt_range = [-9,6]
rxy_range = [5.,20]

; E. Donovan's event.
time_range = time_double(['2008-03-05/06:00','2008-03-05/06:15'])
probes = ['th'+letters('e'),'g'+['11','12']]
mlt_range = [-6,6]
rxy_range = [3.,20]

; Ohtani's event.
time_range = time_double(['2011-03-09/05:00','2011-03-09/08:00'])
probes = ['th'+['a','d','e'],'g'+['13','15']]
mlt_range = [-6,6]
rxy_range = [5.,20]

; Ohtani+2017, Fig 11.
time_range = time_double(['2015-01-15/06:00','2015-01-15/08:00'])
probes = ['th'+['a','d','e'],'g'+['13','14','15','16','17'],'rbsp'+letters('b')]
mlt_range = [-6,6]
rxy_range = [5.,20]

; Ohtani+2017, Fig 7.
time_range = time_double(['2015-01-25/06:00','2015-01-25/08:00'])
probes = ['th'+['a','d','e'],'g'+['13','14','15','16','17'],'rbsp'+letters('b')]
mlt_range = [-6,6]
rxy_range = [5.,20]

; 2014-08-28 first substorm.
time_range = time_double(['2014-08-28/04:30','2014-08-28/06:00'])
probes = ['th'+['a','d'],'g'+['13','15'],'rbsp'+['a']]
mlt_range = [-9,9]
rxy_range = [2.,15]


; 2013-06-07 overview.
time_range = time_double(['2013-06-07/07:00','2013-06-07/10:00'])
probes = ['th'+['d','e'],'g'+['13','15']]
mlt_range = [-9,9]
rxy_range = [3.,20]



; Runov event.
time_range = time_double(['2009-02-27/07:30','2009-02-27/08:30'])
probes = ['th'+['a','b','c','d','e'],'g'+['10','11','12']]
mlt_range = [-6,6]
rxy_range = [3.,30]


; THEMIS event.
time_range = time_double(['2008-02-29/07:50','2008-02-29/09:00'])
probes = ['th'+['a','c','d','e'],'g'+['11','12']]
mlt_range = [-6,6]
rxy_range = [3.,30]

time_range = time_double(['2008-03-20/11:30','2008-03-20/13:00'])
probes = ['th'+['a','b','c','d','e'],'g'+['11','12']]
mlt_range = [-6,6]
rxy_range = [3.,30]


; Rumi Nakamura's event.
time_range = time_double(['2008-02-16/02:00','2008-02-16/03:30'])
probes = ['th'+['a','b','c','d','e'],'g'+['10','11','12']]
mlt_range = [-9,0]
rxy_range = [3.,30]


; Check Goes for my event.
time_range = time_double(['2008-02-29/08:00','2008-02-29/10:00'])
probes = ['th'+['a','c','d','e'],'g'+['10','11','12']]
mlt_range = [-6,0]
rxy_range = [3.,30]


azim_dp_candidate_simple_overview, time_range, probes=probes, $
    mlt_range=mlt_range, rxy_range=rxy_range, panel_x_ratio=2, no_dp_plotted=1, test=1
end


; Done.
;; 2013-06-07 overview.
;
;time_range = time_double(['2013-06-07/02:00','2013-06-07/04:00'])
;probes = ['th'+['a','d'],'g'+['13','15'],'rbsp'+['a','b']]
;mlt_range = [-9,9]
;rxy_range = [3.,20]
;
;azim_dp_candidate_simple_overview, time_range, probes=probes, $
;    mlt_range=mlt_range, rxy_range=rxy_range, panel_x_ratio=1, no_dp_plotted=1, test=0
;
;time_range = time_double(['2013-06-07/04:00','2013-06-07/07:00'])
;probes = ['th'+['a','d','e'],'g'+['13','15'],'rbsp'+['a','b']]
;mlt_range = [-6,9]
;rxy_range = [3.,20]
;
;azim_dp_candidate_simple_overview, time_range, probes=probes, $
;    mlt_range=mlt_range, rxy_range=rxy_range, panel_x_ratio=1, no_dp_plotted=1, test=0
;end
