;+
; Plot AE/SME, ewograms of upward and downward currents, and DP. 
; newer than figewogram_of_dp_and_up_down_current. Here it auto load all spacecraft.
;-

pro fig_ewogram_of_dp_and_up_down_current2, time_range, probes=probes, $
    filename=plot_file, test=test, $
    mlt_range=mlt_range, rxy_range=rxy_range, mlat_range=mlat_range, pdyn=pydn, $
    fix_x_ratio=fix_x_ratio, panel_x_ratio=panel_x_ratio, $
    avail_probes=avail_probes   ; output

    if n_elements(time_range) ne 2 then message, 'Invalid time_range ...'
    duration = total(time_range*[-1,1])
    if duration le 600 then message, 'time_range too short ...'

    if n_elements(mlt_range) ne 2 then mlt_range = [-9,9]
    if n_elements(rxy_range) ne 2 then rxy_range = [4.,30]
    if n_elements(mlat_range) ne 2 then mlat_range = [60,80]

    nprobe = n_elements(probes)
    if nprobe eq 0 then probes = ['th'+letters('e'),'rbsp'+letters('b'),'g'+string(make_bins([10,15],1),format='(I0)')]
    if n_elements(plot_file) eq 0 then plot_file = join_path([homedir(),$
        'fig_ewogram_of_dp_and_up_down_current_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hh'),'_')+'_v01.pdf'])
    if keyword_set(test) then plot_file = 0



;---Plot settings.
    plot_vars = list()
    fig_labels = list()
    


;---Read AE, SME.
    ae_var = 'sm_sme'
    plot_vars.add, ae_var
    fig_labels.add, 'AE'

    ; This automatically read AE and combine AE and SME in one panel.
    if check_if_update(ae_var, time_range) then supermag_read_sme, time_range
    get_data, ae_var, times, data
    yrange = minmax(data)
    foreach ystep, [50,100,200,500,1000] do begin
        ytickv = make_bins(yrange, ystep)
        yticks = n_elements(ytickv)-1
        if yticks le 3 then break
    endforeach
    yrange = minmax(ytickv)
    yminor = 5

    add_setting, ae_var, smart=1, dictionary($
        'display_type', 'stack', $
        'unit', 'nT', $
        'ystyle', 1, $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'colors', sgcolor(['black','red']), $
        'labels', ['AE','SME'] )



;---Read theta, mlt, r_sm.
    plot_vars.add, 'theta_ewo'
    fig_labels.add, 'Dipolari-!C     zation'

    time_step = 10.
    common_times = make_bins(time_range+[0,-1]*time_step, time_step)
    avail_probes = list()
    foreach probe, probes do begin
        prefix = probe+'_'

        theta_var = prefix+'theta'
        if check_if_update(theta_var, time_range) then begin
            azim_dp_read_theta, time_range, probe=probe
            get_data, theta_var, times, data
            if n_elements(times) le 10 then continue
            interp_time, theta_var, common_times
        endif

        mlt_var = prefix+'mlt'
        if check_if_update(mlt_var, time_range) then begin
            azim_dp_read_mlt, time_range, probe=probe
            get_data, mlt_var, times, data
            if n_elements(times) le 10 then continue
            interp_time, mlt_var, common_times
        endif

        r_sm_var = prefix+'r_sm'
        if check_if_update(r_sm_var, time_range) then begin
            azim_dp_read_orbit, time_range, probe=probe
            if tnames(r_sm_var) eq '' then continue
            interp_time, r_sm_var, common_times
        endif

        scaled_theta_var = prefix+'scaled_theta'
        if check_if_update(scaled_theta_var, time_range) then begin
            theta = get_var_data(theta_var)
            mlt = get_var_data(mlt_var)
            scaled_theta = azim_dp_scale_theta(theta, mlt, width=scale_width)
            store_data, scaled_theta_var, common_times, scaled_theta
        endif
        avail_probes.add, probe
    endforeach



;---Read J up down ewogram.
    j_up_ewo_var = 'thg_j_up_ewo'
    plot_vars.add, j_up_ewo_var
    fig_labels.add, 'Upwd!C    current'
    if check_if_update(j_up_ewo_var, time_range) then begin
        themis_read_upward_current_ewo, time_range, mlt_range=mlt_range, mlat_range=mlat_range
    endif

    j_down_ewo_var = 'thg_j_down_ewo'
    plot_vars.add, j_down_ewo_var
    fig_labels.add, 'Downwd!C    current'
    if check_if_update(j_down_ewo_var, time_range) then begin
        themis_read_downward_current_ewo, time_range, mlt_range=mlt_range, mlat_range=mlat_range
    endif



;---Prepare the plot.
    ; Settings.
    probe_xys = [-1,1]
    line_linestyle = 2
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 1
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1
    panel_xsize_max = 10
    panel_xsize_min = 2


;---Get the size of the figure, and position of the panels.
if keyword_set(test) then magnify = 2
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, wdelete=1
    ; Get the sizes of the left panels.
    margins = [10.,4,8,1]
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1,1,1]
    nypanel = n_elements(ypans)
    fig_labels = letters(nypanel)+') '+fig_labels.toarray()
    ypad = 0.4
    ; xsize of the main panel.
    if n_elements(panel_x_ratio) eq 0 then panel_x_ratio = 1 ; inch per hour.
    panel_xsize = abs(total(time_range*[-1,1]))/3600*panel_x_ratio
    if ~keyword_set(fix_x_ratio) then begin
        panel_xsize <= panel_xsize_max
        panel_xsize >= panel_xsize_min
    endif
    fig_xsize = panel_xsize+(total(margins[[0,2]]))*abs_xchsz
    fig_ysize = panel_ysize*total(ypans)+(total(margins[[1,3]])+ypad*(nypanel-1))*abs_ychsz

;---Common x-axis.
    xrange = time_range
    xticklen_chsz = -0.2
    yticklen_chsz = -0.4
    xminor = 6 ; hr.
    if duration ge 7200. then begin
        xstep = 3600.
    endif else if duration ge 3600 then begin
        xstep = 1800.
    endif else if duration ge 1800 then begin
        xstep = 600.
    endif else step = 60.
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


    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz, magnify=magnify
    poss = sgcalcpos(nypanel, ypans=ypans, ypad=ypad, margins=margins)


    ; Label.
    for ii=0,nypanel-1 do begin
        tpos = poss[*,ii]
        tx = label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, fig_labels[ii]
    endfor



;---AE.
    the_var = 'sm_sme'
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xtickformat = '(A1)'
    yys = float(get_var_data(the_var, in=time_range, times=xxs, limits=lim))
    index = where(abs(yys) ge 1e5, count)
    if count ne 0 then yys[index] = !values.f_nan
    yrange = lim.yrange
    ytickv = lim.ytickv
    yticks = lim.yticks
    yminor = lim.yminor
    ytitle = lim.ytitle

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Add data.
    colors = lim.colors
    oplot, xxs, yys[*,0], color=colors[0]
    oplot, xxs, yys[*,1], color=colors[1]
    
    ; Add labels.
    labels = lim.labels
    tx = tpos[2]+xchsz*1
    dy = (tpos[3]-tpos[1])/3
    for ii=0,1 do begin
        ty = tpos[3]-dy*ii-dy
        xyouts, tx,ty,normal=1, labels[ii], color=colors[ii]
    endfor

    ; Add grid.
    grid_color = sgcolor('gray')
    grid_linestyle = 1
    y_constants = [1000]
    x_constants = make_bins(time_range, xstep, /inner)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


;---Panel b. UT-MLT long.
    tpos = poss[*,1]

    ; Add SME2d.
    sme2d_var = 'sm_sme2d'
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    yrange = mlt_range
    foreach ystep, [3,4,6,8] do begin
        ytickv = make_bins(yrange, ystep, inner=1)
        yticks = n_elements(ytickv)-1
        if yticks le 3 then break
    endforeach
    yminor = ystep
    options, sme2d_var, 'ztitle', 'SME (nT)'
    options, sme2d_var, 'zcharsize', 0.8
    options, sme2d_var, 'zposition', cbpos
    options, sme2d_var, 'xticklen', xticklen
    options, sme2d_var, 'yticklen', yticklen
    options, sme2d_var, 'color_table', 60
    options, sme2d_var, 'yrange', yrange
    options, sme2d_var, 'yticks', yticks
    options, sme2d_var, 'ytickv', ytickv
    options, sme2d_var, 'yminor', yminor    
    tplot, sme2d_var, trange=time_range, position=tpos, nouttick=1, noerase=1



    ; Add spacecraft tracks.
    constant = 0

    data_var_suffix = 'scaled_theta'
    ztitle = 'Scaled '+tex2str('theta')+' (deg)'
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

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase


    ntime = n_elements(common_times)
    nprobe = n_elements(avail_probes)
    thetas = fltarr(ntime,nprobe)
    mlts = fltarr(ntime,nprobe)
    rxys = fltarr(ntime,nprobe)
    rsms = fltarr(ntime,3,nprobe)
    foreach probe, avail_probes, probe_id do begin
        prefix = probe+'_'
        the_var = prefix+'scaled_theta'
        if tnames(the_var) eq '' then continue
        thetas[*,probe_id] = get_var_data(the_var)
        mlts[*,probe_id] = get_var_data(prefix+'mlt')
        rsms[*,*,probe_id] = get_var_data(prefix+'r_sm')
        rxys[*,probe_id] = snorm(reform(rsms[*,0:1,probe_id]))
    endforeach

    ; Draw the DP ewogram.
    foreach time, common_times, time_id do begin
        zzs = thetas[time_id,*]
        yys = mlts[time_id,*]
        rxy = rxys[time_id,*]
        rsm = reform(rsms[time_id,*,*])

        index = where(finite(zzs),count)
        if count eq 0 then continue
        yys = yys[index]
        zzs = azim_df_normalize_theta(zzs[index], zrange=tilt_range, ct=ct, /reverse_ct)
        rxy = rxy[index]
        rsm = rsm[index,*]
        
;        if time ge mean(time_range) then stop

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
    
    ; Add probes.
    ncommon_time = n_elements(common_times)
    foreach probe, avail_probes, probe_id do begin
        prefix = probe+'_'
        the_var = prefix+'scaled_theta'
        if tnames(the_var) eq '' then continue
        thetas = get_var_data(the_var)
        mlts = get_var_data(prefix+'mlt')
        rsms = get_var_data(prefix+'r_sm')
        rxys = snorm(reform(rsms[*,0:1]))
        magns = check_if_in_magn(rsms, dynamic_pressure=pdyn)
        
        mask = fltarr(ncommon_time)+1
        mask = mask and (mlts ge mlt_range[0] and mlts le mlt_range[1])
        mask = mask and (rxys ge rxy_range[0] and rxys le rxy_range[1])
        mask = mask and (magns eq 1)
        
        index = (where(mask eq 1, count))[0]
        if count eq 0 then continue
        
        tx = common_times[index]
        ty = mlts[index]
        xyouts, tx,ty,data=1, strupcase(probe), charsize=0.7
    endforeach


    ; Add grid.
    y_constants = make_bins(yrange, ystep, /inner)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle


;    ; Add box.
;    xtickformat = '(A1)'
;    plot, xrange, yrange, position=tpos, $
;        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
;        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
;        /nodata, /noerase

;    ; Draw color bar.
;    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
;    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;---Panel c, EWOgram of upward current.
    tpos = poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    yrange = mlt_range
    foreach ystep, [3,4,6,8] do begin
        ytickv = make_bins(yrange, ystep, inner=1)
        yticks = n_elements(ytickv)-1
        if yticks le 3 then break
    endforeach
    yminor = ystep


    j_up_ewo_var1 = j_up_ewo_var+'1'
    copy_data, j_up_ewo_var, j_up_ewo_var1
    get_data, j_up_ewo_var1, xx,yy,zz
    yy *= 1e-3
    store_data, j_up_ewo_var1, xx,yy,zz
    options, j_up_ewo_var1, 'ztitle', 'Upward current (kA)'
    options, j_up_ewo_var1, 'zrange', [0,50]
    options, j_up_ewo_var1, 'zcharsize', 0.8
    options, j_up_ewo_var1, 'zposition', cbpos
    j_up_ewo_color_table = 62
    options, j_up_ewo_var1, 'xticklen', xticklen
    options, j_up_ewo_var1, 'yticklen', yticklen
    options, j_up_ewo_var1, 'color_table', j_up_ewo_color_table
    options, j_up_ewo_var1, 'yrange', yrange
    options, j_up_ewo_var1, 'yticks', yticks
    options, j_up_ewo_var1, 'ytickv', ytickv
    options, j_up_ewo_var1, 'yminor', yminor
    tplot, j_up_ewo_var1, trange=time_range, position=tpos, /noerase, /nouttick

    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, noerase=1, nodata=1

    ; Add grid.
    grid_color = sgcolor('gray')
    grid_linestyle = 1
    y_constants = make_bins(yrange, ystep, /inner)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle


;---Panel d, EWOgram of downward current.
    tpos = poss[*,3]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    yrange = mlt_range
    foreach ystep, [3,4,6,8] do begin
        ytickv = make_bins(yrange, ystep, inner=1)
        yticks = n_elements(ytickv)-1
        if yticks le 3 then break
    endforeach
    yminor = ystep


    j_down_ewo_var1 = j_down_ewo_var+'1'
    copy_data, j_down_ewo_var, j_down_ewo_var1
    get_data, j_down_ewo_var1, xx,yy,zz
    yy *= 1e-3
    store_data, j_down_ewo_var1, xx,yy,zz
    options, j_down_ewo_var1, 'ztitle', 'Downward current (kA)'
    options, j_down_ewo_var1, 'zrange', [0,50]
    options, j_down_ewo_var1, 'zcharsize', 0.8
    options, j_down_ewo_var1, 'zposition', cbpos
    j_down_ewo_color_table = 57
    options, j_down_ewo_var1, 'xticklen', xticklen
    options, j_down_ewo_var1, 'yticklen', yticklen
    options, j_down_ewo_var1, 'color_table', j_down_ewo_color_table
    options, j_down_ewo_var1, 'yrange', yrange
    options, j_down_ewo_var1, 'yticks', yticks
    options, j_down_ewo_var1, 'ytickv', ytickv
    options, j_down_ewo_var1, 'yminor', yminor    
    tplot, j_down_ewo_var1, trange=time_range, position=tpos, /noerase

    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, noerase=1, nodata=1

    ; Add grid.
    grid_color = sgcolor('gray')
    grid_linestyle = 1
    y_constants = make_bins(yrange, ystep, inner=1)
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle




    if keyword_set(test) then stop
    sgclose
end




time_range = time_double(['2013-06-07/02:30','2013-06-07/07:00'])
probes = ['tha','thd','the','g13','rbspa','rbspb','g15']
mlt_range = [-9,7]

time_range = time_double(['2014-08-28/09:50','2014-08-28/11:10'])
probes = ['tha','thd','the','g13','rbspb','g15']
mlt_range = [-3,9]

time_range = time_double(['2008-02-29/07:50','2008-02-29/09:00'])
probes = ['th'+['a','c','d','e'],'g'+['11','12']]
mlt_range = [-1,1]*6

;time_range = time_double(['2008-03-20/11:30','2008-03-20/13:00'])
;probes = ['th'+['a','b','c','d','e']]
;mlt_range = [-1,1]*6
;
;; Current not obvious.
;time_range = time_double(['2014-12-26/00:30','2014-12-26/02:00'])
;probes = ['g13','g15','tha','the']
;mlt_range = [-12,0]
;
;; Current does not go along with DP, but start at the same time.
;time_range = time_double(['2016-06-15/03:30','2016-06-15/06:00'])
;probes = ['g13','g14','g15','tha','mms1']
;mlt_range = [-9,6]
;
;; Good example.
;time_range = time_double(['2016-10-13/12:00','2016-10-13/13:30'])
;probes = ['rbspa','rbspb','g13','g14','g15','thd']
;mlt_range = [0,9]
;;
;; Interesting example.
;time_range = time_double(['2017-03-31/02:00','2017-03-31/03:30'])
;probes = ['g13','g15','tha','thd','the','rbspa','rbspb']
;mlt_range = [-9,3]

;time_range = time_double(['2017-03-21/07:00','2017-03-21/09:00'])
;probes = ['tha','thd','the','g13','g15']
;mlt_range = [-6,9]
;mlat_range = [70,80]
;test = 1

time_range = time_double(['2019-08-05/23:00','2019-08-06/00:30'])
probes = ['tha','thd','the','g14','g16']
mlt_range = [-9,3]
mlat_range = [55,70]
test = 1

; Mary Hudson's event.
time_range = time_double(['2017-03-21/07:00','2017-03-21/09:00'])
probes = ['tha','thd','the','g13','g15']
mlt_range = [-5,5]
mlat_range = [55,70]
test = 1

; Mary Hudson's event 2.
time_range = time_double(['2017-09-08/03:00','2017-09-08/05:00'])
probes = ['tha','thd','the','g13','g15','mms1']
mlt_range = [-12,2]
mlat_range = [50,80]
test = 1

; Merkin's event.
time_range = time_double(['2016-08-09/06:00','2016-08-09/20:00'])
time_range = time_double(['2016-08-09/09:00','2016-08-09/10:30'])
probes = ['g13','g14','g15','mms1','rbspa','rbspb','tha','thd','the']
mlt_range = [-4,8]
mlat_range = [50,80]
test = 1


; Christine's event.
time_range = time_double(['2008-03-14/06:00','2008-03-14/08:00'])
mlt_range = [-7,3]
probes = ['tha','thd','the','g11','g12','g13']


; Panov 2018 event.
time_range = time_double(['2008-02-15/04:00','2008-02-15/11:00'])
mlt_range = [-1,1]*10

; Sneha's event.
time_range = time_double(['2014-11-16/01:00','2014-11-16/06:00'])
mlt_range = [-1,1]*10
probes = ['th'+['a','d','e'],'rbsp'+['a','b'],'g'+['13','15']]
rxy_range = [3,20]

; Sneha's STEVE event.
time_range = time_double(['2019-03-28/07:00','2019-03-28/10:00'])
mlt_range = [-2,1]
probes = ['th'+['a','e'],'rbsp'+['a','b'],'g'+['17','15','14']]
probes = ['th'+['a','e'],'g'+['15','17']]
rxy_range = [3,20]

fig_ewogram_of_dp_and_up_down_current2, time_range, test=1, mlt_range=mlt_range, rxy_range=rxy_range, $
    probes=probes, avail_probes=avail_probes
end
