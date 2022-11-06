;+
; Plot dipolarization and flow.
;-


;---Settings.
    test = 0
    test_panels = 0
    plot_file = join_path([srootdir(),'fig_2014_0828_flow_and_dp_hgio_22.pdf'])
    if keyword_set(test) then plot_file = 0


    ; For line plots of Pflux, BBFs, etc.
    mid_time_range = time_double(['2014-08-28/10:00','2014-08-28/10:40'])
    mid_probes = 'th'+['a','d','e']
    mid_colors = constant('rgb')
    data_time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])

;---Check inputs.
    if n_elements(project) eq 0 then project = azim_df_load_project()
    if n_elements(mid_time_range) ne 2 then begin
        errmsg = handle_error('No input mid_time_range ...')
        return
    endif

;---Find the event.
    events = azim_df_search_all_events(project=project)
    foreach event, events do if product(event.time_range-mid_time_range) lt 0 then break
    azim_df_load_basic_data, project=project


    ; For add DP.
    dp_info = dictionary($
        'probes', ['thd','the','tha','g15','rbspb','g13'], $
        'dp_times', time_double('2014-08-28/'+['10:10:40','10:13:20','10:15:15','10:20:25','10:44:55','10:47:25']), $
        'omega', 2.1, $
        'angular_width', 10., $
        'plot_time_range', time_double(['2014-08-28/10:07','2014-08-28/10:55']), $
        'plot_mlt_range', [0,7] )
    _2014_0828_10_load_data
    nprobe = n_elements(dp_info.probes)
    dp_mlts = fltarr(nprobe)
    dp_rxys = fltarr(nprobe)
    ndim = 3
    dp_r_sms = fltarr(nprobe,ndim)
    foreach probe, dp_info.probes, probe_id do begin
        time = dp_info.dp_times[probe_id]
        r_gsm = get_var_data(probe+'_r_gsm', at=time)
        r_sm = cotran(r_gsm, time, 'gsm2sm')
        dp_r_sms[probe_id,*] = r_sm
        dp_mlts[probe_id] = pseudo_mlt(r_sm)
        dp_rxys[probe_id] = snorm(r_sm[0:1])
    endforeach
    dp_info['dp_mlts'] = dp_mlts
    dp_info['dp_rxys'] = dp_rxys
    dp_info['dp_r_sms'] = dp_r_sms


    ; Probe colors.
    top_color = 240
    bottom_color = 90
    probe_colors = smkarthm(bottom_color,top_color,nprobe, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)


;---Calc panel positions.
    sgopen, 0, xsize=1, ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ypad = 0.4
    margins = [8,4,8,1]


    second_color = sgcolor('silver')
    test_times = time_double(['2014-08-28/10:10','2014-08-28/10:20'])


    ; Figure.
    fig_xsize = 5
    fig_ysize = 3


    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz
    nmid_pan = 3
    mid_poss = sgcalcpos(nmid_pan, ypad=ypad, margins=margins)

    if keyword_set(test_panels) then begin
        for ii=0,nmid_pan-1 do begin
            tpos = mid_poss[*,ii]
            plot, [0,1],[0,1], $
                xstyle=1, ystyle=1, xtitle='X', ytitle='Y', $
                position=tpos, nodata=1, noerase=1
        endfor

        if keyword_set(test) then stop
        sgclose
        stop
    endif

    xticklen_chsz = -0.20
    yticklen_chsz = -0.40



;---Middle panels.
    time_step = 3
    common_times = make_bins(mid_time_range, time_step)
    ntime = n_elements(common_times)
    nmid_probe = n_elements(mid_probes)

    ; Pflux and BBF.
    foreach probe, mid_probes, probe_id do begin
        prefix = probe+'_'

;        ; Pflux.
;        get_data, prefix+'pf_fac', times, pffac
;        pffac = get_var_data(prefix+'pf_fac', at=common_times)
;        cmap = get_var_data(prefix+'cmap', at=common_times)
;        pf_para_norm = cmap*pffac[*,0]
;        store_data, prefix+'pflux', common_times, pf_para_norm

        ; BBF.
        u_var = prefix+'u_gsm'
        if check_if_update(u_var, data_time_range) then begin
            the_probe = strmid(probe,2,1)
            themis_read_ion_vel, data_time_range, probe=the_probe
        endif
        u_gsm = get_var_data(u_var, at=common_times)
        r_gsm = get_var_data(prefix+'r_gsm', at=common_times)
        bbf = -vec_dot(u_gsm, sunitvec(r_gsm))
        store_data, prefix+'bbf', common_times, bbf
    endforeach




    ; Start to plot.
    xrange = mid_time_range
    xstep = 60*10
    xtickv = make_bins(xrange, xstep)
    xtickn = time_string(xtickv,tformat='hh:mm')
;    xtickn[0] += '!C'+time_string(xtickv[0],tformat='YYYY-MM-DD')
    xtickn[0] = time_string(xtickv[0],tformat='YYYY!CMTH DD')
    xticks = n_elements(xtickv)-1
    xminor = 6


    foreach probe, 'th'+['d','e','a'], probe_id do begin
        prefix = probe+'_'
        tpos = mid_poss[*,probe_id]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        ; theta.
        yys = get_var_data(prefix+'theta', times=xxs)
        case probe of
            'thd': begin
                yrange = [-1,1]*30
                end
            'the': begin
                yrange = [-1,1]*25
                end
            'tha': begin
                yrange = [-1,1]*20
                end
        endcase
        ystep = 5
        ytickv = [-1,0,1]*ceil(yrange[0]/ystep)*ystep
        yticks = n_elements(ytickv)-1
        yminor = 5
        xtickformat='(A1)'
        if probe_id eq nmid_probe-1 then xtickformat=''
        plot, xxs, yys, $
            xstyle=1, xrange=xrange, xlog=0, xtickv=xtickv, xtickname=xtickn, $
            xticks=xticks, xticklen=xticklen, xtickformat=xtickformat, xminor=xminor, $
            ystyle=9, yrange=yrange, ylog=0, $
            ytickv=ytickv, yticks=yticks, yticklen=yticklen, yminor=yminor, $
            position=tpos, noerase=1, nodata=0
;        plots, xrange, [0,0], linestyle=1


        yys = get_var_data(prefix+'bbf', times=xxs)
        ytitle = '(km/s)'
        case probe of
            'thd': yrange = [-1,2]*200
            'the': yrange = [-1,2]*200
            'tha': yrange = [-1,2]*200
        endcase
        axis, yaxis=1, save=1, ystyle=1, yrange=yrange, ylog=0, $
            yticks=yticks, yticklen=yticklen, yminor=yminor, color=second_color
        plots, xxs, yys, color=second_color
        plots, xrange, [0,0], linestyle=1

        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*0.9
        xyouts, tx,ty, /normal, alignment=0, 'a-'+string(probe_id+1,format='(I0)')+'. '+strupcase(probe)
    endforeach

    ytitle = 'Tilt angle theta (deg)'
    tx = tpos[0]-ychsz*1.5
    ty = (mid_poss[1,0]+mid_poss[3,2])*0.5
    xyouts, tx,ty,/normal, alignment=0.5, orientation=90, ytitle

    ytitle = 'Inward ion vel (km/s)'
    tx = tpos[2]+ychsz*1.8
    xyouts, tx,ty,/normal, alignment=0.5, orientation=90, ytitle, color=second_color





    tpos = mid_poss[*,0]
    tpos[1] = mid_poss[1,-1]
    yrange = [0,1]
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
        position=tpos, noerase=1, nodata=1
    foreach tx, test_times do begin
        plots, tx+[0,0], yrange, linestyle=1
    endforeach




    if keyword_set(test) then stop
    sgclose

end
