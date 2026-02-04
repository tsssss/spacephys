function micro_injection_themis_2009_0619_fig_thd_the_v01, test = test

    info = dictionary()
    info['a'] = dictionary( $
        'color', 'red', $
        'time_range', time_double('2009-06-19/' + ['09:00', '13:00']))
    info['b'] = dictionary( $
        'color', 'green', $
        'time_range', time_double('2009-06-19/' + ['01:00', '06:00']))
    info['c'] = dictionary( $
        'color', 'blue', $
        'time_range', time_double('2009-06-19/' + ['01:00', '06:00']))
    info['d'] = dictionary( $
        'color', 'purple', $
        'time_range', time_double('2009-06-19/' + ['01:00', '04:00']))
    info['e'] = dictionary( $
        'color', 'violet', $
        'time_range', time_double('2009-06-19/' + ['01:00', '04:00']))

    xrange = [2, -3]
    yrange = [16.5, 10]
    plot_dir = srootdir()
    base = 'micro_injection_themis_2009_0619_fig_thd_the_v01.pdf'
    plot_file = join_path([plot_dir, base])
    if keyword_set(test) then plot_file = 0
    margins = [8, 4, 2, 1]
    pansize = abs([(xrange[1] - xrange[0]), (yrange[1] - yrange[0])])
    pansize = pansize / pansize[0] * 3.5
    nypan = 1
    nxpan = 2
    xpans = [1, 2]
    poss = panel_pos(plot_file, nypan = nypan, nxpan = nxpan, margins = margins, $
        xpans = xpans, pansize = pansize, fig_size = fig_size)
    sgopen, plot_file, size = fig_size, xchsz = xchsz, ychsz = ychsz

    tpos = poss[*, 0]
    coord = 'gsm'
    xtitle = strupcase(coord) + ' X (Re)'
    ytitle = strupcase(coord) + ' Y (Re)'
    plot, xrange, yrange, $
        xstyle = 5, xrange = xrange, xtitle = xtitle, $
        ystyle = 5, yrange = yrange, ytitle = ytitle, $
        nodata = 1, noerase = 1, position = tpos, iso = 1

    probes = info.keys()
    probes = ['a','d','e']
    gray_probes = ['a']
    dt = 3600d
    pad_time = 600d
    set_circ, fill = 1
    symsize = 0.5
    foreach probe, probes do begin
      my_info = info[probe]
      tr = my_info['time_range']
      color = sgcolor(my_info['color'])
      index = where(gray_probes eq probe, count)
      if count ne 0 then color = sgcolor('silver')
      r_var = themis_read_orbit(tr + [-1, 1] * dt, probe = probe, coord = coord)

      tts = make_bins(tr, pad_time)
      r_vec = var_get_data(r_var, at = tts)
      xxs = r_vec[*, 0]
      yys = r_vec[*, 1]
      plots, xxs, yys, color = color

      tts = make_bins(tr, dt)
      r_vec = var_get_data(r_var, at = tts)
      xxs = r_vec[*, 0]
      yys = r_vec[*, 1]
      msgs = time_string(tts, tformat = 'hh')
      msgs[-1] += ' UT'
      plots, xxs, yys, color = color, psym = 8, symsize = symsize
      ntt = n_elements(tts)
      foreach time, tts, tid do begin
        tmp = convert_coord(xxs[tid], yys[tid], data = 1, to_normal = 1)
        if probe eq 'a' then begin
          tx = tmp[0] - xchsz * 2
          ty = tmp[1] - ychsz * 0.8
        endif else if probe eq 'b' then begin
          tx = tmp[0] - xchsz * 1
          ty = tmp[1] + ychsz * 0.3
        endif else if probe eq 'c' then begin
          tx = tmp[0] - xchsz * 1
          ty = tmp[1] + ychsz * 0.3
        endif else if probe eq 'd' then begin
          tx = tmp[0]
          ty = tmp[1] + ychsz * 0.3
        endif else if probe eq 'e' then begin
          tx = tmp[0]
          ty = tmp[1] - ychsz * 1
        endif
        alignment = 0.5

        msg = msgs[tid]
        xyouts, tx, ty, msg, color = color, normal = 1, alignment = alignment

        if tid eq ntt - 1 then begin
          msg = 'TH-' + strupcase(probe)
          tx = tmp[0] + xchsz * 2.5
          ty = tmp[1] - ychsz * 0.5
          if probe eq 'e' then tx = tmp[0] + xchsz * 3.5
          xyouts, tx, ty, msg, color = color, normal = 1, alignment = 0.5
        endif
      endforeach
    endforeach

;---Magnetopause.
    thick = keyword_set(test) ? 2 : 8
    tts = smkarthm(0, 2 * !dpi, 50, 'n')
    rrs = 2
    xxs = rrs * cos(tts)
    yys = rrs * cos(tts)
    zzs = fltarr(n_elements(xxs))
    time_range = []
    foreach probe, probes do begin
        index = where(gray_probes eq probe, count)
        if count ne 0 then continue
        my_info = info[probe]
        tr = my_info['time_range']
        time_range = [time_range, tr]
    endforeach
    time_range = minmax(time_range)
    pdyn_var = omni_read_sw_p(time_range)
    pdyn = var_get_data(pdyn_var, times=times, in=time_range)
    test_times = []
    max_pdyn = max(pdyn, index)
    test_times = [test_times, times[index]]
    min_pdyn = min(pdyn, index)
    test_times = [test_times, times[index]]

    the_probe = 'd'
    my_info = info[the_probe]
    color = sgcolor(my_info['color'])
    foreach test_time, test_times do begin
        pdyn = var_get_data(pdyn_var, at=test_time)
        mpause_t96, pdyn, xmgnp=xmgnp, ymgnp=ymgnp, zmgnp=zmgnp, $
            xgsm=xxs, ygsm=yys, zgsm=zzs, id=id, distan=distan
        oplot, xmgnp, ymgnp, linestyle=2, color=color
    endforeach

    mean_pdyn = mean([min_pdyn,max_pdyn])
    mpause_t96, mean_pdyn, xmgnp=xmgnp, ymgnp=ymgnp, zmgnp=zmgnp, $
        xgsm=xxs, ygsm=yys, zgsm=zzs, id=id, distan=distan
    index = where_pro(ymgnp, '[]', minmax(yrange))
    xmgnp = xmgnp[index]
    ymgnp = ymgnp[index]
    tmp = max(ymgnp, index)
    tmp = convert_coord(xmgnp[index], ymgnp[index], data = 1, to_normal = 1)
    tx = tmp[0]-xchsz*1
    ty = tmp[1]
    msg = 'Magnetopause'
    xyouts, tx, ty, msg, normal = 1, color = color, alignment = 0.5
    

;---Draw axis.
    uniform_ticklen = -ychsz * fig_size[0] * 0.15
    xticklen = uniform_ticklen / (tpos[3] - tpos[1]) / fig_size[1]
    yticklen = uniform_ticklen / (tpos[2] - tpos[0]) / fig_size[0]

    xstep = 1
    xtickv = make_bins(xrange, xstep, inner = 1)
    xticks = n_elements(xtickv) - 1
    xminor = 2
    ystep = 1
    ytickv = make_bins(yrange, ystep, inner = 1)
    yticks = n_elements(ytickv) - 1
    yminor = 2
    plot, xrange, yrange, $
        xstyle = 1, xrange = xrange, xtitle = xtitle, xminor = xminor, xtickv = xtickv, xticks = xticks, $
        ystyle = 1, yrange = yrange, ytitle = ytitle, yminor = yminor, ytickv = ytickv, yticks = yticks, $
        xticklen = xticklen, yticklen = yticklen, $
        nodata = 1, noerase = 1, position = tpos, iso = 1

    msg = 'a. SC Position'
    tx = tpos[0] + xchsz * 0.5
    ty = tpos[3] - ychsz * 1
    xyouts, tx, ty, normal = 1, msg




;---Right panels.
    probes = ['d','e']
    foreach probe, probes do begin
        prefix = 'th' + probe + '_'
        my_info = info[probe]
        time_range = my_info['time_range']

        ; Load data.
        mission_probe = 'th' + probe

        ; B field and ion velocity.
        b_gsm_var = themis_read_bfield(time_range, probe = probe, errmsg = errmsg, id = 'fgs')
        u_gsm_var = themis_read_ion_vel(time_range, probe = probe, errmsg = errmsg, id = 'peir')

        ; SST.
        datatype = 'psef'
        prefix2 = prefix + datatype + '_'
        en_high_var = prefix2 + 'en_eflux'
        pa_high_var = prefix2 + 'an_eflux_pa'
        if check_if_update(en_high_var, time_range) then begin
            thm_part_load, data_type = datatype, probe = probe, trange = time_range
            thm_part_getspec, data_type = datatype, probe = probe, trange = time_range, outputs = 'energy'
            thm_part_getspec, data_type = datatype, probe = probe, trange = time_range, outputs = 'pa'
            options, en_high_var, requested_time_range = time_range
            options, pa_high_var, requested_time_range = time_range
        endif

        unit = 'eV/cm!E2!N-s-sr-eV'
        zrange = [1e2, 1e7]
        ztickv = [1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
        ; zrange = [1e1,1e5]
        ; ztickv = [1e1,1e2,1e3,1e4,1e5]
        ztickv_log = alog10(ztickv)
        zticks = n_elements(ztickv) - 1
        ztickn = '10!U' + string(ztickv_log, format = '(I0)')
        ztickn[0 : * : 2] = ' '
        options, [en_high_var], color_table = 40, no_interp = 1, $
        ytitle = 'Energy!C(eV)', ztitle = unit, $
        zrange = zrange, zstyle = 1, zlog = 1, ztickv = ztickv, ztickname = ztickn, zticks = zticks, zminor = 9, $
        yrange = [3.1e4, 7.2e5], ystyle = 1, ylog = 1, ytickv = [5e4, 5e5], ytickname = '10!U' + ['4', '5'], yticks = 1, yminor = 9
        options, [pa_high_var], color_table = 40, no_interp = 1, $
        ytitle = 'PA!C(deg)', ztitle = unit, $
        zrange = zrange, zstyle = 1, zlog = 1, ztickv = ztickv, ztickname = ztickn, zticks = zticks, zminor = 9, $
        yrange = [0, 180], ystyle = 1, ylog = 0, ytickv = [30, 90, 150], ytickname = ['30', '90', '150'], yticks = 2, yminor = 6
        foreach var, [en_high_var, pa_high_var] do begin
        add_setting, var, smart = 1, dictionary( $
            'display_type', 'spec')
        endforeach

        ; ESA.
        datatype = 'peef'
        prefix2 = prefix + datatype + '_'
        zrange = [1e5, 1e8]
        ztickv = [1e5, 1e6, 1e7, 1e8]
        ztickn = '10!U' + ['5', '6', '7', '8']
        ; ztickn[0:2:*] = ' '
        zticks = n_elements(ztickv) - 1
        en_low_var = themis_read_en_spec(time_range, probe = probe, species = 'e', id = 'esa_l2')
        options, en_low_var, color_table = 40, no_interp = 1, $
        zrange = zrange, zstyle = 1, zlog = 1, ztickv = ztickv, ztickname = ztickn, zticks = zticks, zminor = 9, $
        yrange = [1.1e1, 2.6e4], ystyle = 1, ylog = 1, ytickv = [1e2, 1e3, 1e4], ytickname = '10!U' + ['2', '3', '4'], yticks = 2, yminor = 9
;    
;        if probe eq 'e' then begin
;            dt = -280d
;            data = var_get_data(en_low_var, val, times=times)
;            en_low_var = var_store(en_low_var, data, times+dt, val)
;        endif
    endforeach

    margins = [6, 0, 8, 0]
    en_low_var = 'th'+probes+'_e_en_spec'
    b_gsm_var = 'th'+probes+'_b_gsm'
    plot_vars = []
    foreach probe, probes do begin
        prefix = 'th'+probe+'_'
        en_low_var = prefix+'e_en_spec'
        b_gsm_var = prefix+'b_gsm'
        plot_vars = [plot_vars, en_low_var,b_gsm_var]
    endforeach
    nplot_var = n_elements(plot_vars)
    suffix = '!C    TH-'+strupcase(probes)
    labels = []
    foreach probe, probes, pid do begin
        labels = [labels,['e- EN low','B GSM']+suffix[pid]]
    endforeach
    panel_labels = letters([0, nplot_var] + 1) + '. ' + labels
    right_pos = poss[*, 1]
    panel_poss = sgcalcpos(nplot_var, margins = margins, region = right_pos)

    uniform_ticklen = -ychsz * fig_size[0] * 0.15
    for pid = 0, nplot_var - 1 do begin
        tpos = panel_poss[*, pid]
        xticklen = uniform_ticklen / (tpos[3] - tpos[1]) / fig_size[1]
        yticklen = uniform_ticklen / (tpos[2] - tpos[0]) / fig_size[0]
        var = plot_vars[pid]
        options, var, xticklen = xticklen, yticklen = yticklen
        is_spec = var_get_setting(var, 'spec')
        if is_spec then begin
            zticklen = -0.5
            options, var, zticklen = zticklen
        endif
    endfor

    options, b_gsm_var, yrange = [-20,25], ytickv = [-10,5,20], yticks=2, yminor=3, constant=[0]

    tplot, plot_vars, trange = time_range, position = panel_poss, noerase = 1
    for pid = 0, nplot_var - 1 do begin
      tpos = panel_poss[*, pid]
      tx = tpos[0] - xchsz * 12
      ty = tpos[3] - ychsz * 0.8
      msg = panel_labels[pid]
      xyouts, tx, ty, msg, normal = 1
    endfor

    ; Add label.
    tpos = panel_poss[*, 0]
    tx = tpos[0] + xchsz * 0.5
    ty = tpos[3] - ychsz * 1
    msg = 'TH-' + strupcase(probe)
    xyouts, tx, ty, msg, normal = 1, color = sgcolor('white')

;    ; Add microinjection times.
;    mi_times = time_double([ $
;      '2009-06-19/09:44:30', $
;      '2009-06-19/09:51:39', $
;      '2009-06-19/10:14:40', $
;      '2009-06-19/10:23:30', $
;      '2009-06-19/10:31:20', $
;      '2009-06-19/10:38:21'])
;    timebar, mi_times, color = sgcolor('red'), linestyle = 2


;---Calc max correlation and dt.
    dt = -280d

    
    pid = where(plot_vars eq 'thd_b_gsm', count)
    if count ne 0 then begin
        bvecs = var_get_data('the_b_gsm', times=times, settings=settings)
        xrange = time_range
        yrange = settings['yrange']
        tpos = panel_poss[*,pid]
        set_axis, plot_vars[pid], xrange=xrange, yrange=yrange, position=tpos
        txs = times+dt
        tys = bvecs[*,2]
        oplot, txs, tys, color=sgcolor('purple')
    endif
    
    linestyle = 1
    color = sgcolor('red')
    dts = [0,dt]
    sample_energys = [1e3,1e4]
    bar_times = time_double([ $
        '2009-06-19/02:45:51', $
        '2009-06-19/02:54:11'])
    foreach bar_time, bar_times do begin
        txs = []
        tys = []
        foreach probe, probes, probe_id do begin
            prefix = 'th'+probe+'_'
            the_var = prefix+'e_en_spec'
            pid = where(plot_vars eq the_var, count)
            if count eq 0 then message, 'Variable not found'
            xrange = time_range
            yrange = var_get_setting(the_var, 'yrange')
            tpos = panel_poss[*,pid]
            set_axis, the_var, xrange=xrange, yrange=yrange, position=tpos
            tx = bar_time-dts[probe_id]
            ty = sample_energys[probe_id]
            tmp = convert_coord(tx,ty, data=1, to_normal=1)
            txs = [txs, tmp[0]]
            tys = [tys, tmp[1]]
        endforeach
        hsize = keyword_set(test)? 4: 80
        arrow, txs[0],tys[0],txs[1],tys[1],normal=1, color=color, $
            linestyle=linestyle, hsize=hsize, solid=1
    endforeach
    ;    foreach probe, probes, probe_id do begin
;        prefix = 'th'+probe+'_'
;        the_var = prefix+'e_en_spec'
;
;        pid = where(plot_vars eq the_var, count)
;        if count ne 0 then begin
;            xrange = time_range
;            yrange = var_get_setting(the_var, 'yrange')
;            tpos = panel_poss[*,pid]
;
;            foreach tx, bar_times do begin
;                x0 = tx
;                y0 = 1e3
;                set_axis, the_var, xrange=xrange, yrange=yrange, position=tpos
;                tmp = convert_coord(x0,y0, data=1, to_normal)
;                
;                y1 = y0*10
;                x1 = tx-dts[probe_id]
;                set_axis, the_var, xrange=xrange, yrange=yrange, position=tpos
;
;                ;plots, tx+[0,0], yrange, linestyle=linestyle, color=color
;                arrow, 
;            endforeach
;        endif
;    endforeach
    

    if keyword_set(test) then stop
    sgclose

    return, plot_file
end

test = 1
print, micro_injection_themis_2009_0619_fig_thd_the_v01(test = test)
end
