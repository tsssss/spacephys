function micro_injection_themis_2009_0619_fig_sc_location_v01, test = test
  compile_opt idl2

  info = dictionary()
  info['a'] = dictionary( $
    'color', 'red', $
    'time_range', time_double('2009-06-19/' + ['09:00', '13:00']))
  info['b'] = dictionary( $
    'color', 'green', $
    'time_range', time_double('2009-06-19/' + ['02:00', '06:00']))
  info['c'] = dictionary( $
    'color', 'blue', $
    'time_range', time_double('2009-06-19/' + ['01:00', '06:00']))
  info['d'] = dictionary( $
    'color', 'purple', $
    'time_range', time_double('2009-06-19/' + ['01:00', '04:00']))
  info['e'] = dictionary( $
    'color', 'orange', $
    'time_range', time_double('2009-06-19/' + ['01:00', '04:00']))

  xrange = [2, -4]
  yrange = [16, 10]
  ; xrange = [15,-5]
  ; yrange = [15,-15]
  plot_dir = srootdir()
  base = 'micro_injection_themis_2009_0619_fig_sc_location_v01.pdf'
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
    nodata = 1, noerase = 1, position = tpos

  probes = info.keys()
  dt = 3600d
  pad_time = 600d
  set_circ, fill = 1
  symsize = 0.5
  foreach probe, probes do begin
    my_info = info[probe]
    tr = my_info['time_range']
    color = sgcolor(my_info['color'])
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

  ; ---Magnetopause.
  thick = keyword_set(test) ? 2 : 8
  tts = smkarthm(0, 2 * !dpi, 50, 'n')
  rrs = 2
  xxs = rrs * cos(tts)
  yys = rrs * cos(tts)
  zzs = fltarr(n_elements(xxs))
  time_range = []
  foreach probe, probes do begin
    my_info = info[probe]
    tr = my_info['time_range']
    time_range = [time_range, tr]
  endforeach
  time_range = minmax(time_range)
  pdyn_var = omni_read_sw_p(time_range)
  pdyn = var_get_data(pdyn_var, times = times)
  max_pdyn = max(pdyn, index)
  test_time = times[index]
  pdyn = var_get_data(pdyn_var, at = test_time)
  mpause_t96, max_pdyn, xmgnp = xmgnp, ymgnp = ymgnp, zmgnp = zmgnp, $
    xgsm = xxs, ygsm = yys, zgsm = zzs, id = id, distan = distan
  oplot, xmgnp, ymgnp, linestyle = 2 ; 0, color=color, thick=thick
  ; index = where_pro(ymgnp, '[]', minmax(yrange))
  ; xmgnp = xmgnp[index]
  ; ymgnp = ymgnp[index]
  ; tmp = max(ymgnp, index)
  ; tmp = convert_coord(xmgnp[index],ymgnp[index], data=1, to_normal=1)
  ; tx = tmp[0]+xchsz*1
  ; ty = tmp[1]+ychsz*0
  ; msg = 'Magnetopause'
  ; xyouts, tx,ty,msg,normal=1

  times = make_bins(time_range, 3600d, inner = 1)
  times = [test_time, time_double(['2009-06-19/02:00'])]
  ; ntime = n_elements(times)
  ; colors = get_color(ntime, color_table=64)
  colors = sgcolor(['red', 'blue'])
  foreach test_time, times, tid do begin
    color = colors[tid]
    pdyn = var_get_data(pdyn_var, at = test_time)
    mpause_t96, pdyn, xmgnp = xmgnp, ymgnp = ymgnp, zmgnp = zmgnp, $
      xgsm = xxs, ygsm = yys, zgsm = zzs, id = id, distan = distan
    oplot, xmgnp, ymgnp, linestyle = 3, color = color, thick = thick

    index = where_pro(ymgnp, '[]', minmax(yrange))
    xmgnp = xmgnp[index]
    ymgnp = ymgnp[index]
    tmp = max(ymgnp, index)
    tmp = convert_coord(xmgnp[index], ymgnp[index], data = 1, to_normal = 1)
    tx = tmp[0] + xchsz * 2
    ty = tmp[1] + ychsz * 0.5
    ; msg = 'Magnetopause!Cat '+time_string(test_time,tformat='hh:mm')+' UT'
    msg = 'M/Pause!Cat ' + time_string(test_time, tformat = 'hh:mm') + ' UT'
    xyouts, tx, ty, msg, normal = 1, color = color
  endforeach

  ; ;---Add circles.
  ; rrs = smkarthm(9,15,2,'dx')
  ; color = sgcolor('silver')
  ; tmp = smkarthm(0,2*!dpi,80,'n')
  ; foreach rr, rrs do begin
  ; oplot, rr*cos(tmp), rr*sin(tmp), color=color, linestyle=0
  ; endforeach

  ; ---Draw axis.
  uniform_ticklen = -ychsz * fig_size[0] * 0.15
  xticklen = uniform_ticklen / (tpos[3] - tpos[1]) / fig_size[1]
  yticklen = uniform_ticklen / (tpos[2] - tpos[0]) / fig_size[0]

  plot, xrange, yrange, $
    xstyle = 1, xrange = xrange, xtitle = xtitle, $
    ystyle = 1, yrange = yrange, ytitle = ytitle, $
    xticklen = xticklen, yticklen = yticklen, $
    nodata = 1, noerase = 1, position = tpos

  msg = 'a. SC Position'
  tx = tpos[0] + xchsz * 0.5
  ty = tpos[3] - ychsz * 1
  xyouts, tx, ty, normal = 1, msg

  ; ---Right panels.
  probe = 'a'
  prefix = 'th' + probe + '_'
  my_info = info[probe]
  time_range = time_double('2009-06-19/' + ['09:00', '11:00'])

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
  zrange = [1e1, 1e5]
  ztickv = [1e1, 1e2, 1e3, 1e4, 1e5]
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

  margins = [6, 0, 8, 0]
  plot_vars = [en_high_var, pa_high_var, en_low_var, b_gsm_var, u_gsm_var]
  nplot_var = n_elements(plot_vars)
  panel_labels = letters([0, nplot_var] + 1) + '. ' + ['e- EN high', 'e- PA high', 'e- EN low', 'B GSM', 'V!Dion!N GSM']
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

  options, b_gsm_var, yrange = [-7, 17], ytickv = [-5, 5, 15], yticks = 2, yminor = 2, constant = [0]
  options, u_gsm_var, yrange = [-1, 1] * 85, ytickv = [-80, 0, 80], yticks = 2, yminor = 4, constant = [0]

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

  ; Add microinjection times.
  mi_times = time_double([ $
    ; '2009-06-19/10:13:39', $
    '2009-06-19/09:44:30', $
    '2009-06-19/09:51:39', $
    '2009-06-19/10:14:40', $
    '2009-06-19/10:23:30', $
    '2009-06-19/10:31:20', $
    ; '2009-06-19/10:46:39', $
    ; '2009-06-19/10:22:38', $
    ; '2009-06-19/10:31:41', $
    '2009-06-19/10:38:21'])
  timebar, mi_times, color = sgcolor('red'), linestyle = 2

  if keyword_set(test) then stop
  sgclose

  return, plot_file
end

compile_opt idl2
test = 0
print, micro_injection_themis_2009_0619_fig_sc_location_v01(test = test)
end
