;+
; Low lattiude TEC enhancement events.
; May be related to the E field of over shielding.
; Check if RBSP observed any E field changes. 
;-

test = 0
time_range = time_double(['2013-03-17/05:00','2013-03-17/15:00'])

probes = ['a','b']
probe_colors = sgcolor(['red','blue'])
probe_psyms = [6,8]
symthick = keyword_set(test)? 2:4
symsize = 0.8
tmp = smkarthm(0,2*!dpi,30, 'n')
txs = cos(tmp)
tys = sin(tmp)
usersym, txs, tys, thick=symthick

event_times = time_double('2013-03-17/'+['05:58','07:40','08:00','08:35','08:50','09:15','09:45','11:15','12:05','12:50'])
nevent = n_elements(event_times)
event_colors = smkarthm(50,250, nevent, 'n')
event_colors = sgcolor(event_colors, ct=40)

orbit_xrange = [2,-6.]
orbit_yrange = [3,-4.]
xspan = abs(total(orbit_xrange*[-1,1]))
yspan = abs(total(orbit_yrange*[-1,1]))


foreach probe, probes do begin
    rbsp_efw_read_l3, time_range, probe=probe
endforeach

plot_vars = []
re = constant('re')
foreach probe, probes do begin
    prefix = 'rbsp'+probe+'_efw_'
    rbspx = strupcase('rbsp'+probe)
    get_data, prefix+'position_gse', times, r_gse
    
    dis_var = prefix+'dis'
    store_data, dis_var, times, snorm(r_gse)/re
    add_setting, dis_var, smart=1, dictionary($
        'short_name', '|R|', $
        'display_type', 'scalar', $
        'unit', 'Re' )
    
    r_sm = cotran(r_gse, times, 'gse2sm')
    for ii=0,2 do r_sm[*,ii] /= re
    r_sm_var = 'rbsp'+probe+'_r_sm'
    store_data, r_sm_var, times, r_sm
    add_setting, r_sm_var, smart=1, dictionary($
        'short_name', 'R', $
        'display_type', 'vector', $
        'unit', 'Re', $
        'coord', 'SM' )
    
    mlt_var = prefix+'mlt'
    add_setting, mlt_var, smart=1, dictionary($
        'short_name', 'MLT', $
        'display_type', 'scalar', $
        'unit', 'h', $
        'ystyle', 1, $
        'yticks', 2, $
        'yminor', 6, $
        'yrange', [-1,1]*12 )
    
    e_var = prefix+'efield_in_corotation_frame_spinfit_edotb_mgse'
    add_setting, e_var, smart=1, dictionary($
        'short_name', 'E', $
        'display_type', 'vector', $
        'coord', 'MGSE', $
        'yrange', [-1,1]*10, $
        'unit', 'mV/m' )
    
    ; convert E from mGSE to SM
    get_data, e_var, times, e_mgse, limits=lim
    e_sm = cotran(e_mgse, times, 'mgse2sm', probe=probe)
    e_var = prefix+'efield_in_corotation_frame_spinfit_edotb_sm'
    store_data, e_var, times, e_sm
    add_setting, e_var, smart=1, dictionary($
        'short_name', 'E', $
        'display_type', 'vector', $
        'coord', 'SM', $
        'yrange', [-1,1]*10, $
        'constant', 0, $
        'unit', 'mV/m' )
    ;plot_vars = [plot_vars, [dis_var,mlt_var,e_var]]
    
    
    ; Read B in SM.
    prefix1 = 'rbsp'+probe+'_'
    rbsp_read_bfield, time_range, probe=probe, coord='sm'
    b_var = prefix1+'b_sm'
    b_sm = get_var_data(b_var, at=times)
    v_sm = vec_cross(e_sm, b_sm)
    coef = 1d3/snorm(b_sm)^2
    for ii=0,2 do v_sm[*,ii] *= coef
    v_sm_var = prefix1+'vexb_sm'
    store_data, v_sm_var, times, v_sm
    add_setting, v_sm_var, smart=1, dictionary($
        'short_name', 'ExB V', $
        'display_type', 'vector', $
        'coord', 'SM', $
        'constant', 0, $
        'yrange', [-1,1]*40, $
        'yticks', 4, $
        'yminor', 2, $
        'unit', 'km/s' )
    plot_vars = [plot_vars, [v_sm_var]]
endforeach

plot_file = join_path([srootdir(),'check_efield_for_tec_event_'+time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_v03.pdf'])
if keyword_set(test) then plot_file = 0

margins = [10,4,8,1]
nvar = n_elements(plot_vars)+1
ypans = [1,fltarr(nvar-1)+0.4]
ypads = fltarr(nvar-1)+0.4
ypads[0] = 5
pansize = [xspan,yspan]
pansize = pansize/pansize[0]*3
poss = panel_pos(plot_file, margins=margins, $
    nxpan=1, nypan=nvar, ypans=ypans, $
    pansize=pansize, panid=0, ypads=ypads, $
    fig_size=fig_size)
    
sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
label_size = 0.8

    
tplot_options, 'labflag', -1
tposs = poss[*,1:nvar-1]
tplot, plot_vars, trange=time_range, position=tposs, noerase=1

pan_letters = letters(3)
foreach plot_var, plot_vars, ii do begin
    tpos = tposs[*,ii]
    tx = tpos[0]-xchsz*8
    ty = tpos[3]-ychsz*0.7
    probe = strmid(plot_var, 4,1)
    msg = pan_letters[ii+1]+') RB-'+strupcase(probe)
    xyouts, tx, ty, normal=1, msg
endforeach

tpos = poss[*,1]
tpos[1] = poss[1,nvar-1]
xrange = time_range
yrange = [0,1]
plot, xrange, yrange, nodata=1, noerase=1, position=tpos, $
    xstyle=5, ystyle=5, xrange=xrange, yrange=yrange

    
foreach event_time, event_times, event_id do begin
    color = event_colors[event_id]

    tmp = convert_coord(event_time,0, data=1, to_normal=1)
    tx = tmp[0]
    tys = tpos[[1,3]]
    plots, tx+[0,0], tys, normal=1, color=color, linestyle=0
    
    msg = string(event_id+1,format='(I0)')
    xyouts, tx,tys[1]+ychsz*0.5, msg, normal=1, alignment=0.5, charsize=label_size, color=color
endforeach


tpos = poss[*,0]
xrange = orbit_xrange
yrange = orbit_yrange
xtitle = 'SM X (Re)'
ytitle = 'SM Y (Re)'
xstep = 2
xminor = 4
xtickv = make_bins(xrange, xstep, inner=1)
ystep = 2
yminor = 4
ytickv = make_bins(yrange, ystep, inner=1)
xticks = n_elements(xtickv)-1
yticks = n_elements(ytickv)-1
xticklen = -0.02
yticklen = -0.02

plot, xrange, yrange, $
    xstyle=5, ystyle=5, $
    xrange=xrange, yrange=yrange, $
    nodata=1, noerase=1, position=tpos

    ; Add earth.
    plots, xrange, [0,0], linestyle=1
    plots, [0,0], yrange, linestyle=1
    tmp = smkarthm(0,2*!dpi, 40, 'n')
    circ_xs = cos(tmp)
    circ_ys = sin(tmp)
    polyfill, circ_xs>0, circ_ys, color=sgcolor('white')
    polyfill, circ_xs<0, circ_ys, color=sgcolor('silver')
    plots, circ_xs, circ_ys, data=1
    

foreach probe, probes, probe_id do begin
    color = probe_colors[probe_id]
    color = sgcolor('black')
    psym = probe_psyms[probe_id]
    r_sm_var = 'rbsp'+probe+'_r_sm'
    r_sm = get_var_data(r_sm_var, in=time_range, times=times)
    plots, r_sm[*,0], r_sm[*,1], color=color, linestyle=1
    
    msg = 'RBSP-'+strupcase(probe)
    tx = tpos[0]+xchsz*1+xchsz*8*(probe_id)
    ty = tpos[3]-ychsz*1
    xyouts, tx+xchsz*1.5,ty, msg, normal=1, color=color
    plots, tx+xchsz*0.5, ty+ychsz*0.3, normal=1, psym=psym, $
        color=color, symsize=symsize, thick=symthick

    event_r_sm = sinterpol(r_sm, times, event_times)
    foreach event_time, event_times, event_id do begin
        event_color = event_colors[event_id]
        plots, event_r_sm[event_id,0], event_r_sm[event_id,1], $
            color=event_color, psym=psym, symsize=symsize, thick=symthick
        tmp = convert_coord(event_r_sm[event_id,0], event_r_sm[event_id,1], data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.0
        ty = tmp[1]+ychsz*0.5
        if probe eq 'b' then ty = tmp[1]-ychsz*1
        msg = strupcase(probe)+string(event_id+1,format='(I0)')
        xyouts, tx,ty,normal=1, msg, alignment=0.5, color=event_color, charsize=label_size
    endforeach
endforeach


plot, xrange, yrange, $
    xstyle=1, ystyle=1, $
    xtitle=xtitle, ytitle=ytitle, $
    xrange=xrange, yrange=yrange, $
    xticklen=xticklen, yticklen=yticklen, $
    xtickv=xtickv, ytickv=ytickv, $
    xticks=xticks, yticks=yticks, $
    xminor=xminor, yminor=yminor, $
    nodata=1, noerase=1, position=tpos

tx = tpos[0]-xchsz*6
ty = tpos[3]-ychsz*0.7
msg = pan_letters[0]+')'
xyouts, tx, ty, normal=1, msg


sgclose

end