;+
; Analysis events using 3-spacecraft timing, then analysis
;-

test = 1
project = azim_df_load_project()

events = azim_df_filter_load_candidate(project=project)
df_count = 0
foreach event, events do foreach df, event.df_list do df_count += 1
events = azim_df_filter_triads(events, project=project)
events = azim_df_filter_duration(events, project=project)
events = azim_df_filter_v2d(events, project=project)


;---Plot.
rsm_xrange = [10,-30]
rsm_yrange = [1,-1]*20
xsize = abs(total(rsm_xrange*[-1,1]))
ysize = abs(total(rsm_yrange*[-1,1]))
fig_xsize = 5
margins = [8.,5,2,2]
sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
sgclose, /wdelete
xpanel_size = fig_xsize-total(margins[[0,2]])*abs_xchsz
ypanel_size = xpanel_size/xsize*ysize
fig_ysize = ypanel_size+total(margins[[1,3]])*abs_ychsz

tpos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
tpos = [tpos[0:1],1-tpos[2:3]]


file = 1
magnify = keyword_set(test)? 1:1
sgopen, file, xsize=fig_xsize, ysize=fig_ysize, /inch, magnify=magnify, xchsz=xchsz, ychsz=ychsz


timing_key = 'df_time'
azim_event_count = 0
good_event_count = 0

foreach event, events do begin
    timing_info = event.timing_info[timing_key]

    ; Get the good triads.
    triad_timing = timing_info.triad_timing
    good_triads = azim_df_filter_check_good_triads(triad_timing)


    good_event_count += 1


    ; Exclude cases that is inconsistent with azimuthal.
    phi_angles = list()
    if_away_from_midnight = 0
    foreach probe_combo, good_triads do begin
        combo_info = triad_timing[probe_combo]
        v_hat = combo_info.v_hat
        rsm_center = combo_info.center_rsm
        phi_angles.add, sang(rsm_center, v_hat, /deg)
        if_away_from_midnight = rsm_center[1]*v_hat[1] ge 0
    endforeach
    phi_angles = phi_angles.toarray()
    angle_range = 90+[-1,1]*30
    if product(mean(phi_angles)-angle_range) ge 0 then continue
    if if_away_from_midnight eq 0 then continue

    azim_event_count += 1


;---Plot settings.
    xticklen_chsz = -0.3
    yticklen_chsz = -0.3
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    label_size = constant('label_size')

    xrange = rsm_xrange
    xstep = 10
    xminor = 5
    xtickv = make_bins(xrange,xstep,/inner)
    xticks = n_elements(xtickv)-1
    xlog = 0
    xtitle = 'SM X (Re)'

    yrange = rsm_yrange
    ystep = 10
    yminor = 5
    ytickv = make_bins(yrange,ystep,/inner)
    yticks = n_elements(ytickv)-1
    ylog = 0
    ytitle = 'SM Y (Re)'

    ; Set up coord.
    plot, xrange,yrange, /nodata, /noerase, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, ytitle=ytitle, yticklen=yticklen, $
        position=tpos

    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    plots, circle_x, circle_y

    ; Add ROI and dash lines.
    pdyn = 10.
    rxy_range = [4,30]
    mlt_range = [-1,1]*9
    rad = constant('rad')
    deg = constant('deg')
    angle_range = (mlt_range*15+180)*rad
    roi_color = sgcolor('gray')
    roi_thick = keyword_set(test)? 2: 8

    magn_test = fltarr(nangle,3)
    magn_test[*,0] = circle_x*30
    magn_test[*,1] = circle_y*30
    tmp = check_if_in_magn(magn_test, magn_pos=magn_pos, dynamic_pressure=pdyn)
    magn_pos = magn_pos[*,[0,2]]
    magn_rxy = snorm(magn_pos)
    magn_angle = atan(magn_pos[*,1],magn_pos[*,0])
    p_mlt = sinterpol(magn_pos, magn_angle, min(angle_range))
    p_rxy = sinterpol(magn_pos, magn_rxy, max(rxy_range))

    ; Magnetopause dash line.
    foreach sector_flag, [-1,1] do begin
        ; -1: pre-midn; 1: post-midn.
        txs = magn_pos[*,0]
        tys = magn_pos[*,1]*sector_flag
        oplot, txs,tys, linestyle=1
    endforeach

    ; MLT dash line.
    rs = [1,max(rxy_range)]
    foreach t, [angle_range,!dpi] do begin
        oplot, cos(t)*rs, sin(t)*rs, linestyle=1
    endforeach

    ; Rxy dash line.
    foreach r, make_bins(rxy_range,10,/inner) do begin
        oplot, circle_x*r, circle_y*r, linestyle=1
    endforeach

    ; Magnetopause.
    index = where(magn_rxy le max(rxy_range) and magn_angle ge min(angle_range))
    magn_pos = magn_pos[index,*]
    magn_xs = [magn_pos[*,0],p_mlt[0],p_rxy[0]]
    magn_ys = [magn_pos[*,1],p_mlt[1],p_rxy[1]]
    index = sort(magn_ys)
    magn_xs = magn_xs[index]
    magn_ys = magn_ys[index]

    foreach sector_flag, [-1,1] do begin
        ; -1: pre-midn; 1: post-midn.
        txs = magn_xs
        tys = magn_ys*sector_flag
        oplot, txs,tys, color=roi_color, thick=roi_thick
    endforeach
    ; MLT.
    mlt_rxy = snorm(p_mlt)
    rs = [min(rxy_range),mlt_rxy]
    foreach t, angle_range do begin
        oplot, cos(t)*rs, sin(t)*rs, color=roi_color, thick=roi_thick
    endforeach
    ; Rxy.
    rxy_angle = atan(p_rxy[1],p_rxy[0])
    foreach r, rxy_range do begin
        if r eq min(rxy_range) then begin
            ts = smkarthm(angle_range[0],angle_range[1],nangle, 'n')
        endif else begin
            ts = smkarthm(rxy_angle,2*!dpi-rxy_angle,nangle, 'n')
        endelse
        oplot, cos(ts)*r, sin(ts)*r, color=roi_color, thick=roi_thick
    endforeach

    ; Add label.
    label_x = tpos[0]+xchsz*1
    label_y = tpos[3]-ychsz*1
    msg = 'Region of interest (ROI)!C(1) Magnetopause @10 nPa;!C(2) Rxy [4,30] Re;!C(3) MLT [-9,9] hr'
    xyouts, label_x,label_y,/normal, msg


    ; Add arrow.
    timing_type = 'df_time'
    triad_angle_range = [15.,165]
    triad_tdiff_range = [70.,1e10]
    nvertex = 3
    dis_scale = 1.5
    vel_scale = 20
    hsize = keyword_set(test)? 6: 100
    arrow_solid = 0
    fill_alpha = 0.5
    fill_color = sgcolor('white_smoke')
    fill_color2 = sgcolor('light_cyan')
    fill_index = [0,1,2,0]
    fill_psym = 8
    fill_symsize = 0.5
    tmp = 20
    tmp = findgen(tmp)/(tmp-1)*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    usersym, xs, ys, /fill

    ; Add scale.
    scale_x = tpos[2]-xchsz*2
    scale_y = tpos[3]-ychsz*1
    tx = scale_x
    ty = scale_y-ychsz*0.5
    tmp = convert_coord(tx,ty, /normal, /to_data)
    tmp = convert_coord(tmp[0]-dis_scale,tmp[1], /data, /to_normal)
    dx = abs(tmp[0]-tx)
    xs = tx+[-1,0]*dx
    plots, xs, ty+[0,0], /normal
    foreach tx, xs do plots, tx+[0,0], ty+[-1,1]*ychsz*0.1, /normal
    tx = mean(xs)
    ty = scale_y
    xyouts, tx,ty,/normal, alignment=0.5, sgnum2str(vel_scale)+' km/s'


    ; The edge.
    foreach probe_combo, good_triads do begin
        combo_info = triad_timing[probe_combo]
        rsms = combo_info.rsms
        xs = rsms[fill_index,0]
        ys = rsms[fill_index,1]
        ;polyfill, xs, ys, transparent=fill_alpha, color=fill_color
        plots, xs, ys, color=fill_color
    endforeach

    ; The vertex.
    foreach probe_combo, good_triads do begin
        combo_info = triad_timing[probe_combo]
        rsms = combo_info.rsms
        xs = rsms[fill_index,0]
        ys = rsms[fill_index,1]
        ;polyfill, xs, ys, transparent=fill_alpha, color=fill_color
        plots, xs, ys, color=fill_color2, psym=fill_psym, symsize=fill_symsize
    endforeach

;    print, time_string(event.time_range)
;    stop

    ; Add arrow
    foreach probe_combo, good_triads do begin
        combo_info = triad_timing[probe_combo]
        rsm_center = combo_info.center_rsm
        v_mag = combo_info.v_mag
        v_hat = combo_info.v_hat

        x0 = rsm_center[0]
        y0 = rsm_center[1]
        scale = dis_scale/vel_scale*v_mag
        x1 = x0+v_hat[0]*scale
        y1 = y0+v_hat[1]*scale
        arrow, x0,y0,x1,y1,/data, solid=arrow_solid, hsize=hsize
    endforeach


endforeach

print, 'Events with consistent propagation: ', good_event_count
print, 'Azim events: ', azim_event_count


if keyword_set(test) then stop
sgclose


end
