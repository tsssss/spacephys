;+
; Plot the absolute value of the angle between the 2D velocity and x_sm versus MLT.
;-

pro azim_prop_plot_angle_vs_mlt, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events

    test = 0

;---Plot settings.
    margins = [6,4,7,1]
    fig_xsize = 3
    fig_ysize = 2
    yticklen = -0.02
    half_ychsz = 0.35
    full_ychsz = 0.7
    line_skip = 0.35
    label_size = 0.7
    xpad = 10
    ypad = 0.2
    mlt_panel = 2
    str_phi = '!9'+string(102b)+'!X'
    str_delta = '!9'+string(68b)+'!X'
    label_xpos = 5

    deg = project.constant.deg
    rad = project.constant.rad
    re = project.constant.re

    vel_mlt = hash()
    all_mlt = []
    all_angle = []
    all_vmag = []
    all_rxy = []
    v2d_angle = project.v2d_angle
    v2d_dt = project.v2d_dt
    foreach event_id, events.keys() do begin
        combos = project.events[event_id].probe_combos
        event_vmag = []
        omega_2ds = []
        event_rxy = []
        foreach combo, combos do begin
            key = strjoin(combo[sort(combo)],'_')
            combo_info = (project.events[event_id])[key]
            v_mag = combo_info.v_mag
            v_hat = combo_info.v_hat
            rsm_center = combo_info.rsm_center

            ; Filter by time difference.
            times = combo_info.times
            time_diffs = abs([times[0]-times[1],times[1]-times[2],times[2]-times[0]])
            if min(time_diffs) lt v2d_dt then continue

            ; Filter by geometry.
            triangle_angles = combo_info.triangle_angles
            if max(triangle_angles) gt 180-v2d_angle then continue
            if min(triangle_angles) lt v2d_angle then continue


            if ~vel_mlt.haskey(event_id) then vel_mlt[event_id] = dictionary()
            ; make the angle in [-180,180].
            angle_from_xsm = atan(v_hat[1],v_hat[0])*deg
            if angle_from_xsm gt 180 then angle_from_xsm -= 360
            if angle_from_xsm lt -180 then angle_from_xsm += 360

            event_info = project.events[event_id]
            ref_times = dblarr(n_elements(combo))
            foreach probe, combo, ii do ref_times[ii] = event_info[probe].ref_time
            ut_center = mean(ref_times)
            rmag = cotran([rsm_center,0], ut_center, 'sm2mag')
            mlon = atan(rmag[1],rmag[0])*deg
            mlt = mlon2mlt(mlon, ut_center)
            (vel_mlt[event_id])[strjoin(combo,'_')] = {mlt:mlt,theta:angle_from_xsm}
            all_mlt = [all_mlt,mlt]
            all_angle = [all_angle,angle_from_xsm]

            angle_from_azimuth = combo_info.angle_from_azimuth
            omega_2d = v_mag/snorm(rsm_center)/re*deg*cos(angle_from_azimuth*rad)*60
            omega_2ds = [omega_2ds, omega_2d]
            event_vmag = [event_vmag,v_mag]
            event_rxy = [event_rxy,snorm(rsm_center)]
        endforeach
;        if keyword_set(test) then begin
;            if n_elements(event_rxy) lt 2 then continue
;            tpos = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz)
;            plot, event_rxy, event_vmag, psym=1, symsize=0.5, /ynozero, position=tpos
;            tx = tpos[0]
;            ty = tpos[3]+ychsz*0.5
;            xyouts, tx,ty,/normal, event_id
;            stop
;        endif
    endforeach

    ;print, 'omega 2D = '+sgnum2str(mean(omega_2ds),ndec=3)+'+/-'+$
    ;    sgnum2str(stddev(omega_2ds),ndec=3)+' deg/min'


;---Add plot.
    if keyword_set(test) then begin
        file = test
        magnify = 2
    endif else begin
        file = join_path([project.plot_dir,'vel_analysis','fig_mlt_vs_angle.pdf'])
        magnify = 1
    endelse

    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify
    tpos = sgcalcpos(1, margins=margins, xchsz=xchsz,ychsz=ychsz)

    xrange = [-1,1]*9
    xstep = 3
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = 3
    xticklen = -0.02
    xtitle = 'MLT (hr)'

    yrange = xrange*10
    ystep = 45
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 3
    yticklen = -0.02
    ytitle = str_phi+' (deg)'

    psym = 1

    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtitle=xtitle, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        position=tpos, /nodata, /noerase


;---3 schemes.
    schemes = hash($
        'azim.', dictionary($
            'pre_mlt', [-12.,0], $
            'pre_phi', [-90.,90], $
            'post_mlt', [0.,12], $
            'post_phi', [-90.,90], $
            'color', sgcolor('blue')), $
        'radial', dictionary($
            'pre_mlt', [-12.,0], $
            'pre_phi', [-180.,0], $
            'post_mlt', [0.,12], $
            'post_phi', [0.,180], $
            'color', sgcolor('green')), $
        'sunward', dictionary($
            'pre_mlt', [-12.,0], $
            'pre_phi', [0.,0], $
            'post_mlt', [0.,12], $
            'post_phi', [0.,0], $
            'color', sgcolor('red')))
        
    linestyle = 1
    foreach scheme, schemes do begin
        color = scheme.color
        foreach key, ['pre','post'] do begin
            txs = scheme[key+'_mlt']
            tys = scheme[key+'_phi']
            oplot, txs,tys, linestyle=linestyle, color=color
        endforeach
    endforeach
    
;---Plot all data.
    txs = all_mlt
    tys = all_angle
    index = where(txs le 0)
    oplot, txs[index], tys[index], psym=psym, symsize=label_size*0.5
    index = where(txs gt 0)
    oplot, txs[index], tys[index], psym=psym, symsize=label_size*0.5
    
    
;---Calculate the coef.
    names = schemes.keys()
    errors = list()
    foreach name, names do begin
        scheme = schemes[name]
        color = scheme.color
        xx = all_mlt
        yy = all_angle
        yfit = all_angle
        foreach key, ['pre','post'] do begin
            tx = scheme[key+'_mlt']
            ty = scheme[key+'_phi']
            index = where_pro(xx, tx)
            yfit[index] = interpol(ty,tx, xx[index])
        endforeach
        
        error = total((yy-yfit)^2)/total((yy-mean(yy))^2)
        schemes[name].error = error
        errors.add, error
    endforeach
    errors = errors.toarray()
    

;---Print to plot.
    tx = tpos[2]+xchsz*0.5
    ty = tpos[3]-ychsz*1*label_size
    color = sgcolor('black')
    xyouts, tx,ty,/normal, 'Model', color=color, charsize=label_size
    tx = tpos[2]+xchsz*0.5+xchsz*4
    xyouts, tx,ty,/normal, 'Err', color=color, charsize=label_size
    
    
    index = sort(errors)
    names = names[index]
    foreach name, names, ii do begin
        tx = tpos[2]+xchsz*0.5
        ty = tpos[3]-ychsz*(ii+2)*label_size
        color = schemes[name].color
        error = schemes[name].error
        xyouts, tx,ty,/normal, name, color=color, charsize=label_size
        
        tx = tpos[2]+xchsz*0.5+xchsz*4
        xyouts, tx,ty,/normal, sgnum2str(error,ndec=1), color=color, charsize=label_size
    endforeach

    if keyword_set(test) then stop
    sgclose

end
