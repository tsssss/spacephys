;+
; Plot spacecraft location.
;-

function micro_injection_fig_sc_location_v01, input_event_id, test=test, get_name=get_name, update=update, errmsg=errmsg

    errmsg = ''
    retval = !null
    version = 'v01'
    project= micro_injection_load_project()
    project_id = project.id

    if n_elements(input_event_id) eq 2 then begin
        time_range = time_double(input_event_id)
        event_id = time_string(time_range[0],tformat='YYYY_MMDD_hh')
    endif else begin
        event_id = input_event_id
    endelse
    event = project_get_event(project, id=event_id)
    time_range = event.time_range
    if n_elements(event) eq 0 then message, 'Inconsistency ...'

    if n_elements(plot_dir) eq 0 then plot_dir = event.plot_dir
    base = project_id+'_fig_sc_location_'+event_id+'_'+version+'.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(get_name) then return, plot_file
    print, plot_file
    if keyword_set(update) then file_delete, plot_file, allow_nonexist=1
    if keyword_set(test) then begin
        plot_file = 0
    endif else begin
        if file_test(plot_file) eq 1 then begin
            print, plot_file+' exists, skip ...'
            return, plot_file
        endif
    endelse


;---Load data.
    default_coord = 'sm'
    probes = string(findgen(4)+1,format='(I0)')
    nprobe = n_elements(probes)
    colors = sgcolor(['red','green','blue','purple'])
    labels = strupcase('mms'+probes)
    comps = constant('xyz')
    ncomp = n_elements(comps)
    
    foreach phys_quant, ['orbit'] do begin
        vars = list()
        foreach probe, probes do begin
            vars.add, lets_read(phys_quant, time_range, source=['mms',probe], coord=default_coord)
        endforeach
        
        ; recombine according to component.
        vars = vars.toarray()
        for ii=0,ncomp-1 do begin
            var = 'mms_'+phys_quant+'_'+comps[ii]
            times = get_var_time(vars[0])
            ntime = n_elements(times)
            data = fltarr(ntime,nprobe)
            for jj=0,nprobe-1 do data[*,jj] = (get_var_data(vars[jj],at=times, limits=lim))[*,ii]
            store_data, var, times, data, limits=lim
            options, var, labels=labels, colors=colors
        endfor

        if phys_quant eq 'orbit' then begin
            for ii=0,ncomp-1 do begin
                var = 'mms_'+phys_quant+'_'+comps[ii]
                get_data, var, times, data

                var = 'mms_d'+phys_quant+'_'+comps[ii]
                del_data, var
                data0 = data[*,0]
                for jj=0,nprobe-1 do data[*,jj] -= data0
                ;data[*,0] = !values.f_nan
                data *= constant('re')
                store_data, var, times, data
                options, var, labels=labels, colors=colors, ytitle='(km)', labflag=-1
            endfor
        endif
    endforeach


;---Make the plot.
    prefix = 'mms1_'
    plot_vars = prefix+['b','u','r']+'_'+default_coord
    plot_tr = time_range
    tickinterval = 10*60d
    xrange = [15,-15]
    yrange = [15,-15]
    zrange = [-10,10]
    xrange = [10,-2]
    yrange = [10,-2]
    yrange = [12,-2]
    zrange = [-10,2]
    xtitle = strupcase(default_coord)+' X (Re)'
    ytitle = strupcase(default_coord)+' Y (Re)'
    ztitle = strupcase(default_coord)+' Z (Re)'
    
    ; in Re.
    re = constant('re')
    r_coord = get_var_data('mms1_r_'+default_coord)
    dxrange = minmax(get_var_data('mms_orbit_x'))
    dyrange = minmax(get_var_data('mms_orbit_y'))
    dzrange = minmax(get_var_data('mms_orbit_z'))
    del_step = 1e3  ; km.
    delx = total(minmax(dxrange*[-1,1]))
    dely = total(minmax(dyrange*[-1,1]))
    delz = total(minmax(dzrange*[-1,1]))
    dxrange = reverse(dxrange)
    dyrange = reverse(dyrange)
    
    abs_xsize = 5

    ypans = abs([total(zrange*[-1,1]),total(yrange*[-1,1])])
    pansize = abs([total(xrange*[-1,1]),total(zrange*[-1,1])])
    pansize = pansize/pansize[0]*abs_xsize
    margins = [10,4,2,2]
    poss = panel_pos(ypans=ypans,pansize=pansize, fig_size=fig_size, margins=margins)
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    abs_ticklen = -0.3*ychsz*fig_size[1]
    abs_xticklen = abs_ticklen
    abs_yticklen = abs_ticklen
    
    ; Field lines.
    model_time = mean(time_range)
    snapshot_time = model_time
    fline_color = sgcolor('silver')
    psym = 8
    symsize = 0.5
    tmp = smkarthm(0,2*!dpi,10,'n')
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs, tys, fill=1
    

;---XZ plane.
    tpos = poss[*,0]
    xtickformat = '(A1)'
    xticklen = abs_xticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = abs_yticklen/(tpos[2]-tpos[0])/fig_size[0]
    set_axis, xrange=xrange, yrange=zrange, position=tpos, iso=1
    tmp = lets_add_earth()
    foreach probe, probes[0], pid do begin
        prefix = 'mms'+probe+'_'
        var = prefix+'r_'+default_coord
        r_coord = get_var_data(var, times=times)
        oplot, r_coord[*,0], r_coord[*,2], color=colors[pid]
        
        foreach tid, [0,ntime-1] do begin
            tx = r_coord[tid,0]
            ty = r_coord[tid,2]
            plots, tx,ty, psym=psym, symsize=symsize, color=colors[pid]

            tmp = convert_coord(tx,ty,data=1,to_normal=1)
            dx =-xchsz*0.6
            dy =-ychsz*0.3
            alignment = 1
            tx = tmp[0]+dx
            ty = tmp[1]+dy
            msg = time_string(times[tid],tformat='hh:mm')+' UT'
            xyouts, tx,ty,msg, normal=1, alignment=alignment, color=colors[pid]
        endforeach
    endforeach


    ; Field lines.
    mlts = pseudo_mlt(r_coord)
    the_mlt = mean(mlts)

    xstep = 1
    xs = make_bins([1,xrange[1]], xstep, inner=1)
    npoint = n_elements(xs)
    ndim = 3
    
    zstep = 1
    zs = make_bins([zrange[0],-1], zstep, inner=1)
    npoint = n_elements(zs)
    r_sm2s = fltarr(npoint,ndim)
    r_sm2s[*,0] = mean(r_coord[*,0])
    r_sm2s[*,0] = 4
    ;the_angle = (the_mlt*15)*constant('rad')
    ;r_sm2s[*,0] =-4*cos(the_angle)
    ;r_sm2s[*,1] =-4*sin(the_angle)
    r_sm2s[*,2] = zs
    r_sms = r_sm2s
    
    r_gsms = cotran_pro(r_sms, model_time, coord_msg=['sm','gsm'])
    flines = lets_get_bfield_lines(model_time, r_gsms)
    foreach fline, flines do begin
        f_sm = cotran_pro(fline, model_time, coord_msg=['gsm','sm'])
        ;srotate, f_sm, -(the_angle-!dpi), 2
        oplot, f_sm[*,0], f_sm[*,2], color=fline_color
    endforeach    
    draw_axis, xrange=xrange, yrange=zrange, noxtitle=1, $
        xtitle=xtitle, ytitle=ztitle, xstep=5, ystep=5, position=tpos
    
    
    ; zoom in panel.
;    pansize = [delx,delz]
;    abs_delx = abs_xsize*0.5
;    pansize = pansize/pansize[0]*total(tpos[[0,2]]*[-1,1])*0.3
;    zpos = tpos
;    zpos[0] = tpos[0]+xchsz*4
;    zpos[1] = tpos[1]+ychsz*3
;    zpos[2] = zpos[0]+pansize[0]
;    zpos[3] = zpos[1]+pansize[1]*fig_size[0]/fig_size[1]
;    
;    set_axis, xrange=dxrange, yrange=dzrange, position=zpos, iso=1
;    foreach probe, probes, pid do begin
;        prefix = 'mms'+probe+'_'
;        var = prefix+'r_'+default_coord
;        r_coord = get_var_data(var, times=times)
;        oplot, r_coord[*,0], r_coord[*,2], color=colors[pid]
;    endforeach
;    draw_axis, xrange=dxrange, yrange=dzrange, noxtitle=0, $
;        xtitle='', ytitle='', xstep=0.1, ystep=0.1, position=zpos

    
    color = sgcolor('black')
    tmp = sinterpol(r_coord, times, snapshot_time)
    tx = tmp[0]
    ty = tmp[2]
    plots, tx,ty, psym=psym, symsize=symsize, color=color

    xxs = get_var_data('mms_orbit_x', at=snapshot_time)
    yys = get_var_data('mms_orbit_z', at=snapshot_time)
    xx0 = xxs[0]
    yy0 = yys[0]
    re = constant('re')
    xxs = (xxs-xx0)*re
    yys = (yys-yy0)*re
    dstep = 1e2
    dxrange = minmax(xxs)
    dyrange = minmax(yys)
    dxrange = mean(dxrange)+[-1,1]*1.01*total(dxrange*[-1,1])*0.5
    dyrange = mean(dyrange)+[-1,1]*1.01*total(dyrange*[-1,1])*0.5
    dxrange = reverse(minmax(make_bins(dxrange,dstep)))
    dyrange = minmax(make_bins(dyrange,dstep))
    delx = abs(total(dxrange*[-1,1]))
    dely = abs(total(dyrange*[-1,1]))
    pansize = [delx,dely]
    pansize = pansize/pansize[0]*total(tpos[[0,2]]*[-1,1])*0.3
    zpos = tpos
    zpos[0] = tpos[0]+(tpos[2]-tpos[0])*0.5+xchsz*8
    zpos[1] = tpos[1]+ychsz*4
    zpos[2] = zpos[0]+pansize[0]
    zpos[3] = zpos[1]+pansize[1]*fig_size[0]/fig_size[1]
    
    tmp = sinterpol(r_coord, times, snapshot_time)
    tx = tmp[0]
    ty = tmp[2]
    tmp = convert_coord(tx,ty,data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]
    plots, [tx,zpos[0]],[ty,zpos[3]], normal=1
        
    set_axis, xrange=dxrange, yrange=dyrange, position=zpos, iso=1
    tx = zpos[0]+xchsz*0.5
    ty = zpos[3]-ychsz*1
    msg = 'a-1) '+time_string(snapshot_time,tformat='hh:mm')+' UT'
    xyouts, tx,ty,msg, normal=1
    foreach probe, probes, pid do begin
        plots, xxs[pid],yys[pid], color=colors[pid], psym=psym
        tx = xxs[pid]
        ty = yys[pid]
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]-xchsz*2
        if probe eq '3' then tx = tmp[0]-xchsz*1
        ty = tmp[1]-ychsz*1.1
        msg = strupcase('mms'+probe)
        xyouts, tx,ty,msg, color=colors[pid], normal=1
        
        if probe eq '1' then begin
            pos = sinterpol(r_coord, times, snapshot_time)
            msg = strjoin(strtrim(string(pos[[0,2]],format='(F4.1)'),2),',')
            msg = '('+msg+') Re'
            tx = tmp[0]
            ty = tmp[1]+ychsz*0.5
            xyouts, tx,ty,msg, normal=1, color=colors[pid], alignment=0.5
        endif
    endforeach
    draw_axis, xrange=dxrange, yrange=dyrange, position=zpos, $
        xtitle='SM X (km)', ytitle='SM Z (km)', xstep=dstep, ystep=dstep
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = 'a) XZ Plane'
    xyouts, tx,ty,msg, normal=1


;---XY plane.
    tpos = poss[*,1]
    xtickformat = ''
    xticklen = abs_xticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = abs_yticklen/(tpos[2]-tpos[0])/fig_size[0]
    set_axis, xrange=xrange, yrange=yrange, position=tpos, iso=1
    foreach probe, probes[0], pid do begin
        prefix = 'mms'+probe+'_'
        var = prefix+'r_'+default_coord
        r_coord = get_var_data(var, times=times)
        oplot, r_coord[*,0], r_coord[*,1], color=colors[pid]
        
        foreach tid, [0,ntime-1] do begin
            tx = r_coord[tid,0]
            ty = r_coord[tid,1]
            plots, tx,ty, psym=psym, symsize=symsize, color=colors[pid]

            tmp = convert_coord(tx,ty,data=1,to_normal=1)
            dx =-xchsz*0.6
            dy =-ychsz*0.3
            alignment = 1
            tx = tmp[0]+dx
            ty = tmp[1]+dy
            msg = time_string(times[tid],tformat='hh:mm')+' UT'
            xyouts, tx,ty,msg, normal=1, alignment=alignment, color=colors[pid]
        endforeach
    endforeach


    angle_range = [0,360]
    angle_step = 30
    angles = make_bins(angle_range, angle_step, inner=1)
    npoint = n_elements(angles)
    ndim = 3
    r_sms = fltarr(npoint,ndim)
    rr = mean(snorm(r_coord))
    r_sms[*,0] = rr*cos(angles)
    r_sms[*,1] = rr*sin(angles)
    r_gsms = cotran_pro(r_sms, model_time, coord_msg=['sm','gsm'])
    flines = lets_get_bfield_lines(model_time, r_gsms)

    foreach fline, flines do begin
        f_sm = cotran_pro(fline, model_time, coord_msg=['gsm','sm'])
        index = where(f_sm[*,2] lt 0, count)
        if count eq 0 then continue
        oplot, f_sm[index,0], f_sm[index,1], color=fline_color, linestyle=2
    endforeach
    tmp = lets_add_earth()
    foreach fline, flines do begin
        f_sm = cotran_pro(fline, model_time, coord_msg=['gsm','sm'])
        index = where(f_sm[*,2] ge 0, count)
        if count eq 0 then continue
        oplot, f_sm[index,0], f_sm[index,1], color=fline_color
    endforeach
    
    draw_axis, xrange=xrange, yrange=yrange, noxtitle=0, $
        xtitle=xtitle, ytitle=ytitle, xstep=5, ystep=5, position=tpos


;    ; zoom in panel.
    color = sgcolor('black')
    tmp = sinterpol(r_coord, times, snapshot_time)
    tx = tmp[0]
    ty = tmp[1]
    plots, tx,ty, psym=psym, symsize=symsize, color=color
    
    xxs = get_var_data('mms_orbit_x', at=snapshot_time)
    yys = get_var_data('mms_orbit_y', at=snapshot_time)
    xx0 = xxs[0]
    yy0 = yys[0]
    re = constant('re')
    xxs = (xxs-xx0)*re
    yys = (yys-yy0)*re
    dstep = 1e2
    dxrange = minmax(xxs)
    dyrange = minmax(yys)
    dxrange = mean(dxrange)+[-1,1]*1.01*total(dxrange*[-1,1])*0.5
    dyrange = mean(dyrange)+[-1,1]*1.01*total(dyrange*[-1,1])*0.5
    dxrange = reverse(minmax(make_bins(dxrange,dstep)))
    dyrange = minmax(make_bins(dyrange,dstep))
    delx = abs(total(dxrange*[-1,1]))
    dely = abs(total(dyrange*[-1,1]))
    pansize = [delx,dely]
    pansize = pansize/pansize[0]*total(tpos[[0,2]]*[-1,1])*0.3
    zpos = tpos
    zpos[0] = tpos[0]+(tpos[2]-tpos[0])*0.1+xchsz*8
    zpos[1] = tpos[1]+(tpos[3]-tpos[1])*0.2+ychsz*4
    zpos[2] = zpos[0]+pansize[0]
    zpos[3] = zpos[1]+pansize[1]*fig_size[0]/fig_size[1]
    
    tmp = sinterpol(r_coord, times, snapshot_time)
    tx = tmp[0]
    ty = tmp[1]
    tmp = convert_coord(tx,ty,data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]
    plots, [tx,zpos[2]],[ty,zpos[1]], normal=1
    
    set_axis, xrange=dxrange, yrange=dyrange, position=zpos, iso=1
    tx = zpos[0]+xchsz*0.5
    ty = zpos[3]-ychsz*1
    msg = 'b-1) '+time_string(snapshot_time,tformat='hh:mm')+' UT'
    xyouts, tx,ty,msg, normal=1
    foreach probe, probes, pid do begin
        plots, xxs[pid],yys[pid], color=colors[pid], psym=psym
        tx = xxs[pid]
        ty = yys[pid]
        tmp = convert_coord(tx,ty,data=1,to_normal=1)
        tx = tmp[0]-xchsz*2
        ;if probe eq '3' then tx = tmp[0]-xchsz*1
        ty = tmp[1]-ychsz*1.1
        msg = strupcase('mms'+probe)
        xyouts, tx,ty,msg, color=colors[pid], normal=1
    
        if probe eq '1' then begin
            pos = sinterpol(r_coord, times, snapshot_time)
            msg = strjoin(strtrim(string(pos[[0,1]],format='(F4.1)'),2),',')
            msg = '('+msg+') Re'
            tx = tmp[0]
            ty = tmp[1]+ychsz*0.5
            xyouts, tx,ty,msg, normal=1, color=colors[pid], alignment=0.5
        endif
    endforeach
    draw_axis, xrange=dxrange, yrange=dyrange, position=zpos, $
        xtitle='SM X (km)', ytitle='SM Y (km)', xstep=dstep, ystep=dstep
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = 'b) XY Plane'
    xyouts, tx,ty,msg, normal=1
    

    if keyword_set(test) then stop
    sgclose

    return, plot_file
    
end


event_id = '2015_0901_18'
print, micro_injection_fig_sc_location_v01(event_id, test=0)
stop

event_id = '2015_0901_11'
print, micro_injection_fig_sc_location_v01(event_id, test=0)
end