;+
; Plot spacecraft location.
;-

function micro_injection_fig_sc_location_v02, input_event_id, test=test, get_name=get_name, update=update, errmsg=errmsg

    errmsg = ''
    retval = !null
    version = 'v02'
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


    ; Model related vars.
    probe = '4'
    prefix = 'mms'+probe+'_'
    mission_probe = 'mms'+probe

    ; fields.
    ; B field related vars.
    field_time_range = time_range+[-1,1]*30.*60
    b_gsm_var = lets_read_this(func='mms_read_bfield', $
        field_time_range, probe=mission_probe, coord=default_coord, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load B field data ...'
        return, retval
    endif
    e_gsm_var = lets_read_this(func='mms_read_efield', $
        field_time_range, probe=mission_probe, coord=default_coord, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load E field data ...'
        return, retval
    endif

    ; Orbit related vars.
    orbit_time_range = time_range+[-1,1]*30*60
    r_gsm_var = lets_read_this(func='mms_read_orbit', $
        orbit_time_range, probe=mission_probe, coord=default_coord)
    print, 'Loading '+r_gsm_var+' ...'
    mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)
    foreach var, mlat_vars.values() do print, 'Loading '+var+' ...'

    external_models = ['t89','t96','t01','t04s']
    internal_models = ['dipole','igrf']
    hemispheres = ['north','south']

    foreach external_model, external_models do begin
        foreach internal_model, internal_models do begin
            ; The B at sc position.
            suffix = '_'+internal_model+'_'+external_model
            bmod_var = prefix+'bmod_gsm'+suffix
            if tnames(bmod_var) ne '' then del_data, bmod_var
            bmod_var = lets_read_geopack_bfield(var_info=bmod_var, $
                orbit_var=r_gsm_var, time_var=orbit_time_var, $
                internal_model=internal_model, external_model=external_model, save_to=data_file, update=update)
            print, 'Loading '+bmod_var+' ...'
        endforeach
    endforeach
    
    
    ; B model.
    external_model = 't89'
    b0_window = 15.*60
    bmod_var = prefix+'bmod_gsm_igrf_t89'
    b_vars = lets_decompose_bfield(b0_window=b0_window, b_var=b_gsm_var, bmod_var=bmod_var)
    b0_gsm_var = b_vars['b0']
    b_elev_var = lets_calc_vec_elev(b_gsm_var, coord='sm')
    bmod_elev_var = lets_calc_vec_elev(bmod_var, coord='sm', var_info=prefix+'bmod_elev')
    db_elev_var = lets_subtract_vars(b_elev_var, bmod_elev_var, save_to=prefix+'db_elev')
    options, db_elev_var, constant=0, yrange=[-1,1]*90

    b0_var = b_vars['b0']
    b1_var = b_vars['b1']

    ion_vel_var = lets_read_this(func='mms_read_ion_vel', $
        time_range, probe=mission_probe, errmsg=errmsg)

    ; Convert to FAC.
    fac_coord = 'mms_fac'
    options, [r_gsm_var,b0_var], mission='mms'
    q_fac_var = lets_define_fac(r_var=r_gsm_var, b_var=b0_var, fac_coord=fac_coord)
    coord_msgs = [default_coord,fac_coord]
    fac_vars = list()
    foreach var, [b1_var,e_gsm_var,ion_vel_var] do begin
        options, var, mission='mms'
        fac_vars.add, lets_cotran(coord_msgs, input=var, q_var=q_fac_var)
    endforeach
    
    get_data, prefix+'b1_mms_fac', times, b1_fac, limits=lim
    b0_mag = snorm(get_var_data(b0_var, at=times))
    b1_fac[*,0] += b0_mag
    store_data, prefix+'b_mms_fac', times, b1_fac, limits=lim
    

    ; Vexb.    
    b0_mag = snorm(get_var_data(b0_var, times=times))
    e_fac = get_var_data(prefix+'e_mms_fac', at=times)
    ntime = n_elements(times)
    ndim = 3
    b_fac = fltarr(ntime,ndim)
    b_fac[*,0] = b0_mag
    vexb_fac = vec_cross(e_fac,b_fac)
    coef = 1e3/b0_mag^2
    for ii=0,ndim-1 do vexb_fac[*,ii] *= coef
    var = prefix+'vexb_mms_fac'
    store_data, var, times, vexb_fac
    add_setting, var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'V', $
        'unit', 'km/s' )



;---Make the plot.
    prefix = 'mms4_'
    plot_vars = prefix+['b','u','r']+'_'+default_coord
    plot_tr = time_range
    tickinterval = 10*60d
    xrange = [15,-15]
    yrange = [15,-15]
    zrange = [-10,10]
    xrange = [13,-5]
    yrange = [15,-2]
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
    
    abs_xsize = 4

    ypans = abs(total(yrange*[-1,1]))
    pansize = abs([total(xrange*[-1,1]),total(yrange*[-1,1])])
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


;---XY plane.
    tpos = poss[*,0]
    xtickformat = ''
    xticklen = abs_xticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = abs_yticklen/(tpos[2]-tpos[0])/fig_size[0]
    set_axis, xrange=xrange, yrange=yrange, position=tpos, iso=1
    

;---Add Model B field lines.
    probe = '4'
    prefix = 'mms'+probe+'_'
    r_var = prefix+'r_'+default_coord

    angle_range = [0,360]
    angle_step = 30
    angles = make_bins(angle_range, angle_step, inner=1)
    npoint = n_elements(angles)
    ndim = 3
    r_sms = fltarr(npoint,ndim)
    r_coord = get_var_data(r_var, times=times)
    rr = mean(snorm(r_coord))
    r_sms[*,0] = rr*cos(angles)
    r_sms[*,1] = rr*sin(angles)
    r_gsms = cotran_pro(r_sms, model_time, coord_msg=[default_coord,'gsm'])
    flines = lets_get_bfield_lines(model_time, r_gsms, external_model=external_model)

    foreach fline, flines do begin
        f_sm = cotran_pro(fline, model_time, coord_msg=['gsm',default_coord])
        index = where(f_sm[*,2] lt 0, count)
        if count eq 0 then continue
        oplot, f_sm[index,0], f_sm[index,1], color=fline_color, linestyle=2
    endforeach
    tmp = lets_add_earth()
    foreach fline, flines do begin
        f_sm = cotran_pro(fline, model_time, coord_msg=['gsm',default_coord])
        index = where(f_sm[*,2] ge 0, count)
        if count eq 0 then continue
        oplot, f_sm[index,0], f_sm[index,1], color=fline_color
    endforeach

    draw_axis, xrange=xrange, yrange=yrange, noxtitle=0, $
        xtitle=xtitle, ytitle=ytitle, xstep=5, ystep=5, position=tpos

;---magnetopause.
    tts = smkarthm(0,2*!dpi,50,'n')
    rrs = 2
    xxs = rrs*cos(tts)
    yys = rrs*cos(tts)
    ;    xxs = r_coord[*,0]
    ;    yys = r_coord[*,1]
    zzs = fltarr(n_elements(xxs))
    pdyn_var = omni_read_sw_p(time_range)
    test_time = mean(time_range)
    pdyn = get_var_data(pdyn_var, at=test_time)
    mpause_t96, pdyn, xmgnp=xmgnp, ymgnp=ymgnp, zmgnp=zmgnp, $
        xgsm=xxs, ygsm=yys, zgsm=zzs, id=id, distan=distan
    oplot, xmgnp, ymgnp, linestyle=1
    
    probe = '4'
    prefix = 'mms'+probe+'_'
    r_var = prefix+'r_'+default_coord
    f_var = lets_trace_to_equator(orbit_var=r_var, external_model=external_model)
    r_coord = get_var_data(r_var, times=times)
    r_color = sgcolor('misty_rose')
    r_color = sgcolor('red')
    oplot, r_coord[*,0], r_coord[*,1], color=r_color
    f_color = sgcolor('magenta')
    f_coord = get_var_data(f_var, at=times)
    oplot, f_coord[*,0], f_coord[*,1], color=f_color
    
    the_times = make_bins(minmax(times),3600,inner=1)
    foreach the_time, the_times, tid do begin

    ;---SC footpoint.
        tx = interpol(f_coord[*,0], times, the_time)
        ty = interpol(f_coord[*,1], times, the_time)

        if tid mod 3 eq 0 then begin
            plots, tx,ty, psym=psym, symsize=symsize*1.5, color=f_color

            tmp = convert_coord(tx,ty,data=1,to_normal=1)
            dx =-xchsz*0.6
            dy =-ychsz*0.3
            alignment = 1
            if product(the_time-time_double(['2015-09-01/18:00','2015-09-01/20:00'])) le 0 then begin
                dy =-ychsz*1
                alignment = 0.5
            endif
            tx = tmp[0]+dx
            ty = tmp[1]+dy
            msg = time_string(the_time,tformat='hh:mm')+' UT'
        endif else begin
            plots, tx,ty, psym=psym, symsize=symsize, color=f_color
        endelse
        
    ;---SC location.
        tx = interpol(r_coord[*,0], times, the_time)
        ty = interpol(r_coord[*,1], times, the_time)

        if tid mod 3 eq 0 then begin
            plots, tx,ty, psym=psym, symsize=symsize*1.5, color=r_color

            tmp = convert_coord(tx,ty,data=1,to_normal=1)
            dx =-xchsz*0.6
            dy =-ychsz*0.3
            alignment = 1
            if product(the_time-time_double(['2015-09-01/18:00','2015-09-01/20:00'])) le 0 then begin
                dy =-ychsz*1
                alignment = 0.5
            endif
            tx = tmp[0]+dx
            ty = tmp[1]+dy
            msg = time_string(the_time,tformat='hh:mm')+' UT'
            xyouts, tx,ty,msg, normal=1, alignment=alignment, color=r_color
        endif else begin
            plots, tx,ty, psym=psym, symsize=symsize, color=r_color
        endelse
    endforeach



;---Draw velocity vectors along path.
    vel_color = sgcolor('green')
    data_scale = 10         ; km.
    norm_scale = ychsz*5  ; ychsz.
    vel_scale = norm_scale/data_scale
    orbit_time_range = time_double(['2015-09-01/18:00','2015-09-01/21:30'])
    ;orbit_time_range = time_range
    min_u_mag = 75
    the_r_var = f_var
    ion_vel_var = prefix+'vexb_mms_fac'
    ;ion_vel_var = prefix+'u_mms_fac'
    times = make_bins(orbit_time_range, 15, inner=1)
    u_fac = get_var_data(ion_vel_var, at=times)
    r_coord = get_var_data(the_r_var, at=times)
    foreach time, times, tid do begin
        ;if time mod 60 ne 0 then continue
        x0 = r_coord[tid,0]
        y0 = r_coord[tid,1]
        u_west = u_fac[tid,1]
        u_out = u_fac[tid,2]
        if snorm([u_west,u_out]) le min_u_mag then continue    ; only show large vel.
        out_hat = sunitvec([x0,y0,0])
        b_hat = [0,0,1]
        west_hat = vec_cross(out_hat,b_hat)
        x1 = x0+(west_hat[0]*u_west+out_hat[0]*u_out)*vel_scale
        y1 = y0+(west_hat[1]*u_west+out_hat[1]*u_out)*vel_scale
        plots, [x0,x1],[y0,y1], color=vel_color
        plots, x1, y1, color=vel_color;, normal=1
    endforeach

    ; Add scale.
    nn = 15
    tx = tpos[0]+xchsz*10
    ty = tpos[1]+ychsz*1.5
    tmp = convert_coord(tx,ty, normal=1, to_data=1)
    tx = tmp[0]
    ty = tmp[1]
    txs = tx+[0,norm_scale*nn]
    tys = ty+[0,0]
    plots, txs, tys, color=vel_color
    foreach tx,txs do begin
        tmp = convert_coord(tx,tys[0],data=1,to_normal=1)
        plots, tmp[0]+[0,0],tmp[1]+[-1,1]*ychsz*0.1, normal=1
    endforeach
    msg = 'v = '+string(data_scale*nn,format='(I0)')+' km/s'
    tmp = convert_coord(mean(txs),ty, data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]+ychsz*0.5
    xyouts, tx,ty, msg,normal=1, charsize=label_size, alignment=0.5, color=vel_color
    

    
;---label.
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = 'a) XY Plane'
    xyouts, tx,ty,msg, normal=1




    if keyword_set(test) then stop
    sgclose

    return, plot_file
    
end



event_id = '2015_0901_10'
print, micro_injection_fig_sc_location_v02(event_id, test=0, update=1)
end
