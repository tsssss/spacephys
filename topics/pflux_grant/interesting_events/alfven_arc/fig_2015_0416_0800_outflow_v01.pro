;+
; Test tracing results assuming ion conics and ion beam.
;-

function get_mass0, species

    case species of
        'e': mass0 = 1d/1836
        'p': mass0 = 1d
        'o': mass0 = 16d
        'he': mass0 = 4d
    endcase
    mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.
    return, mass0
    
end

function fig_2015_0416_0800_outflow_v01, event_info=event_info

test = 0

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
    base_name = 'fig_'+id+'_outflow_'+version+'.pdf'

;---RBSP settings.
    rbsp_info = event_info.rbsp.rbspa
    time_range = event_info.time_range
    time_range = time_double(['2015-04-16/08:00','2015-04-16/08:30'])
    prefix = rbsp_info['prefix']
    probe = rbsp_info['probe']


    
    rad = constant('rad')
    deg = constant('deg')
    re = constant('re')
    bar_thick = (keyword_set(test))? 2: 6
    label_size = 0.8
    symsize = 0.5

    ct_oxygen = 64
    ct_proton = 63
    zrange_proton = [1e3,1e6]
    zrange_oxygen = [1e3,1e6]
    
    o_pa_spec_var = rbsp_read_pa_spec(time_range, species='o', probe=probe, update=0, energy_range=[100,5e4])
    vars = o_pa_spec_var
    options, vars, 'color_table', ct_oxygen
    options, vars, 'zrange', zrange_oxygen
    options, vars, 'yrange', [0,180]
    options, vars, 'ytitle', 'PA!C(deg)'
    foreach var, vars do begin
        get_data, var, times, data, val
        index = where(data eq 0, count)
        if count ne 0 then begin
            data[index] = 0.01
            store_data, var, times, data, val
        endif
    endforeach
    
    p_pa_spec_var = rbsp_read_pa_spec(time_range, species='p', probe=probe, update=0, energy_range=[100,5e4])
    vars = p_pa_spec_var
    options, vars, 'color_table', ct_proton
    options, vars, 'zrange', zrange_proton
    options, vars, 'yrange', [0,180]
    options, vars, 'ytitle', 'PA!C(deg)'
    foreach var, vars do begin
        get_data, var, times, data, val
        index = where(data eq 0, count)
        if count ne 0 then begin
            data[index] = 0.01
            store_data, var, times, data, val
        endif
    endforeach

    o_spec_var = prefix+'o_en_spec_'+['anti','perp','para']
    vars = o_spec_var
    options, vars, 'color_table', ct_oxygen
    options, vars, 'zrange', zrange_oxygen
    options, vars, 'yrange', [10,5e4]
    foreach var, vars do begin
        get_data, var, times, data, val
        index = where(data eq 0, count)
        if count ne 0 then begin
            data[index] = 0.01
            store_data, var, times, data, val
        endif
    endforeach

    
    p_spec_var = prefix+'p_en_spec_'+['anti','perp','para']
    vars = p_spec_var
    options, vars, 'color_table', ct_proton
    options, vars, 'zrange', zrange_proton
    options, vars, 'yrange', [10,5e4]
    foreach var, vars do begin
        get_data, var, times, data, val
        index = where(data eq 0, count)
        if count ne 0 then begin
            data[index] = 0.01
            store_data, var, times, data, val
        endif
    endforeach
    
    var = prefix+'o_en_spec_para'
    log_ytickv = [2,3,4]
    ytickv = 10d^log_ytickv
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    yticks = n_elements(ytickv)-1
    yminor = 9
    add_setting, var, dictionary($
        'ytitle', 'Energy!C(eV)', $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'ytickname', ytickn )

        
    
    plot_vars = [prefix+'pfdot0_fac_spinfit_phasef_map',$
        prefix+'o_eflux_map', $
        o_pa_spec_var, prefix+'o_en_spec_para']
    nvar = n_elements(plot_vars)
    fig_letters = letters(nvar+1)
    fig_labels = fig_letters[0:nvar-1]+') '+['S','F!DO!N','O+ PA','O+ EN']
    fig_size = [5,5]
    plot_dir = event_info.plot_dir
    if n_elements(plot_file) eq 0 then plot_file = join_path([plot_dir,base_name])
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    margins = [12,4,10,1]
    all_poss = sgcalcpos(nvar+1, ypad=[fltarr(nvar-1)+0.4,3], margins=margins, xpad=total(margins[[0,2]]))
    ;poss = reform(all_poss[*,0,0:nvar-1])
    ;tplot, plot_vars, trange=time_range, position=poss

    uniform_ticklen = -ychsz*0.15*fig_size[0]


    
    
    
;---Prepare the fline_info for tracing.
    model_time = time_double('2015-04-16/08:10')
    external_model = 't04s'
    internal_model = 'dipole'
    r_gsm_var = prefix+'r_gsm'
    r_gsm = get_var_data(r_gsm_var, at=model_time)
    par_var = external_model+'_par'
    if check_if_update(par_var, time_range) then par_var = geopack_read_par(time_range, model=external_model)
    par = get_var_data(par_var, at=model_time)

    time = model_time
    model = external_model
    fline_info = dictionary()
    foreach trace_dir, [-1,1] do begin
        msg = (trace_dir eq -1)? 'north': 'south'
    
        ps = geopack_recalc(time)
        xp = r_gsm[0]
        yp = r_gsm[1]
        zp = r_gsm[2]
    
        geopack_trace, xp,yp,zp, trace_dir, par, $
            xf,yf,zf, r0=r0, refine=1, ionosphere=1, $
            t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm, igrf=igrf, $
            fline=fline
    
        ntrace = n_elements(fline[*,0])
        ndim = 3
        bline = fltarr(ntrace,ndim)
        for ii=0,ntrace-1 do begin
            rx = fline[ii,0]
            ry = fline[ii,1]
            rz = fline[ii,2]
    
            if keyword_set(igrf) then begin
                geopack_igrf_gsm, rx,ry,rz, bx,by,bz
            endif else begin
                geopack_dip, rx,ry,rz, bx,by,bz
            endelse
    
            if model eq 'dip' or model eq 'dipole' or model eq 'igrf' then begin
                dbx = -1
                dby = -1
                dbz = -1
            endif else begin
                routine = 'geopack_'+model
                if model eq 't04s' then routine = 'geopack_ts04'
                call_procedure, routine, par, rx,ry,rz, dbx,dby,dbz
            endelse
    
            bline[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endfor
    
        fline_info[msg] = dictionary($
            'fline', fline, $
            'bline', bline )
    endforeach
    
    
    
;---Find the best time and distance.
    test_diss = smkarthm(1d,4,0.1,'dx')
    test_times = smkarthm(time_double('2015-04-16/08:07:00'),time_double('2015-04-16/08:11:00'), 10, 'dx')
    ntest_dis = n_elements(test_diss)
    ntest_time = n_elements(test_times)
    beam_corr = fltarr(ntest_time,ntest_dis)
    
    test_var = prefix+'o_en_spec_para'
    get_data, test_var, times, fluxs, energys, limits=lim
    species = 'o'
    mass0 = get_mass0(species)
    test_energy_range = [1000,30000]
    energy_index = where_pro(energys[0,*], '[]', test_energy_range, count=ntest_energy)
    norm_fluxs = alog10(fluxs)
    foreach time, times, time_id do begin
        norm_fluxs[time_id,*] /= stddev(norm_fluxs[time_id,*])   ; to normalize the flux b/c flux at low energy dominates.
    endforeach
    energys = reform(energys[0,*])

    target_time_range = time_double(['2015-04-16/08:08:30','2015-04-16/08:20:30'])
    target_time_range = time_double(['2015-04-16/08:09:00','2015-04-16/08:30:30'])
    time_indexs = where_pro(times, '[]', target_time_range)
    target_times = times[time_indexs]
    ntarget = n_elements(target_times)
    target_energys = fltarr(ntarget)
    target_energy_indexs = where_pro(energys, '[]', [100,30000])    ; used to get the target energys.
    the_energys = energys[target_energy_indexs]
    target_fluxs = fltarr(ntarget)
    foreach time_index, time_indexs, target_id do begin
        the_fluxs = reform(fluxs[time_index,target_energy_indexs])
        target_fluxs[target_id] = max(the_fluxs, the_energy_index)
        target_energys[target_id] = the_energys[the_energy_index]
    endforeach
    target_energys[0] = 2.5e4   ; this one somehow doesn't work.

    
    ; Mapping coef.
    external_model = 't89'
    internal_model = 'dipole'
    r_var = prefix+'r_gsm'
    bmod_var = geopack_read_bfield(r_var=r_var, models=external_model, igrf=0, suffix='_'+internal_model, t89_par=2, coord='gsm')
    bfmod_var = prefix+'bf_gsm_t89_dipole'
    if check_if_update(bfmod_var) then begin
        vinfo = geopack_trace_to_ionosphere(r_var, models=external_model, $
            igrf=0, south=1, north=0, refine=1, suffix='_'+internal_model)
    endif
    cmap = snorm(get_var_data(bfmod_var, times=uts))/snorm(get_var_data(bmod_var, at=uts))
    cmap_var = prefix+'cmap'
    store_data, cmap_var, uts, cmap
    
    ; flux unit is '#/cm!E2!N-s-sr-keV'
    theta_pa = 30d
    sr = !dpi*sin(theta_pa*constant('rad'))^2
    nenergy_bin = 3     ; # of energy_bins around
    eflux = target_fluxs*target_energys^2*1e-3*sr*nenergy_bin  ; in eV/cm^2-s
    eflux *= 1.6e-19*1e4*1e3  ; convert from eV/cm^2-s to mW/m^2.
    cmap = get_var_data(cmap_var, at=target_times)
    var = prefix+'o_eflux_map'
    store_data, var, target_times, eflux*cmap
    yrange = [0,50]
    ystep = 20
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 2
    add_setting, var, smart=1, dictionary($
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'display_type', 'scalar', $
        'short_name', ' ', $
        'unit', 'mW/m!U2!N' )
    
    var = prefix+'pfdot0_fac_spinfit_phasef_map'
    yrange = [-60,30]
    ytickv = [-60,-20,20]
    yticks = n_elements(ytickv)-1
    yminor = 2
    fac_labels = [tex2str('perp')+','+['north','west'],tex2str('parallel')]
    add_setting, var, dictionary($
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'constant', 0, $
        'ytitle', '(mW/m!U2!N)', $
        'labels', 'S!D'+fac_labels+'!N' )
    
    ; Trace backwards in time and along field line.
    flines = fline_info.south.fline
    blines = fline_info.south.bline
    bmags = snorm(blines)
    fline_diss = snorm(flines)
    nfline = n_elements(fline_diss)
    dtimes = times[time_indexs+1]-target_times
    max_target_time = time_double('2015-04-16/08:13')
    
    
    ; Beam.
    trace_times = dblarr(ntarget,nfline)
    foreach fline_dis, fline_diss, fline_id do begin
        if fline_id eq 0 then begin
            trace_times[*,fline_id] = target_times
            continue
        endif

        vmags = sqrt(2*target_energys/mass0)*1e-3 ; in km/s.
        vpara_steps = vmags

        r_step = reform(flines[fline_id,*]-flines[fline_id-1,*])
        dr_step = snorm(r_step)*constant('re')      ; in km.
        dt_steps = dr_step/vpara_steps
        trace_times[*,fline_id] = trace_times[*,fline_id-1]-dt_steps
    endforeach

    wanted_index = where(target_times lt max_target_time)
    trace_dtimes = dblarr(nfline)
    for ii=0,nfline-1 do begin
        trace_dtimes[ii] = stddev(trace_times[wanted_index,ii])
    endfor
    err = min(trace_dtimes, abs=1, trace_index)
    trace_time = mean(trace_times[wanted_index,trace_index])
    trace_dis = fline_diss[trace_index]
    trace_tr = trace_time+[-1,1]*err
    bar_color = sgcolor('red')
    print, time_string(trace_time)+'+/-'+string(err,format='(F6.1)')+' sec'
    
    
    
    var = prefix+'o_eflux_map'
    get_data, var, times, efluxs, limits=lim
    new_times = trace_times[*,trace_index]
    time_step = 22d
    uts = make_bins(time_range, time_step, inner=1)
    ntime = n_elements(uts)
    new_efluxs = fltarr(ntime);+!values.f_nan
    for ii=1,ntime-1 do begin
        index = where_pro(new_times, '[)', uts[ii-1:ii], count=count)
        if count eq 0 then continue
        new_efluxs[ii-1] = total(efluxs[index])
    endfor
    old_efluxs = fltarr(ntime);+!values.f_nan
    for ii=1,ntime-1 do begin
        index = where_pro(times, '[)', uts[ii-1:ii], count=count)
        if count eq 0 then continue
        old_efluxs[ii-1] = total(efluxs[index])
    endfor
    store_data, var, uts, [[old_efluxs],[new_efluxs]]
    options, var, 'labels', ['O+ eflux','Traced']
    options, var, 'colors', sgcolor(['black','orange'])
;    pid = where(plot_vars eq var)
;    tpos = poss[*,pid]
;    xrange = time_range
;    yrange = get_setting(var, 'yrange')
;    plot, xrange, yrange, $
;        xstyle=5, xrange=xrange, $
;        ystyle=5, yrange=yrange, $
;        nodata=1, noerase=1, position=tpos
;    oplot, uts, new_efluxs, color=target_color
;
;    dy = (tpos[3]-tpos[1])/3
;
;    tx = tpos[2]+xchsz*1
;    ty0 = tpos[1]-ychsz*0.3
;    msg = 'Traced'
;    ty = ty0+dy*1
;    xyouts, tx,ty,msg, normal=1, color=target_color
;    msg = 'O+ eflux'
;    ty = ty0+dy*2
;    xyouts, tx,ty,msg, normal=1
    
    
    
    
    
    ; Plot results.
    down_pos = all_poss[*,nvar]
    down_poss = sgcalcpos(1,2, xpans=[6,4], xpad=1, position=down_pos)
    tpos = down_poss[*,0]
    xrange = time_range
    xrange[1] = xrange[0]+(tpos[2]-tpos[0])/(down_pos[2]-down_pos[0])*(xrange[1]-xrange[0])
    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
    xtitle = ' '
    ytitle = '|R| (Re)'
    ystep = 2
    yrange = [1,6]
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 2
    xstep = 600
    xtickv = make_bins(xrange, xstep, inner=1)
    xtickn = time_string(xtickv,tformat='hhmm')
    xtickn[0] += '!C'+time_string(xtickv[0],tformat='YYYY MTH DD')
    xticks = n_elements(xtickv)-1
    xminor = 10
    target_color = sgcolor('orange')

    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
        ystyle=1, yrange=yrange, ytitle=ytitle, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, $
        position=tpos, nodata=1, noerase=1
    foreach target_energy, target_energys, target_id do begin
        target_time = target_times[target_id]
        if target_time ge max_target_time then continue
        plots, trace_times[target_id,*], fline_diss, color=target_color, linestyle=0
    endforeach
    plots, xrange, trace_dis+[0,0], linestyle=1
    foreach tx, trace_tr do plots, tx+[0,0], yrange, linestyle=0, color=bar_color
    fig_label = fig_letters[nvar]+'-1)'
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = fig_label
    xyouts, tx,ty,msg, normal=1

    msg = 'Trace O+!Cbackward in!Ctime & altitude'
    tx = tpos[2]-xchsz*1
    ty = tpos[1]+ychsz*3
    xyouts, tx,ty,msg, normal=1, charsize=label_size, alignment=1
    


    tpos = down_poss[*,1]
    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
    xtitle = 'Time spread (sec)'
    ytitle = ' '
    xrange = [0,100]
    
    plot, trace_dtimes, fline_diss, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, ytitle=ytitle, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, ytickformat='(A1)', $
        noerase=1, position=tpos
    plots, xrange, trace_dis+[0,0], linestyle=1
    tx = tpos[2]-xchsz*1
    ty = tpos[1]+ychsz*3
    msg = 'Mininum!Cspread!C~'+string(trace_dis,format='(F3.1)')+' Re'
    xyouts, tx,ty,msg, normal=1, charsize=label_size, alignment=1

    fig_label = fig_letters[nvar]+'-2)'
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = fig_label
    xyouts, tx,ty,msg, normal=1
    

    poss = reform(all_poss[*,0:nvar-1])
    foreach var, plot_vars, var_id do begin
        tpos = poss[*,var_id]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, 'xticklen', xticklen
        options, 'yticklen', yticklen
    endforeach
    
    tplot, plot_vars, trange=time_range, position=poss, noerase=1, single_line_uttick=1
    timebar, trace_tr, color=bar_color
    for pid=0,nvar-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]-ychsz*0.7
        msg = fig_labels[pid]
        xyouts, tx,ty,msg, normal=1
    endfor
    
    test_var = prefix+'o_en_spec_para'
    get_data, test_var, times, fluxs, energys, limits=lim
    pid = where(plot_vars eq test_var, count)
    tpos = poss[*,pid]
    yrange = lim.yrange
    xrange = time_range
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, xlog=0, $
        ystyle=5, yrange=yrange, ylog=1, $
        nodata=1, noerase=1, position=tpos
    plots, target_times, target_energys, psym=1, symsize=0.5, color=target_color
    msg = 'PA [0,45] deg'
    tx = tpos[0]+xchsz*0.5
    ty = tpos[1]+ychsz*0.3
    xyouts, tx,ty,normal=1, msg, charsize=label_size
    
    tx0 = max_target_time+60
    txs = [min(target_times),tx0]
    ty0 = tpos[1]+ychsz*1.3
    foreach tx,txs,ii do begin
        tmp = convert_coord(tx,10, data=1, to_normal=1)
        txs[ii] = tmp[0]
        plots, tmp[0]+[0,0], ty0+[-1,1]*ychsz*0.2, normal=1
    endforeach
    tys = ty0+[0,0]
    plots, txs,tys, normal=1
    msg = 'Para acc'
    tx = mean(txs)
    ty = tys[0]-ychsz*1
    xyouts, tx,ty,msg, normal=1, charsize=label_size, alignment=0.5
    
    txs = [tx0,max(target_times)]
    foreach tx,txs,ii do begin
        tmp = convert_coord(tx,10, data=1, to_normal=1)
        txs[ii] = tmp[0]
        plots, tmp[0]+[0,0], ty0+[-1,1]*ychsz*0.2, normal=1
    endforeach
    tys = ty0+[0,0]
    plots, txs,tys, normal=1
    msg = 'Para acc or Perp heating'
    tx = mean(txs)
    ty = tys[0]-ychsz*1
    xyouts, tx,ty,msg, normal=1, charsize=label_size, alignment=0.5

    

;;---Conics.
;    tpos = down_poss[*,0]
;    xrange = time_range
;    xrange[1] = xrange[0]+(tpos[2]-tpos[0])/(down_pos[2]-down_pos[0])*(xrange[1]-xrange[0])
;    yrange = [1,6]
;    
;    plot, xrange, yrange, $
;        xstyle=5, xrange=xrange, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
;        ystyle=5, yrange=yrange, ytitle=ytitle, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, $
;        position=tpos, nodata=1, noerase=1
;        
;    pitch_angle = 15d*rad
;    ;pitch_angle = 30d*rad
;    mu = sin(pitch_angle)^2/bmags[0]
;    trace_times = dblarr(ntarget,nfline)
;    foreach fline_dis, fline_diss, fline_id do begin
;        if fline_id eq 0 then begin
;            trace_times[*,fline_id] = target_times
;            continue
;        endif
;
;        vmags = sqrt(2*target_energys/mass0)*1e-3 ; in km/s.
;        bmag = mean(bmags[fline_id-1:fline_id])
;        pitch_angle_step = asin(sqrt(mu*bmag))
;        vpara_steps = vmags*cos(pitch_angle_step)
;
;        r_step = reform(flines[fline_id,*]-flines[fline_id-1,*])
;        dr_step = snorm(r_step)*constant('re')      ; in km.
;        dt_steps = dr_step/vpara_steps
;        trace_times[*,fline_id] = trace_times[*,fline_id-1]-dt_steps
;    endforeach
;
;    conic_color = sgcolor('green')
;    foreach target_energy, target_energys, target_id do begin
;        target_time = target_times[target_id]
;        if target_time ge max_target_time then continue
;        plots, trace_times[target_id,*], fline_diss, color=conic_color, linestyle=0
;    endforeach


    if keyword_set(test) then stop
    sgclose
    return, plot_file

end

print, fig_2015_0416_0800_outflow_v01(event_info=event_info)
end