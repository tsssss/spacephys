;+
; Check ULF wave properties.
;-


function micro_injection_fig_ulf_v02, input_event_id, probe=probe, $
    plot_dir=plot_dir, test=test, get_name=get_name, update=update

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
    base = project_id+'_fig_ulf_'+event_id+'_mms'+probe+'_'+version+'.pdf'
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
    prefix = 'mms'+probe+'_'
    default_coord = 'gsm'
    mission_probe = 'mms'+probe

    ; omni IMF Bz.
    imf_var = omni_read_sw_b(time_range)
    options, imf_var, constant=0, yrange=[-1,1]*8, ytickv=[-1,0,1]*5, yminor=5, yticks=2


    ; Particle related vars.
    ele_pad_var = lets_read_this(func='mms_read_pad_ele_all', $
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load e thermal or kev data ...'
        return, retval
    endif

    ion_pad_var = lets_read_this(func='mms_read_pad_ion_all', $
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load ion thermal data ...'
        return, retval
    endif
    ion_kev_pad_var = lets_read_this(func='mms_read_pad_ion_kev', $
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load ion kev data ...'
        return, retval
    endif

    
    ; en spec.
    en_spec_vars = list()
    pa_spec_vars = list()
    pad_vars = prefix+['e_pad_thermal','e_pad_kev','ion_pad_thermal','p_pad_kev']
    foreach pad_var, pad_vars do begin
        en_spec_vars.add, pad_get_en_spec(pad_var=pad_var)
        var = pad_get_pa_spec(pad_var=pad_var)
        pa_spec_vars.add, var
        options, var, 'energy_range', minmax(get_var_setting(pad_var,'en_centers'))
    endforeach
    en_spec_vars = en_spec_vars.toarray()
    pa_spec_vars = pa_spec_vars.toarray()
    
    ele_pa_en_low = 80
    ele_pa_en = 800
    pad_var = prefix+'e_pad_thermal'
    ele_en_max = max(get_var_setting(pad_var, 'en_centers'))
    energy_range = [ele_pa_en_low,ele_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_low')
    options, tmp, 'energy_range', energy_range
    
    energy_range = [ele_pa_en,ele_en_max]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_high')
    options, tmp, 'energy_range', energy_range
    
    ion_pa_en_low = 80
    ion_pa_en = 4000
    pad_var = prefix+'ion_pad_thermal'
    ion_en_max = max(get_var_setting(pad_var, 'en_centers'))
    energy_range = [ion_pa_en_low,ion_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_low')
    options, tmp, 'energy_range', energy_range

    energy_range = [ion_pa_en,ion_en_max]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_high')
    options, tmp, 'energy_range', energy_range
    

    log_ytickv = [2,3,4]
    yticks = n_elements(log_ytickv)-1
    zrange = [1e3,1e8]
    log_ztickv = [3,4,5,6,7,8]
    zticks = n_elements(log_ztickv)-1
    ztickn = '10!U'+string(log_ztickv,format='(I0)')
    ztickn[0:*:2] = ' '
    vars = prefix+'e_en_spec_thermal'
    options, vars, zrange=zrange, constant=[ele_pa_en_low,ele_pa_en], $
        ytickv=10d^log_ytickv, yticks=yticks, ytickname='10!U'+string(log_ytickv,format='(I0)'), $
        ztickv=10d^log_ztickv, zticks=zticks, ztickname=ztickn, zminor=9, yminor=9
        
    
    vars = prefix+'e_en_spec_kev'
    log_ytickv = [4,5]
    yfactor = 5
    yrange = yfactor*10d^log_ytickv
    ytickn = string(yfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_ytickv,format='(I0)')
    ytickv = 1e5
    yticks = n_elements(ytickv)-1
    ytickn = '10!U5'
    log_zrange = [1,5]-1
    zfactor = 5
    zrange = zfactor*10d^log_zrange
    log_ztickv = make_bins(log_zrange,1, inner=1)
    ztickv = 10d^log_ztickv
    zticks = n_elements(ztickv)-1
    ztickn = string(zfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_ztickv,format='(I0)')
    index = where(log_ztickv eq 1, count)
    if count ne 0 then ztickn[index] = '10'
    ztickn[1:*:2] = ' '
    options, vars, zrange=zrange, yrange=yrange, $
        ytickv=ytickv, yticks=yticks, ytickname=ytickn, $
        ztickv=ztickv, zticks=zticks, ztickname=ztickn, zminor=9, yminor=9

    vars = prefix+'ion_pa_spec_thermal_low'
    zfactor = 1.3
    log_zrange = [4,5]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, vars, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    
    vars = prefix+'ion_pa_spec_thermal_high'
    zfactor = 8
    log_zrange = [4,5]-1
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, vars, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    vars = prefix+'ion_en_spec_thermal'
    log_zrange = [4,6]
    zfactor = 1
    zrange = zfactor*10d^log_zrange
    log_yrange = [1,4]
    yfactor = 2
    yrange = yfactor*10d^log_yrange
    log_ytickv = make_bins(minmax(alog10(yrange)),1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    options, vars, constant=[ion_pa_en_low,ion_pa_en], $
        yrange=yrange, ytickv=ytickv, yticks=yticks, ytickname=ytickn, yminor=9, $
        zrange=zrange
    
    var = prefix+'e_pa_spec_kev'
    zfactor = 3.5
    log_zrange = [2,3]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, var, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    var = prefix+'e_pa_spec_thermal_high'
    zfactor = 2.8
    log_zrange = [5,6]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, var, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    var = prefix+'e_pa_spec_thermal_low'
    zfactor = 1.5
    log_zrange = [6,7]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, var, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    
    var1 = mms_read_density_fpi(time_range, probe=probe, species='e')
    var2 = mms_read_density_fpi(time_range, probe=probe, species='p')
    dens_combo_var = stplot_merge([var1,var2], output=prefix+'density_fpi')
    options, [var1,var2], ylog=1
    options, dens_combo_var, ylog=1
    dens_var = var1
    var1 = mms_read_temperature_fpi(time_range, probe=probe, species='e')
    var2 = mms_read_temperature_fpi(time_range, probe=probe, species='p')
    options, [var1,var2], ylog=1
    
    ion_vel_var = lets_read_this(func='mms_read_ion_vel', $
        time_range, probe=mission_probe, errmsg=errmsg)
    
    
;---fields.
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

    ; Model related vars.
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
    

    ; Convert to FAC.
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
    
    ; displacement.
    dr_fac = vexb_fac
    dr_fac[*] = 0
    index = where(finite(snorm(vexb_fac)), count)
    vexb_fac = sinterpol(vexb_fac[index,*], times[index], times)
    index = where(finite(snorm(vexb_fac),nan=1), count)
    vexb_fac[index,*] = 0
    ;for ii=0,ndim-1 do vexb_fac[*,ii] = interpol(vexb_fac[index,ii], times[index], times)
    dt = sdatarate(times)
    foreach time, times, tid do begin
        if tid eq 0 then continue
        dr_fac[tid,*] = dr_fac[tid-1,*]+vexb_fac[tid,*]*dt
    endforeach
    re1 = 1d/constant('re')
    window = 60*60d
    width = window/dt
    for ii=0,ndim-1 do begin
        dr_fac[*,ii] *= re1
    endfor
    dr_mag = snorm(dr_fac)
    dr_mag -= smooth(dr_mag, width)
    var = prefix+'dr_mag'
    store_data, var, times, dr_mag
    add_setting, var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', '|dR|', $
        'unit', 'Re' )

    var = prefix+'dr_mms_fac'
    store_data, var, times, dr_fac
    add_setting, var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'dR', $
        'unit', 'Re' )
    
    
    ; Settings.
    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    vars = prefix+['b_mms_fac','e_mms_fac','vexb_mms_fac','u_mms_fac','dr_mms_fac']
    options, vars, 'labels', fac_labels
    
    var = prefix+'b_mms_fac'
    options, var, yrange=[-5,45], ytickv=[0,20,40], yminor=4, yticks=2, constant=[0,20,40]
    var = prefix+'e_mms_fac'
    options, var, yrange=[-1,1]*4, ytickv=[-1,0,1]*3, yminor=3, yticks=2, constant=[0]
    var = prefix+['vexb_mms_fac','u_mms_fac']
    options, var, yrange=[-1,1]*180, ytickv=[-1,0,1]*150, yminor=5, yticks=2, constant=[-1,0,1]*100
    var = prefix+'b_gsm'
    get_data, var, times, vec, limits=lim
    mag = snorm(vec)
    store_data, var+'_plot', times, [[vec],[mag]], limits=lim
    ;options, var, yrange=[-10,42], ytickv=[0,15,30], yminor=3, yticks=2, constant=[0,15,30], $
    options, var+'_plot', yrange=[-8,30], ytickv=[0,15,30], yminor=3, yticks=2, constant=[0,15,30], $
        labels=['B!D'+['x','y','z'],'|B|'], colors=[constant('rgb'),sgcolor('black')]
    

;---Pressure.
    pmag_var = lets_calc_pmag(b_var=b_gsm_var)
    tavg_var = prefix+'p_tavg'
    get_data, prefix+'p_t_fpi', times, temp
    store_data, tavg_var, times, temp[*,0]
    pth_var = lets_calc_pthermal(n_var=prefix+'p_density_fpi',t_var=tavg_var, var_info=prefix+'p_pth')
    tavg_var = prefix+'e_tavg'
    get_data, prefix+'e_t_fpi', times, temp
    store_data, tavg_var, times, temp[*,0]
    eth_var = lets_calc_pthermal(n_var=prefix+'e_density_fpi',t_var=tavg_var, var_info=prefix+'e_pth')
    pcombo_var = prefix+'p_combo'
    pth = get_var_data(pth_var, times=times)
    eth = get_var_data(eth_var, at=times)
    pth_total = pth+eth
    pb = get_var_data(pmag_var, at=times)
    ptot = pth_total+pb
    store_data, pcombo_var, times, [[ptot],[pth_total],[pb]]
    add_setting, pcombo_var, smart=1, dictionary($
        'display_type', 'stack', $
        'unit', 'nPa', $
        'ylog', 0, $
        'labels', 'P!D'+['tot','th','B'], $
        'colors', constant('rgb') )
    options, pcombo_var, yrange=[0,0.8], ystyle=1, ytickv=[0,0.4,0.8], yticks=2, yminor=4
    
    colors = sgcolor(['red','blue'])
    labels = ['Ele','Ion']
    ; density
    dens_var = stplot_merge(prefix+['e','p']+'_density_fpi', $
        output=prefix+'dens_combo', colors=colors, labels=labels)
    options, dens_var, yrange=[0.2,2], ylog=1, $
        ytickv=[0.2,2], yticks=1, yminor=9, ytickname=['0.2','2'], ytitle='(cm!U-3!N)'
    
    temp_var = stplot_merge(prefix+['e','p']+'_tavg', $
        output=prefix+'tavg_combo', colors=colors, labels=labels)
    options, temp_var, yrange=[1.2e2,2e4], ylog=1, $
        ytickv=[1e3,1e4], yticks=1, yminor=9, ytickname='10!U'+['3','4'], ytitle='(eV)'
    options, prefix+'e_mms_fac', ytitle='(mV/m)'
    options, prefix+'p_combo', ytitle='(nPa)'
    options, prefix+['vexb_mms_fac','u_mms_fac'], ytitle='(km/s)'
    options, prefix+'dr_mag', yrange=[-1,1]*1.4, ytickv=[-1,0,1]*1, yminor=5, yticks=2, constant=[-1,0,1]*1

;---Plot settings.
    pa_suffix = ''
    plot_vars = [$
        imf_var, $
        prefix+[$
        ;'e_en_spec_kev','e_en_spec_thermal', $
        'b_gsm_plot','e_mms_fac','vexb_mms_fac','u_mms_fac','dr_mag', $
        'dens_combo','tavg_combo','p_combo']]
    panel_labels = [$
        'IMF', $
        ;'Ele EN high','Ele EN',$
        'B GSM','E FAC','V!DExB!N','V!Dion!N','dR', $
        'N','T','P']
    nplot_var = n_elements(plot_vars)
    fig_letters = letters(nplot_var)
    
    ypads = [4,0.4+fltarr(nplot_var-2)]
    pansize = [6,0.75]
    margins = [10,4,10,1]
    poss = panel_pos(plot_file, nypan=nplot_var, ypads=ypads, $
        pansize=pansize, fig_size=fig_size, margins=margins)
    
    plot_tr = list()
    plot_tr.add, time_double(['2015-09-01/15:00','2015-09-01/17:00'])
    plot_tr.add, time_double(['2015-09-01/18:00','2015-09-01/21:30'])
    nxpan = n_elements(plot_tr)
    xpans = fltarr(nxpan)
    foreach tr, plot_tr, tid do begin
        xpans[tid] = total(tr*[-1,1])
    endforeach
    pos0 = sgcalcpos(1,nxpan, xpans=xpans, position=poss[*,1], xpad=5+fltarr(nxpan-1))
    
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    if n_elements(abs_ticklen) eq 0 then abs_ticklen = -0.3*ychsz*fig_size[1]
    if n_elements(abs_xticklen) eq 0 then abs_xticklen = abs_ticklen
    if n_elements(abs_yticklen) eq 0 then abs_yticklen = abs_ticklen

    
    tpos = poss[*,0]
    xticklen = abs_xticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = abs_yticklen/(tpos[2]-tpos[0])/fig_size[0]
    options, imf_var, xticklen=xticklen, yticklen=yticklen
    tplot_options, 'tickinterval', 3600*1
    tplot, imf_var, trange=time_range, noerase=1, position=tpos
    msg = fig_letters[0]+') '+panel_labels[0]
    tx = tpos[0]-xchsz*6
    ty = tpos[3]-ychsz*0.8
    xyouts, tx,ty,msg, normal=1
    
    thick = (keyword_set(test))? 2: 8
    set_axis, position=tpos, xrange=time_range, yrange=[0,1]
    foreach tr, plot_tr, tid do begin
        txs = tr
        tys = 0
        plots, txs, tys, thick=thick, data=1
    endforeach
    
    
    tplot_options, 'tickinterval', 1800
    bottom_poss = poss[*,1:nplot_var-1]
    bottom_vars = plot_vars[1:nplot_var-1]
    my_letters = fig_letters[1:nplot_var-1]
    my_labels = panel_labels[1:nplot_var-1]
    
    ; save limits.
    plot_lims = list()
    foreach var, bottom_vars do plot_lims.add, get_var_setting(var)
    
    
    foreach tr, plot_tr, tid do begin
        my_pos = bottom_poss
        my_pos[0,*] = pos0[0,tid]+xchsz*2
        my_pos[2,*] = pos0[2,tid]+xchsz*2
        
        foreach var, bottom_vars, vid do begin
            lim = plot_lims[vid]
            tpos = my_pos[*,vid]
            xticklen = abs_xticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = abs_yticklen/(tpos[2]-tpos[0])/fig_size[0]
            
            ; ytitle shown in the first column.
            ytitle = (tid eq 0)? lim.ytitle: ' '
            ytickformat = (tid eq 0)? '': '(A1)'
            ; labels shown in the last column.
            labels = (tid eq nxpan-1)? lim.labels: strarr(n_elements(lim.labels))
            options, var, ytitle=ytitle, ytickformat=ytickformat, labels=labels, $
                xticklen=xticklen, yticklen=yticklen
        endforeach
        
       
        tplot, bottom_vars, position=my_pos, noerase=1, trange=tr, novtitle=1
        foreach var, bottom_vars, vid do begin
            tpos = my_pos[*,vid]

            msg = my_letters[vid]+'-'+string(tid+1,format='(I0)')+') '
            if tid eq 0 then msg += my_labels[vid]
            
            tx = tpos[0]-xchsz*10
            if tid ne 0 then tx = tpos[0]-xchsz*3
            
            ty = tpos[3]-ychsz*0.8
            xyouts, tx,ty,msg, normal=1
        endforeach
    endforeach
    
    foreach var, bottom_vars, vid do begin
        lim = (plot_lims[vid]).tostruct()
        foreach key, tag_names(lim), kid do options, var, key, lim.(kid)
    endforeach

    
    if keyword_set(test) then stop
    sgclose

    return, plot_file
    
end

probe = '4'
event_id = '2015_0901_18'
event_id = '2015_0901_10'
print, micro_injection_fig_ulf_v02(event_id, probe=probe, test=0, update=1)
end