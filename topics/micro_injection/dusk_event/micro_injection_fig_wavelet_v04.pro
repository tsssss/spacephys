;+
; Figure is the same as v03, just to update the codes and added comments.
;-
function micro_injection_fig_wavelet_v04, input_event_id, probe=probe, $
    plot_dir=plot_dir, test=test, get_name=get_name, update=update

    errmsg = ''
    retval = !null
    version = 'v04'
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
    base = project_id+'_fig_wavelet_'+event_id+'_mms'+probe+'_'+version+'.pdf'
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
    ion_pad_var = lets_read_this(func='mms_read_pad_ion_thermal_fpi', $
        time_range, probe=mission_probe, errmsg=errmsg)
    ion_pad_var = rename_var(ion_pad_var, output=prefix+'p_pad_thermal')


;---Settings.
    unit_type = 'nflux'
    ele_pa_en_low = 80
    ele_pa_en = 800
    ion_pa_en_low = 80
    ion_pa_en = 4000
    ion_pa_high = 2.5e4
    ele_pa_high = ion_pa_high

    
;---Derived data.
    suffix = '_'+['thermal','kev','all']
    pad_vars = prefix+['e_pad'+suffix,'p_pad'+suffix]
    unit0 = '/cm!U2!N-s-sr-keV'
    foreach pad_var, pad_vars do begin
        var1 = pad_var+'_flux'

        ; Adjust unit.
        data = get_var_data(pad_var, times=times, limits=lim)
        vals = lim.en_centers
        foreach val, vals, vid do begin
            if unit_type eq 'nflux' then begin
                ; do nothing.
                unit = '#'+unit0
            endif else if unit_type eq 'eflux' then begin
                data[*,*,vid] *= val*1e-3
                unit = 'keV'+unit0
            endif else if unit_type eq 'xflux' then begin
                data[*,*,vid] *= sqrt(val*1e-3)
                unit = 'keV!U0.5!N'+unit0
            endif
        endforeach

        ; adjust flux for ion.
        if pad_var eq prefix+'p_pad_all' then begin
            index = where(vals ge 2.5e4)
            data[*,*,index] *= 5
        endif

        ; remove some useless energy bins
        index = where(vals le 1.5e5)
        data = data[*,*,index]
        vals = vals[index]

        ; save data.
        store_data, var1, times, data, limits=lim
        options, var1, unit=unit, en_centers=vals
    endforeach


    ; en spec.
    en_spec_vars = list()
    pa_spec_vars = list()
    foreach pad_var, pad_vars do begin
        en_spec_vars.add, pad_get_en_spec(pad_var=pad_var)
        var = pad_get_pa_spec(pad_var=pad_var)
        pa_spec_vars.add, var
        options, var, 'energy_range', minmax(get_var_setting(pad_var,'en_centers'))
    endforeach
    en_spec_vars = en_spec_vars.toarray()
    pa_spec_vars = pa_spec_vars.toarray()

    ; pa spec.
    ; electron low and mid.
    pad_var = prefix+'e_pad_thermal'
    ele_en_max = max(get_var_setting(pad_var, 'en_centers'))
    energy_range = [ele_pa_en_low,ele_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_low')
    options, tmp, 'energy_range', energy_range
    ; electron high.
    energy_range = [ele_pa_en,ele_en_max]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_high')
    options, tmp, 'energy_range', energy_range
    
    ; ion low and mid.
    pad_var = prefix+'ion_pad_thermal'
    ion_en_max = max(get_var_setting(pad_var, 'en_centers'))
    energy_range = [ion_pa_en_low,ion_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_low')
    options, tmp, 'energy_range', energy_range
    ; ion high.
    energy_range = [ion_pa_en,ion_en_max]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_high')
    options, tmp, 'energy_range', energy_range
    

;---Settings for en and pa spec.
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
    
    vars = prefix+'p_en_spec_thermal'
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
    
    
;---Load more data
    mission = 'mms'
    default_coord = 'gsm'
    external_model = 't89'
    internal_model = 'igrf'
    fac_coord = mission+'_fac'
    fac_labels = ['b','w','o']
    b0_window = 1200d
    colors = sgcolor(['red','green','blue','purple'])
    comps = constant('xyz')
    ncomp = n_elements(comps)
    
    probe = '4'
    prefix = mission+probe+'_'
    source = ['mms',probe]
    r_var = lets_read('orbit', time_range, source=source, coord=default_coord)
    b_var = lets_read('bfield', time_range, source=source, coord=default_coord)
    e_var = lets_read('efield', time_range, source=source, coord=default_coord)
    u_var = lets_read('ion_vel', time_range, source=source, coord=default_coord)
    mission_probe = mission+probe
    ele_kev_en_spec_var = lets_read_this(func='mms_read_en_spec_ele', $
        time_range, probe=mission_probe, id='kev', errmsg=errmsg)
    fluxs = get_var_data(ele_kev_en_spec_var, en_bins, times=times, settings=settings)


    vars = prefix+'e_en_spec_kev'
    log_ytickv = [4,5]
    yfactor = 3
    yrange = yfactor*10d^log_ytickv
    yrange = [3e4,5e5]
    ytickn = string(yfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_ytickv,format='(I0)')
    ytickv = 1e5
    yticks = n_elements(ytickv)-1
    ytickn = '10!U5'
    log_zrange = [1,5]-1
    zfactor = 8
    zrange = zfactor*10d^log_zrange
    log_ztickv = make_bins(log_zrange,1, inner=1)
    ztickv = 10d^log_ztickv
    zticks = n_elements(ztickv)-1
    ztickn = string(zfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_ztickv,format='(I0)')
    index = where(log_ztickv eq 1, count)
    if count ne 0 then ztickn[index] = '10'
    ztickn[1:*:2] = ' '
    ct = 64
    ct = 56
    ct = 40
    options, vars, zrange=zrange, yrange=yrange, $
        ytickv=ytickv, yticks=yticks, ytickname=ytickn, $
        ztickv=ztickv, zticks=zticks, ztickname=ztickn, zminor=9, yminor=9, color_table=ct

    
    ele_flux_var = prefix+'ele_flux'
    test_energys = [5e4,8e4]
    test_energys = [5e4,1e5]
    ntest_energy = n_elements(test_energys)
    fluxs = get_var_data(ele_kev_en_spec_var, energy_bins, times=times, settings=settings)
    ntime = n_elements(times)
    the_fluxs = dblarr(ntime,ntest_energy)
    the_energys = dblarr(ntest_energy)
    foreach energy, test_energys, eid do begin
        tmp = min(energy_bins-energy, abs=1, energy_index)
        the_fluxs[*,eid] = fluxs[*,energy_index]
        the_energys[eid] = energy_bins[energy_index]
    endforeach
    window = 1200d
    width = window/sdatarate(times)
    foreach energy, test_energys, eid do begin
        data = alog10(the_fluxs[*,eid])
        data -= smooth(data,nan=1, width, edge_truncate=1)
        the_fluxs[*,eid] = data
    endforeach
    store_data, ele_flux_var, times, the_fluxs
    yrange = [-1,1]*0.8
    ytickv = [-1,0,1]*0.5
    yticks = n_elements(ytickv)-1
    yminor = 5
    add_setting, ele_flux_var, smart=1, dictionary($
        'display_type', 'stack', $
        'short_name', 'F', $
        'labels', string(the_energys*1e-3,format='(I0)')+' keV', $
        'ylog', 0, $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'constant', ytickv, $
        'ytitle', 'Detr.!CLog!D10!Nflux' )
    
    ; Calc wavelet.
    scale_info = {s0:2d, s1:4000, dj:1d/8, ns:0d }
    foreach tid, [0,1] do begin
        flux_var = ele_flux_var+string(tid+1,format='(I0)')
        store_data, flux_var, times, the_fluxs[*,tid]
        spec_var = stplot_mor_new(flux_var, scale_info=scale_info)
        data = get_var_data(spec_var, freqs)
        store_data, spec_var, times, data, freqs*1e3
        unit = settings.unit
        ztitle = 'Ele flux [Log!D10!Nflux]!U2'
        yrange = [0.3,30]
        options, spec_var, zrange=[1e-5,1e0], color_table=40, $
            ytitle='Freq!C(mHz)', yrange=yrange, ylog=1, ztitle=ztitle
    endforeach
    
    s0 = scale_info.s0
    s1 = scale_info.s1
    dj = scale_info.dj
    j1 = floor(alog(s1/s0)/alog(2)/dj)  ; # of powers-of-two with dj
    s1 = s0*2d^(dj*j1)            
    ns = j1+1
    j1 = ns-1
    w0 = 6d
    cdelta = 0.776d     ; constant for w0=6, for normalization.
    dr0 = sdatarate(times)
    
    fa = the_fluxs[*,0]
    fb = the_fluxs[*,1]
    mora = wavelet(fa, dr0, pad=1, s0=s0, dj=dj, j=j1, $
        mother='Morlet', param=w0, $
        period = ps, scale=ss, coi=coi)
    morb = wavelet(fb, dr0, pad=1, s0=s0, dj=dj, j=j1, $
        mother='Morlet', param=w0, $
        period = ps, scale=ss, coi=coi)
    morab = mora*conj(morb)
    phase = atan(imaginary(morab)/real_part(morab))*constant('deg')
    phase_var = prefix+'ele_kev_phase'
    store_data, phase_var, times, phase, freqs*1e3
    unit = 'deg'
    ztitle = 'Phase ('+unit+')'
    yrange = [0.3,30]
    options, phase_var, color_table=70, $
        ytitle='Freq!C(mHz)', yrange=yrange, ylog=1, ystyle=1, $
        ztitle=ztitle, spec=1, zlog=0, zrange=[-1,1]*90



    bmod_var = lets_read_geopack_bfield(orbit_var=r_var, external_model=external_model, internal_model=internal_model)
    b_vars = lets_decompose_bfield(b_var=b_var, b0_window=b0_window, bmod_var=bmod_var)
    b0_var = b_vars['b0']
    b1_var = b_vars['b1']

    ; Convert to FAC.
    q_fac_var = lets_define_fac(r_var=r_var, b_var=b0_var, fac_coord=fac_coord)
    coord_msgs = [default_coord,fac_coord]
    fac_vars = list()
    foreach var, [b1_var,e_var,u_var] do begin
        options, var, mission='mms'
        fac_vars.add, lets_cotran(coord_msgs, input=var, q_var=q_fac_var)
    endforeach

    ; Calc wavelet.
    wanted_comp = 'w'
    wanted_index = where_pro(fac_labels, 'eq', wanted_comp)
    ct = 40
    zrange = [1e-2,1e3]
    foreach var, fac_vars do begin
        mor_var = var+'_mor'
        if tnames(mor_var) ne '' then continue
        vec = get_var_data(var, times=times, settings=settings)
        dat = vec[*,wanted_index]
        field_var = var+'_tmp'
        store_data, field_var, times, dat
        spec_var = stplot_mor_new(field_var, scale_info=scale_info)
        get_data, spec_var, times, data, freqs
        store_data, spec_var, times, data, freqs*1e3
        unit = settings.unit
        short_name = settings.short_name
        add_setting, spec_var, smart=1, dictionary($
            'requested_time_range', time_range, $
            'no_interp', 1, $
            'display_type', 'spec', $
            'unit', unit, $
            'ytitle', 'Freq (mHz)', $
            'yrange', minmax(freqs*1e3), $
            'ylog', 1, $
            'zlog', 1, $
            'short_name', short_name )
        ztitle = short_name+'('+unit+')!U2!N'
        options, spec_var, ztitle=ztitle, color_table=ct, zrange=zrange
        tmp = strpos(var, 'u_mms_fac')
        if tmp[0] ne -1 then options, spec_var, zrange=zrange*10
        tmp = rename_var(spec_var,output=var+'_mor')
    endforeach

    


;---Settings for plots.
    plot_vars = prefix+['b1','e','u']+'_mms_fac_mor'
    labels = ['dB','E','Ion Vel']+'!D'+tex2str('perp')+',west!N'
    
    vars = prefix+['b1','e','u']+'_mms_fac'
    options, vars, 'labels', ['||',tex2str('perp')+','+['west','out']]
    


    vars = prefix+['ele_flux1_mor','ele_kev_phase','b1_mms_fac_mor','e_mms_fac_mor','u_mms_fac_mor']
    options, vars, ystyle=1, ylog=1, ytitle='Freq!C(mHz)', spec=1, display_type='spec', yrange=[0.3,30]
    vars = prefix+'u_mms_fac_mor'
    options, vars, ztitle='U (km/s)!U2!N'
    
    var = prefix+'ele_kev_phase'
    get_data, var, times, mors, freqs
    freq_range = [1.8,4] ; mHz
    freq_index = where_pro(freqs, '[]', freq_range, count=count)
    if count eq 0 then message, 'Inconsistency ...'
    ntime = n_elements(times)
    phase_diff = fltarr(ntime)
    ;phase_diff = mean(mors[*,freq_index],dimension=2)
    for tid=0,ntime-1 do begin
        phase_diff[tid] = median(mors[tid,freq_index])
    endfor
    window = 600
    time_step = sdatarate(times)
    width = window/time_step
    phase_diff = smooth(phase_diff, width, edge_mirror=1)
    vars = prefix+'phase_diff'
    store_data, vars, times, phase_diff
    ytickv = [-1,0,1]*50
    yrange = [-1,1]*70
    add_setting, vars, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'Phase Diff', $
        'unit', 'deg', $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', 2, $
        'yminor', 5, $
        'constant', ytickv )


;---Make plot.
    plot_info = orderedhash()

    var = prefix+'e_en_spec_kev'
    plot_info[var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'Ele High', $
        'setting', dictionary() )
    var = prefix+'ele_flux'
    plot_info[var] = dictionary($
        'routine', 'plot_line', $
        'panel_label_text', 'Ele flux', $
        'setting', dictionary() )
    var = prefix+'phase_diff'
    plot_info[var] = dictionary($
        'routine', 'plot_line', $
        'panel_label_text', 'Phase', $
        'setting', dictionary() )

;    plot_vars = [prefix+['e_en_spec_kev','ele_flux','phase_diff','ele_kev_phase','ele_flux1_mor']]
;    labels = ['Ele High','Log Flux','Phase','Phase Mor','Ele Mor']
    plot_vars = plot_info.keys()
    nplot_var = n_elements(plot_vars)

    panel_letters = letters(nplot_var)
    foreach plot_var, plot_vars, pid do begin
        my_info = plot_info[plot_var]
        if ~my_info.haskey('ypan') then my_info['ypan'] = 1d
        if ~my_info.haskey('panel_label_text') then my_info['panel_label_text'] = ' '
        if ~my_info.haskey('panel_letter') then my_info['panel_letter'] = panel_letters[pid]
        if ~my_info.haskey('panel_label_msg') then my_info['panel_label_msg'] = my_info['panel_letter']+') '+my_info['panel_label_text']
    endforeach
    
    var_labels = prefix+['mlat','dis','mlt']
    var_labels = !null
    nvar_label = n_elements(var_labels)
    options, prefix+'mlat', ytitle='MLat (deg)'
    options, prefix+'dis', ytitle='|R| (Re)'
    options, prefix+'mlt', ytitle='MLT (h)'
    var = prefix+'mlt'
    get_data, var, times, data
    index = where(data le 0, count)
    if count ne 0 then begin
        data[index] += 24
        store_data, var, times, data
    endif

    margins = [12,3.5+nvar_label,8,2]
    ypans = []
    foreach plot_var, plot_vars, pid do begin
        ypans = [ypans,(plot_info[plot_var])['ypan']]
    endforeach
    ; need 1 row for freq spec and zoom in flux.
    ; need 1 row for pa spec.
    nypan = nplot_var+2
    ypans = [ypans,1.2,5.8]
    ypads = [0.4+fltarr(nplot_var-1),4,1]
    all_poss = panel_pos(plot_file, nypan=nypan, fig_size=fig_size, $
        ypans=ypans, pansize=[6,0.8], ypads=ypads, margins=margins)
    plot_poss = all_poss[*,0:nplot_var-1]
    
    ; Use positions to determine [x,y]ticklen.
    abs_ticklen = 0.3
    foreach plot_var, plot_vars, pid do begin
        my_info = plot_info[plot_var]
        my_info['position'] = plot_poss[*,pid]
        my_info['abs_ticklen'] = abs_ticklen
        ;plot_info[plot_var] = my_info
    endforeach
    
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    plot_tr = time_range
    tickinterval = 3600
    tplot_options, 'tickinterval', tickinterval

    foreach plot_var, plot_vars, pid do begin
        my_info = plot_info[plot_var]
    
        my_pos = my_info['position']
        plot_routine = my_info['routine']
        plot_setting = my_info['setting']
        plot_setting['position'] = my_pos
        plot_setting['noerase'] = (pid eq 0)? 0: 1
        plot_setting['xtickformat'] = (pid eq nplot_var-1)? '': '(A1)'
        plot_setting['novtitle'] = (pid eq nplot_var-1)? 0: 1
        plot_setting['time_range'] = plot_tr
        plot_setting['tickinterval'] = tickinterval
        plot_setting['panel_label_pos'] = [xchsz*1,my_pos[3]-ychsz*0.7]
        plot_setting['panel_label_msg'] = my_info['panel_label_msg']
        plot_setting['var_labels'] = var_labels
        plot_setting['vlab_margin'] = margins[0]-1
        plot_setting['trange'] = plot_tr
        foreach key, plot_setting.keys() do begin
            index = strpos(key,'tick_setting')
            if index[0] ne -1 then begin
                tick_setting = plot_setting[key]
                foreach comp, constant('xyz') do begin
                    the_comp = comp+'tickv'
                    if tick_setting.haskey(the_comp) then begin
                        tick_setting[comp+'ticks'] = n_elements(tick_setting[the_comp])-1
                    endif
                endforeach
            endif
        endforeach

        plot_setting = plot_setting.tostruct()        
        tmp = call_function(plot_routine, plot_var, _extra=plot_setting)
        
        if pid eq 0 then begin
            msg = strupcase('mms'+probe)
            color = sgcolor('black')
            tpos = my_pos
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            xyouts, tx,ty,msg, normal=1, color=color
        endif
    endforeach


    tplot_options, get_options=opts
    str_element, opts, 'tickinterval', delete=1
    
    
;---Freq panel.
    the_pos = all_poss[*,nplot_var]
    bottom_poss = sgcalcpos(1,2, position=the_pos, xpans=[1,1.5], xpad=15)
    vars = prefix+['ele_flux1_mor',['e']+'_mms_fac_mor']
    msgs = ['Ele flux','E field']
    colors = sgcolor(['red','blue'])
    xrange = [0.2,50]
    yrange = [0,4.5]
    ytickv = [0,1,2,3,4]
    yticks = n_elements(ytickv)-1
    yminor = 2
    xtitle = 'Freq (mHz)'
    ytitle = 'X/stddev(X)'
    tpos = bottom_poss[*,0]
    abs_xticklen = -abs_ticklen*ychsz*fig_size[1]
    abs_yticklen = -abs_xticklen
    xticklen = abs_xticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = abs_yticklen/(tpos[2]-tpos[0])/fig_size[0]
    label_size = 0.8
    plot, xrange, yrange, $
        xstyle=5, xlog=1, xrange=xrange, xtitle=xtitle, $
        ystyle=5, ylog=0, yrange=yrange, ytitle=ytitle, $
        position=tpos, noerase=1, nodata=1
    tx = tpos[0]-xchsz*10
    ty = tpos[1]-ychsz*1.5+xticklen*(tpos[3]-tpos[1])
    xyouts, tx,ty, xtitle, alignment=0, normal=1
    
    foreach var, vars, vid do begin
        color = colors[vid]
        get_data, var, times, mors, freqs, limits=lims
        cwt_info = lims.cwt_info
        txs = cwt_info.fs*1e3
        tys = cwt_info.psd
        ty0 = stddev(tys)

        if vid eq 0 then begin
            max_txs = list()
            foreach tmp, txs, tid do begin
                if tid eq 0 then continue
                if tid eq n_elements(txs)-1 then continue
                if tys[tid] ge tys[tid-1] and tys[tid] ge tys[tid+1] then max_txs.add, tid
            endforeach
            max_txs = txs[max_txs.toarray()]
            foreach max_tx, max_txs, ttid do begin
                oplot, max_tx+[0,0], yrange, linestyle=1, color=color
                tx = max_tx
                ty = yrange[0]
                tmp = convert_coord(tx,ty, data=1, to_normal=1)
                tx = tmp[0]
                ty = tmp[1]-ychsz*1
                msg = string(max_tx,format='(F4.1)')
                if ttid eq 0 then msg += ' mHz'
                xyouts, tx,ty,msg, normal=1, alignment=0.5, charsize=label_size, color=color
                
                the_p = 1d3/max_tx/60
                msg = string(the_p,format='(I0)')
                if ttid eq 0 then msg += ' min'
                ty = tpos[3]-ychsz*label_size
                xyouts, tx,ty,msg, normal=1, alignment=0.5, charsize=label_size, color=color
            endforeach
        endif

        oplot, txs,tys/ty0, color=color
        tx = tpos[2]-xchsz*1
        ty = tpos[3]-ychsz*(1+vid)
        msg = msgs[vid]
        xyouts, tx,ty,msg, alignment=1, color=color, normal=1
    endforeach


    plot, xrange, yrange, $
        xstyle=1, xlog=1, xrange=xrange, xtitle=' ', $
        ystyle=1, ylog=0, yrange=yrange, ytitle=ytitle, ytickv=ytickv, yticks=yticks, yminor=yminor, $
        position=tpos, noerase=1, nodata=1, $
        xticklen=xticklen, yticklen=yticklen
    tx = tpos[0]-xchsz*9
    ty = tpos[3]-ychsz*0.7
    msg = 'd) PSD'
    xyouts, tx,ty,msg, normal=1
    
    
;---Phase diff zoom in.
    tickinterval = 600
    tplot_options, 'tickinterval', tickinterval
    zoom_tr = time_double(['2015-09-01/20:20','2015-09-01/20:50'])
    zoom_tr = time_double(['2015-09-01/19:10','2015-09-01/19:40'])
    zoom_tr = time_double(['2015-09-01/18:00','2015-09-01/18:30'])
    tpos = bottom_poss[*,1]
    plot_var = prefix+'ele_flux'
    my_info = plot_info[plot_var]
    my_pos = tpos
    plot_routine = my_info['routine']
    plot_setting = my_info['setting']
    plot_setting['position'] = my_pos
    pid = nplot_var-1
    plot_setting['noerase'] = (pid eq 0)? 0: 1
    plot_setting['xtickformat'] = (pid eq nplot_var-1)? '': '(A1)'
    plot_setting['novtitle'] = (pid eq nplot_var-1)? 0: 1
    plot_setting['time_range'] = zoom_tr
    plot_setting['tickinterval'] = tickinterval
    plot_setting['panel_label_pos'] = [my_pos[0]-xchsz*10,my_pos[3]-ychsz*0.7]
    plot_setting['panel_label_msg'] = my_info['panel_label_msg']
    plot_setting['var_labels'] = var_labels
    plot_setting['vlab_margin'] = 8
    plot_setting['trange'] = zoom_tr
    plot_setting['panel_label_msg'] = 'e) Zoom-in'
    plot_setting['single_line_uttick'] = 1
    plot_setting = plot_setting.tostruct()        
    tmp = call_function(plot_routine, plot_var, _extra=plot_setting)
    
    
    top_pos = plot_poss[*,-1]
    top_pos[3] = plot_poss[3,0]
    yrange = [0,1]
    xrange = plot_tr
    set_axis, position=top_pos, xrange=xrange, yrange=yrange
    foreach tx, zoom_tr, tid do begin
        plots, tx+[0,0], yrange, linestyle=2
        tmp = convert_coord(tx,yrange[0], data=1, to_normal=1)
        txs = [tmp[0],tpos[0]]
        if tid eq 1 then txs = [tmp[0],tpos[2]]
        tys = [tmp[1],tpos[3]]
        plots, txs,tys, normal=1, linestyle=2
    endforeach
    
    
    

    ; Add 5 min label.
    min_label_color = sgcolor('orange_red')
    tpos = bottom_poss[*,1]
    set_axis, position=tpos, xrange=zoom_tr, yrange=[0,1]
    ;tx = time_double('2015-09-01/20:39')
    tx = time_double('2015-09-01/19:21')
    tx = time_double('2015-09-01/18:12:40')
    del_x = 300
    txs = tx+[0,del_x]
    foreach tx,txs,tid do begin
        tmp = convert_coord(tx,yrange[1], data=1, to_normal=1)
        txs[tid] = tmp[0]
    endforeach
    tys = tpos[3]-ychsz*2
    plots, txs, tys, normal=1, color=min_label_color
    tx = mean(txs)
    ty = tys[0]+ychsz*0.2
    msg = string(del_x/60,format='(I0)')+' min'
    xyouts, tx,ty,normal=1, msg, alignment=0.5, color=min_label_color
    foreach tx, txs, tid do begin
        ttxs = tx+[0,0]
        ttys = tys[0]+[-1,1]*ychsz*0.2
        plots, ttxs,ttys, normal=1, color=min_label_color
    endforeach
    
    ; 22 min label.
    tpos = all_poss[*,1]
    tpos = all_poss[*,0]
    xrange = plot_tr
    yrange = [0,1]
    set_axis, position=tpos, xrange=xrange, yrange=yrange
    del_x = 22*60d
    tx = time_double('2015-09-01/16:25')
    tx = time_double('2015-09-01/17:58')
    txs = tx+[0,del_x]
    foreach tx,txs,tid do begin
        tmp = convert_coord(tx,yrange[1], data=1, to_normal=1)
        txs[tid] = tmp[0]
    endforeach
    tys = tpos[3]-ychsz*1.0
    plots, txs, tys, normal=1, color=min_label_color
    tx = mean(txs)
    ;ty = tys[0]-ychsz*1
    ty = tys[0]+ychsz*0.2
    msg = string(del_x/60,format='(I0)')+' min'
    xyouts, tx,ty,normal=1, msg, alignment=0.5, color=min_label_color
    foreach tx, txs, tid do begin
        ttxs = tx+[0,0]
        ttys = tys[0]+[-1,1]*ychsz*0.2
        plots, ttxs,ttys, normal=1, color=min_label_color
    endforeach
    
    ; 41 min label.
    del_x = 41*60d
    tx = time_double('2015-09-01/13:35')
    txs = tx+[0,del_x]
    foreach tx,txs,tid do begin
        tmp = convert_coord(tx,yrange[1], data=1, to_normal=1)
        txs[tid] = tmp[0]
    endforeach
    tys = tpos[3]-ychsz*1.0
    plots, txs, tys, normal=1, color=min_label_color
    tx = mean(txs)
    ty = tys[0]+ychsz*0.2
    msg = string(del_x/60,format='(I0)')+' min'
    xyouts, tx,ty,normal=1, msg, alignment=0.5, color=min_label_color
    foreach tx, txs, tid do begin
        ttxs = tx+[0,0]
        ttys = tys[0]+[-1,1]*ychsz*0.2
        plots, ttxs,ttys, normal=1, color=min_label_color
    endforeach


;---Add pitch angle 2D.
    ; Settings.
    ele_pad_var = prefix+'e_pad_all_flux'
    ion_pad_var = prefix+'p_pad_all_flux'
    options, [ele_pad_var,ion_pad_var], 'mission', 'mms'
    ct = 40
    ;ct = 64
    
    if unit_type eq 'nflux' then begin
        ele_zrange = [1e1,5e8]*0.9
        ion_zrange = [1e3,1e5]
    endif else if unit_type eq 'eflux' then begin
        ; for keV/xxx.
        ele_zrange = [1e3,1e7]*2
        ion_zrange = [1e2,1e6]*40
    endif else if unit_type eq 'xflux' then begin
        ; for keV^0.5/xxx.
        ele_zrange = [1e2,1e7]*4.5
        ion_zrange = [1e3,1e5]*5
    endif

    ncolor = 25
    letters = ['f','g']
    pad_times = [$
        '2015-09-01/19:58:00',$
        '2015-09-01/20:00:20',$
        '2015-09-01/20:02:40' ]
;        '2015-09-01/20:30:20',$
;        '2015-09-01/20:38:20' ]
    pad_times = [$
        '2015-09-01/19:21:00', $
        '2015-09-01/19:23:30', $
        '2015-09-01/19:26:00' ]
    pad_times = [$
        '2015-09-01/18:12:40', $
        '2015-09-01/18:15:20', $
        '2015-09-01/18:18:00' ]

    event['pad_times'] = time_double(pad_times)
    pad_times = event.pad_times
    npad_time = n_elements(pad_times)
    pa_nxpan = n_elements(pad_times)
    pa_nypan = 2

    ; add pad times to zoom in panel.
    thick = keyword_set(test)? 3:10
    tpos = bottom_poss[*,1]
    set_axis, position=tpos, xrange=zoom_tr, yrange=[0,1]
    tys = tpos[1]+[-1,1]*ychsz*0.3
    foreach tx,pad_times,tid do begin
        tmp = convert_coord(tx,yrange[1], data=1, to_normal=1)
        plots, tmp[0]+[0,0], tys, normal=1, thick=thick
        ty = mean(tys)-ychsz*1.2
        msg = 't'+string(tid+1,format='(I0)')
        xyouts, tmp[0],ty,msg, normal=1, alignment=0.5
    endforeach


    pa_pos = all_poss[*,nplot_var+1]
    pa_pos = [0,0,1,all_poss[3,nplot_var+1]-ychsz*1]
    pa_margins = [8,4,8,1.5]
    xpad = 1
    ypad = 1
    poss = sgcalcpos(pa_nypan,pa_nxpan, region=pa_pos, margins=pa_margins, xpad=1, ypad=0.5)


    for tid=0,npad_time-1 do begin
        pad_time = pad_times[tid]
        if keyword_set(test) then pad_plot_file = 0
        ;if file_test(pad_plot_file) eq 1 then continue
        
        no_cb = (tid eq npad_time-1)? 0: 1
        ytitle = (tid eq 0)? 'Perp E (eV)': ' '
        ytickformat = (tid eq 0)? '': '(A1)'
        xtitle = ' '
        xtickformat = '(A1)'
        title = strupcase('mms-'+probe)+' '+time_string(pad_time)+' UT'
        constants = [ele_pa_en_low,ele_pa_en,ele_pa_high]
        tmp = plot_pad_polygon(ele_pad_var, test=0, plot_times=pad_time, $
            ytitle=ytitle, ytickformat=ytickformat, $
            xtitle=xtitle, xtickformat=xtickformat, $
            zrange=ele_zrange, color_table=ct, ncolor=ncolor, $
            constants=constants, title=title, $
            position=poss[*,tid,0], no_colorbar=no_cb)
        
        
        xtitle = 'Para E (eV)'
        xtickformat = ''
        title = ''
        constants = [ion_pa_en_low,ion_pa_en,ion_pa_high]
        tmp = plot_pad_polygon(ion_pad_var, test=0, plot_times=pad_time, $
            ytitle=ytitle, ytickformat=ytickformat, $
            xtitle=xtitle, xtickformat=xtickformat, $
            zrange=ion_zrange, color_table=ct, ncolor=ncolor, $
            constants=constants, title=title, $
            position=poss[*,tid,1], no_colorbar=no_cb)
        ;tmp = plot_pad_polygon(ion_kev_pad_var, test=test, plot_times=pad_time)     
        
        
    ;---Add label for FEEPS and FPI.
        if tid eq 1 then begin
            color = sgcolor('red')
            foreach pid, [0,1] do begin
                tpos = poss[*,tid,pid]
                xrange = [-1,1]
                yrange = [-1,1]
                set_axis, position=tpos, xrange=xrange, yrange=yrange
                
                tr = 0.6
                tx =-tr
                ty = tr
                msg = 'FEEPS'
                xyouts, tx,ty,msg, color=color, data=1, alignment=0.5
                
                tr = 0.25
                tx =-tr
                ty = tr
                msg = 'FPI'
                xyouts, tx,ty,msg, color=color, data=1, alignment=0.5
            endforeach 
        endif
        
    ;---Add label for species.
        if tid eq 0 then begin
            msgs = ['Ele','Ion']
            foreach pid, [0,1] do begin
                tpos = poss[*,tid,pid]
                tx = tpos[0]-xchsz*6
                ty = tpos[3]-ychsz*0.8
                msg = msgs[pid]
                xyouts, tx,ty,msg, normal=1
            endforeach
        endif
        
        

    ;---Add label.
        foreach pid, [0,1] do begin
            tpos = poss[*,tid,pid]

            num = tid+1
            num_str = string(num,format='(I0)')
            msg = letters[pid]+'-'+num_str+')'
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            xyouts, tx,ty,msg, normal=1
        endforeach
    endfor
    
    
;---Add label for low and high energy electron.
    hsize = (keyword_set(test))? thick*2: thick*20
    tpos = poss[*,0,0]
    xrange = [-1,1]
    yrange = [-1,1]
    set_axis, position=tpos, xrange=xrange, yrange=yrange
    ty0 = 0d
    tx0 = -0.75
    txs = tx0+[0,0]
    tys = ty0+[0,0]
    foreach tx,txs,tid do begin
        tmp = convert_coord(tx,ty0,data=1,to_normal=1)
        txs[tid] = tmp[0]
        tys[tid] = tmp[1]
    endforeach
    tys[1] = tys[0]-ychsz*2
    ;plots, txs,tys, normal=1
    arrow, txs[1],tys[1],txs[0],tys[0], normal=1, solid=1, thick=thick*0.6, hsize=hsize
    tx = txs[1]
    ty = tys[1]-ychsz*0.9
    msg = 'High'
    xyouts, tx,ty,msg, normal=1, alignment=0.5
    
    
    tpos = poss[*,1,0]
    xrange = [-1,1]
    yrange = [-1,1]
    set_axis, position=tpos, xrange=xrange, yrange=yrange
    ty0 = 0d
    tx0 = -0.2
    txs = tx0+[0,0]
    tys = ty0+[0,0]
    foreach tx,txs,tid do begin
        tmp = convert_coord(tx,ty0,data=1,to_normal=1)
        txs[tid] = tmp[0]
        tys[tid] = tmp[1]
    endforeach
    tys[1] = tys[0]-ychsz*2
    ;plots, txs,tys, normal=1
    arrow, txs[1],tys[1],txs[0],tys[0], normal=1, solid=1, thick=thick*0.6, hsize=hsize
    tx = txs[1]
    ty = tys[1]-ychsz*0.9
    msg = 'Low'
    xyouts, tx,ty,msg, normal=1, alignment=0.5

    
    if keyword_set(test) then stop
    sgclose


    return, plot_file

end

probe = '4'
event_id = '2015_0901_10'
print, micro_injection_fig_wavelet_v04(event_id, probe=probe, test=0, update=1)
end