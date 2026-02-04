;+
; Do wavelet transform for various quantities.
;-

function micro_injection_fig_wavelet_v02, input_event_id, test=test, get_name=get_name, update=update, errmsg=errmsg

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
    base = project_id+'_fig_wavelet_'+event_id+'_'+version+'.pdf'
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
    mission = 'mms'
    default_coord = 'gsm'
    external_model = 't89'
    internal_model = 'igrf'
    fac_coord = mission+'_fac'
    fac_labels = ['b','w','o']
    b0_window = 1200d
    probes = string(findgen(4)+1,format='(I0)')
    probes = '4'
    nprobe = n_elements(probes)
    colors = sgcolor(['red','green','blue','purple'])
    labels = strupcase('mms'+probes)
    comps = constant('xyz')
    ncomp = n_elements(comps)
    
    foreach probe, probes do begin
        prefix = mission+probe+'_'
        source = ['mms',probe]
        r_var = lets_read('orbit', time_range, source=source, coord=default_coord)
        b_var = lets_read('bfield', time_range, source=source, coord=default_coord)
        e_var = lets_read('efield', time_range, source=source, coord=default_coord)
        u_var = lets_read('ion_vel', time_range, source=source, coord=default_coord)
        mission_probe = mission+probe
        ele_kev_en_spec_var = lets_read_this(func='mms_read_en_spec_ele', $
            time_range, probe=mission_probe, id='kev', errmsg=errmsg)
        
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
        add_setting, ele_flux_var, smart=1, dictionary($
            'display_type', 'stack', $
            'short_name', 'F', $
            'labels', string(the_energys*1e-3,format='(I0)')+' keV', $
            'ylog', 0, $
            'ytitle', 'Detrended!CLog!D10!Nflux' )
        
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
    endforeach


;---Make plot.
;tplot, prefix+['e_en_spec_kev','ele_flux','ele_flux1_mor','ele_kev_phase']


    prefix = 'mms4_'
    plot_vars = prefix+['b1','e','u']+'_mms_fac_mor'
    labels = ['dB','E','Ion Vel']+'!D'+tex2str('perp')+',west!N'
    
    vars = prefix+['b1','e','u']+'_mms_fac'
    options, vars, 'labels', ['||',tex2str('perp')+','+['west','out']]

    plot_vars = prefix+['b1','e','u']+'_mms_fac_mor'
    tmp = ['','_mor']
    tmp = '_mor'
    plot_vars = prefix+[$
        'b1_mms_fac'+tmp, $
        'e_mms_fac'+tmp, $
        'u_mms_fac'+tmp ]
    labels = ['dB','E','Ion Vel']+'!D'+tex2str('perp')+',west!N'
    comp = ''
    tmp = [' Morlet']
    labels = [$
        'dB'+comp+tmp, $
        'E'+comp+tmp, $
        'Ion Vel'+comp+tmp ]
    plot_vars = [prefix+['e_en_spec_kev','ele_flux','ele_flux1_mor','ele_kev_phase'], $
        plot_vars]
    ;labels = ['Ele High','Ele Flux','Ele Flux!C    Morlet','Phase',labels]
    labels = ['Ele High','Log Flux','Ele Flux','Phase',labels]
    
    vars = prefix+['ele_flux1_mor','ele_kev_phase','b1_mms_fac_mor','e_mms_fac_mor','u_mms_fac_mor']
    options, vars, ystyle=1, ylog=1, ytitle='Freq!C(mHz)', spec=1, display_type='spec', yrange=[0.3,30]
    vars = prefix+'u_mms_fac_mor'
    options, vars, ztitle='U (km/s)!U2!N'

    plot_tr = time_range

    tplot_options, 'tickinterval', 3600
    tmp = sgplot(plot_vars, filename=plot_file, xrange=time_range, panel_labels=labels)
    if keyword_set(test) then stop
    sgclose

    return, plot_file


end

input_event_id = '2015_0901_11'
test = 0
input_event_id = '2015_0901_10'
print, micro_injection_fig_wavelet_v02(input_event_id, test=test, update=1)
end