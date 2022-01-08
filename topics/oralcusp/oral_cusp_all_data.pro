
pro int_eflux, vname, ilatname, eflux_int, tr = tr
    re = 6378.137d & h0 = 100 ; km, altitude.
    c = (re+h0)*!dpi/180d
    
    get_data, vname, t0, eflux
    get_data, ilatname, tmp, ilat
    ilat = sinterpol(ilat, tmp, t0)
    nrec = n_elements(t0)
    
    idx = keyword_set(tr)? where(t0 ge tr[0] and t0 le tr[1]): indgen(nrec)
    eflux = eflux[idx] & ilat = ilat[idx] & nrec = n_elements(eflux)
    v1 = 0.5*(eflux[0:nrec-2]+eflux[1:nrec-1])
    v2 = abs(ilat[1:nrec-1]-ilat[0:nrec-2])*c
    eflux_int = total(v1*v2)
    print, vname, string(eflux_int,format='(f10.2)'), ' W/m'
end

pro oral_cusp_all_data

    store_data, '*', /delete

    ; **** polar data.
    pre = 'po_'
    trpo = time_double(['1998-09-25/05:00','1998-09-25/06:00'])
    cusptrpo = time_double(['1998-09-25/05:27','1998-09-25/05:40'])

    polar_sdt_prep_poynting_flux, sdiskdir('Works')+'/data/cusp/po_sdt_fld_1998_0925_05.sdt'
    tref = time_double('1998-09-25/05:00')
    read_polar_ke, sdiskdir('Works')+'/confs/agu/agu2013/data/polar_hyd_kei.dat', pre+'kei', tref
    read_polar_ke, sdiskdir('Works')+'/confs/agu/agu2013/data/polar_hyd_kee.dat', pre+'kee', tref

    ; poynting flux, map and integrate.
    stplot_calc_poynting_flux, pre+'de_fac', pre+'db_fac', pre+'pf_fac', $
        method = 'mat', filter = [10,40,180,1500], scaleinfo = [0.5,2000,60]
    vars = tnames(pre+['pf_fac_mat','pf_fac_matf?']) & nvar = n_elements(vars)
    for i = 0, nvar-1 do stplot_index, vars[i], 2, newname = vars[i]+'_para'
    vars = [pre+['kei','kee'],vars+'_para']
    nvar = n_elements(vars) & inteflux = dblarr(nvar)
    for i = 0, nvar-1 do begin
        smap2iono, vars[i], pre+'dis', newname = vars[i]+'_map'
        int_eflux, vars[i]+'_map', pre+'ilat', tmp, tr = cusptrpo & inteflux[i] = tmp
    endfor
    store_data, pre+'int_eflux', dblarr(nvar), inteflux
    options, pre+'int_eflux', 'labels', vars
    
    ; **** fast data.
    pre = 'fa_'
    trfa = time_double(['1998-09-25/04:27','1998-09-25/04:32'])
    cusptrfa = time_double(['1998-09-25/04:27:55','1998-09-25/04:29:35'])

    tplot_restore, filename = sdiskdir('Works')+'/data/cusp/fa_sdt_esa_1998_0925_08278.tplot'
    fast_sdt_prep_poynting_flux, sdiskdir('Works')+'/data/cusp/fa_sdt_fld_1998_0925_08278.tplot'
    
    ; rename stuff.
    stplot_renew, 'ion_eflux', newname = pre+'kei'
    stplot_renew, 'ele_eflux', newname = pre+'kee'
    
    ; poynting flux, map and integrate.
    stplot_calc_poynting_flux, pre+'de_fac', pre+'db_fac', pre+'pf_fac', $
        method = 'mat', filter = [3,10,30,100], scaleinfo = [0.5,250,60]
    vars = tnames(pre+['pf_fac_mat','pf_fac_matf?']) & nvar = n_elements(vars)
    for i = 0, nvar-1 do stplot_index, vars[i], 2, newname = vars[i]+'_para'
    vars = [pre+['kei','kee'],vars+'_para']
    nvar = n_elements(vars) & inteflux = dblarr(nvar)
    for i = 0, nvar-1 do begin
        smap2iono, vars[i], pre+'dis', newname = vars[i]+'_map'
        int_eflux, vars[i]+'_map', pre+'ilat', tmp, tr = cusptrfa & inteflux[i] = tmp
    endfor
    store_data, pre+'int_eflux', dblarr(nvar), inteflux
    options, pre+'int_eflux', 'labels', vars
    
    ; print eflux.
    vars = tnames('*_int_eflux')
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], tmp, inteflux, alimits = tmp2
        for j = 0, n_elements(tmp)-1 do $
            print, tmp2.labels[j], string(inteflux[j]*1e-3,format='(f10.2)'), ' kW/m'
    endfor

    ; save data to file.
    tplot_save, ['po_*','fa_*'], filename = shomedir()+'/oral_cusp_data'
end