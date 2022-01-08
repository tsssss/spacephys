pro test_on_poynting_flux_bands, eventid

    ; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')
    
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    info = cusp_read_conjun_list(logfile, event = eventid)
    if size(info,/type) ne 8 then begin
        print, 'no id found ...'
        return
    endif
    
    ; directory to save output data and pdfs.
    orootdir = shomedir()+'/cusp/test'
    if ~file_test(orootdir,/directory) then file_mkdir, orootdir
    
    potr = info.polar.plot_time
    fatr = info.fast.plot_time
    potrcusp = info.polar.cusp_time
    fatrcusp = info.fast.cusp_time
    if potr[1] lt potr[0] then potr[1]+= 86400d
    if fatr[1] lt fatr[0] then fatr[1]+= 86400d
    if potrcusp[1] lt potrcusp[0] then potrcusp[1]+= 86400d
    if fatrcusp[1] lt fatrcusp[0] then fatrcusp[1]+= 86400d
    id = info.id
    
    pofilt0 = info.polar.filters        ; original filters.
    pofilts = [0,pofilt0,1e31]          ; include 0 and max filters.
    pofilts = pofilts[sort(pofilts)]
    pofilts = pofilts[uniq(pofilts)]
    ponfilt = n_elements(pofilts)-1     ; # of bands.
    pofstrs = 'f'+string(indgen(ponfilt),format='(I0)')    ; string version.
    
    ; **** load and calc all needed data.
    ; po_b, po_ilat, po_mlt, po_dis, po_[de,db]_fac, po_[de,db,pf]_mat,
    ; po_[de,db,pf]_matf[0,1,...], po_pf_[mat,direct], po_pf_[mat,direct]_map
    black = 255
    get_data, 'po_pf_fac_mat', t0
    if n_elements(t0) eq 1 then reload = 1
    if keyword_set(reload) then begin
        store_data, '*', /delete
        ; only load Polar dEv and dBp.
        deidx = 0   ; dEv.
        dbidx = 1   ; dBp.
        fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'
        tinfo = info.polar
        polar_sdt_prep_poynting_flux, fn, e56 = tinfo.e56, tr = potr, $
            cusptr = potrcusp, orootdir = orootdir+'/'+id
        vars = ['po_b0_spc','po_spc2fac','po_db_spc','po_de_spc','po_mlat']
        store_data, vars, /delete
        vars = 'po_de_fac'
        stplot_index, vars, deidx, newname = vars
        options, vars, 'ytitle', '(mV/m)'
        options, vars, 'labels', 'dEv'
        options, vars, 'colors', black
        vars = 'po_db_fac'
        stplot_index, vars, dbidx, newname = vars
        options, vars, 'ytitle', '(nT)'
        options, vars, 'labels', 'dBp'
        options, vars, 'colors', black
        
        ; poynting flux.
        ; po_[de,db]_fac_mat, 3-D fields used in Poynting flux calc.
        ; po_[de,db]_fac_matf[1,2,3,...], 3-D fields in freq bands.
        ; po_de_fac_comp[deidx+1]_mat, po_db_fac_comp[dbidx+1]_mat, MAT spec.
        dename = 'po_de_fac' & dbname = 'po_db_fac' & pfname = 'po_pf_fac'
        stplot_calc_poynting_flux, dename, dbname, pfname, $
            filter = pofilts, idfilter = pofstrs, scaleinfo = tinfo.scaleinfo
        stplot_total, tnames('po_pf_fac_matf?'), newname = 'po_pf_fac_mat'
        
        vars = 'po_pf_fac_mat'
        smap2iono, vars, 'po_dis', newname = vars+'_map'
        cusp_int_eflux, vars+'_map', 'po_ilat', tr = potrcusp
        
        ; **** mulitply fields directly and integrate.
        get_data, 'po_de_fac', t0, de
        get_data, 'po_db_fac', t0, db
        store_data, 'po_pf_fac_direct', t0, spoynt(de, db)
        
        vars = 'po_pf_fac_direct'
        smap2iono, vars, 'po_dis', newname = vars+'_map'
        cusp_int_eflux, vars, 'po_ilat', tr = potrcusp
        cusp_int_eflux, vars+'_map', 'po_ilat', tr = potrcusp
        
        vars = [tnames('po_pf_fac_matf?'),'po_pf_fac_mat','po_pf_fac_direct']
        for i = 0, n_elements(vars)-1 do $
            cusp_int_eflux, vars[i], 'po_ilat', tr = potrcusp

    endif
    
    ; **** plot settings.
    charsz = 1
    !p.font = 1
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    tplot_options, 'charsize', charsz*0.9
    tplot_options, 'xcharsize', charsz*0.7
    tplot_options, 'ycharsize', charsz*0.8
    tplot_options, 'zcharsize', charsz*0.6
    time_stamp, /off
    labs = ['po_ilat','po_mlt','po_dis']
    ofn = orootdir+'/test_on_poynting_flux_bands.eps'
    sgpsopen, ofn, xsize = 8, ysize = 10, /inch
    sgindexcolor, 43
    black = 0
    red = 6
    erase
    ymax = 0.95 & ymin = 0.15
    podepos = [0.1,ymin,0.25,ymax]
    podbpos = [0.4,ymin,0.55,ymax]
    popfpos = [0.7,ymin,0.85,ymax]

    options, 'po_de_fac', 'labels', 'original'
    options, 'po_de_fac_mat', 'labels', 'MAT total'
    vars = ['po_de_fac','po_de_fac_mat',tnames('po_de_fac_matf?')]
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position = podepos)
    options, vars, 'ytitle', '(mV/m)'
    options, vars, 'colors', black
    tplot, vars, var_label = labs, position = pos, /noerase, title = 'dEv'
    timebar, potrcusp, color = red, thick = 2

    options, 'po_db_fac', 'labels', 'original'
    options, 'po_db_fac_mat', 'labels', 'MAT total'
    vars = ['po_db_fac','po_db_fac_mat',tnames('po_db_fac_matf?')]
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position = podbpos)
    options, vars, 'ytitle', '(nT)'
    options, vars, 'colors', black
    tplot, vars, var_label = labs, position = pos, /noerase, title = 'dBp'
    timebar, potrcusp, color = red, thick = 2

    vars = 'po_pf_fac_direct'
    get_data, vars, t0, tmp, dat
    options, vars, 'labels', 'direct product!C  Int = '+snum2str(round(dat))+' W/m'
    vars = 'po_pf_fac_mat'
    get_data, vars, t0, tmp, dat
    options, vars, 'labels', 'use MAT!C  Int = '+snum2str(round(dat))+' W/m'
    vars = ['po_pf_fac_direct','po_pf_fac_mat',tnames('po_pf_fac_matf?')]
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position = popfpos)
    options, vars, 'ytitle', '(mV/m!U2!N)'
    options, vars, 'colors', black
    tplot, vars, var_label = labs, position = pos, /noerase, title = 'S!D||!N'
    timebar, potrcusp, color = red, thick = 2

    sgpsclose, /pdf
end

id = '1998_0925_05'
test_on_poynting_flux_bands, id
end