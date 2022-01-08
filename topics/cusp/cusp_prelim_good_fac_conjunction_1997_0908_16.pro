;+
; This event has good correlation in FAC shape from the field data.
; How about the Poynting flux asscociated with the FAC?
;-
pro cusp_prelim_good_fac_conjunction_1997_0908_16, reload = reload

; **** basic info.
    id = '1997_0908_16'
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')

    potr = time_double(['1997-09-08/15:58','1997-09-08/16:19'])
    potrcusp = time_double(['1997-09-08/16:06','1997-09-08/16:12'])
    poscinfo = [1,600,60]
    pofilts = [25.8,1e5]

    fatr = time_double(['1997-09-08/15:32','1997-09-08/15:40'])
    fatrcusp = time_double(['1997-09-08/15:34:53','1997-09-08/15:37:14'])
    fascinfo = [1,300,60]
    fafilts = [15,1e5]

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
    tplot_options, 'constant', 0
    red = 6
    posu = [0.15,0.60,0.90,0.95]
    posd = [0.15,0.10,0.90,0.45]
    
    if keyword_set(reload) then begin
        store_data, '*', /delete

        ; **** load polar and fast field data.
        fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'
        polar_sdt_prep_poynting_flux, fn
        vars = ['po_b','po_b0_spc','po_spc2fac','po_db_spc','po_de_spc','po_mlat']
        store_data, vars, /delete
        
        fn = rootdir+'/data/cusp/fa_sdt_fld_'+id+'_*.tplot'
        fn = file_search(fn)
        fast_sdt_prep_poynting_flux, fn
        vars = ['fa_pos','fa_vel','fa_alt','alt','fa_b0_gei','fa_b0_gei','fa_b']
        store_data, vars, /delete
        
        ; **** poynting flux.
        dename = 'po_de_fac' & dbname = 'po_db_fac' & pfname = 'po_pf_fac'
        stplot_calc_poynting_flux, dename, dbname, pfname, $
            filter = pofilts, scaleinfo = poscinfo
        smap2iono, 'po_pf_fac_mat', 'po_dis', newname = 'po_pf_fac_map'
        cusp_int_eflux, 'po_pf_fac_map', 'po_ilat', tr = potrcusp
        
        dename = 'fa_de_fac' & dbname = 'fa_db_fac' & pfname = 'fa_pf_fac'
        stplot_calc_poynting_flux, dename, dbname, pfname, $
            filter = fafilts, scaleinfo = fascinfo
        smap2iono, 'fa_pf_fac_mat', 'fa_dis', newname = 'fa_pf_fac_map'
        cusp_int_eflux, 'fa_pf_fac_map', 'fa_ilat', tr = fatrcusp
        
        ; **** map polar and fast quantities to uniform ilat.
        
    endif

    vars = ['po_pf_fac_map','po_de_fac_mat','po_db_fac_mat', $
        'fa_pf_fac_map','fa_de_fac_mat','fa_db_fac_mat']
    options, vars, 'colors', [6,4,2]
    options, vars, 'labels', ['v','p','b']
    vars = ['po_','fa_']+'de_fac_mat'
    options, vars, 'ytitle', '(mV/m)'
    vars = ['po_','fa_']+'db_fac_mat'
    options, vars, 'ytitle', '(nT)'
    vars = ['po_','fa_']+'pf_fac_map'
    options, vars, 'ytitle', '(mW/m!U2!N)'
        
    ofn = shomedir()+'/cusp_prelim_good_fac_conjunction_1997_0908_16.eps'
    sgpsopen, ofn, xsize = 5, ysize = 8, /inch

    sgindexcolor, 43
    !p.color = 0
    
    pre = 'po_'
    labs = pre+['mlt','dis','ilat']
    titl = 'Polar FAC Poynting flux, dE, dB'
    vars = pre+['pf_fac_map','de_fac_mat','db_fac_mat']
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, position = posu)
    tplot, vars, trange = potr, /noerase, position = poss, $
        var_label = labs, /short_vlabel, title = titl
    timebar, potrcusp, color = red
    
    pre = 'fa_'
    labs = pre+['mlt','dis','ilat']
    titl = 'FAST FAC Poynting flux, dE, dB'
    vars = pre+['pf_fac_map','de_fac_mat','db_fac_mat']
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, position = posd)
    tplot, vars, trange = fatr, /noerase, position = poss, $
        var_label = labs, /short_vlabel, title = titl
    timebar, fatrcusp, color = red

    sgpsclose, /pdf
end

cusp_prelim_good_fac_conjunction_1997_0908_16;, /reload
end
