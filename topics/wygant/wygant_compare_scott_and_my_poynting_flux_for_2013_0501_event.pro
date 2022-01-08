
pro wygant_compare_scott_and_my_poynting_flux_for_2013_0501_event

    rootdir = shomedir()+'/Google Drive/works/data/wygant/rbspb_2013_0501_pflux_scott'
    tr = time_string(['2013-05-01/07:30','2013-05-01/07:50'])
    tr = time_string(['2013-05-01/11:00','2013-05-01/11:40'])
    tr = time_string(['2013-05-01/06:55','2013-05-01/08:15'])
    
    sgindexcolor, 43
    tplot_options, 'labflag', 1
    store_data, '*', /delete
    posl = [0.1,0.15,0.40,0.95]
    posr = [0.6,0.15,0.90,0.95]
    
    ; load original field.
    
    
    ; load scott's results.
    bands = ['11-to-55','65-to-295','295-to-1801']
    fns = rootdir+'/RBSPb_mapped_pflux_2013-05-01_'+bands+'-sec.dat'
    nfn = n_elements(fns)
    for i = 0, nfn-1 do begin
        labels = strjoin(strsplit(bands[i],'-to-',/extract),'-')+'s'
        labels = 'Scott!C  S!D||!N mapped!C  '+labels
        restore, filename = fns[i]
        store_data, 'scotts_pflux_f'+string(i+1,format='(I0)'), etimes, s_para_mapped, $
            limits = {labels:labels, ytitle:'(mW/m!U2!N)'}
    endfor
    vars = 'scotts_pflux'
    stplot_total, 'scotts_pflux_f'+string(indgen(nfn)+1,format='(I0)'), newname = vars
    options, vars, 'labels', 'Scott!C  S!D||!N mapped!C  total'

    ; load sheng's results.
    fn = shomedir()+'/Google Drive/works/data/rbsp_thg/rbsp_efw_fld_2013_0501_04.tplot'
    tplot_restore, filename = fn
    stplot_renew, 'rbspb_pf_fac_mat', newname = 'shengs_pflux'
    vars = 'rbspb_'+['e_gse','b_mod_gse','pos_gsm','pf_fac_mat']
    store_data, vars, /delete
    
    ; mat info.
    scaleinfo = [10,1000,40]
    dezr = [-0.2,0.2]
    dbzr = [-4,4]
    filters = [1000,207,55,10]
    dename = 'rbspb_de_fac'
    dbname = 'rbspb_db_fac'
    pfname = 'shengs_pflux'
    stplot_calc_poynting_flux, dename, dbname, pfname, $
        filter = filters, scaleinfo = scaleinfo
    
    ; pick out the parallel component.
    vars = pfname+'_mat'+['','f1','f2','f3']
    newnames = 'shengs_pflux'+['','_f1','_f2','_f3']
    for i = 0, n_elements(vars)-1 do begin
        stplot_index, vars[i], 0, newname = newnames[i]
        store_data, vars[i], /delete
    endfor
    options, newnames, 'colors', -1
    options, newnames, 'ytitle', '(mW/m!U2!N)'
    
    ; map to ionosphere.
    get_data, 'rbspb_fpt_lonlat', t0, dat
    store_data, 'rbspb_mlat', t0, dat[*,1]
    for i = 0, n_elements(newnames)-1 do begin
        smap2iono, newnames[i], 'rbspb_mlat', b = 'rbspb_b_gse'
        get_data, newnames[i], limits = tmp
        options, newnames[i], 'labels', 'Sheng!C  S!D||!N mapped!C  '+tmp.labels
    endfor
    options, 'shengs_pflux', 'labels', 'Sheng!C  S!D||!N in situ!C  total'
    vars = 'rbspb_'+['fpt_lonlat*','b_gse','bmod_gse','mlat']
    store_data, vars, /delete
    
    ; decompose dE, dB.
    stplot_split, 'rbspb_de_fac'
    stplot_split, 'rbspb_db_fac'
    
    ; delete vars not used.
    vars = 'rbspb_'+['de_fac_'+['comp?_mat','matf?','mat'],'db_fac_'+['comp?_mat','matf?','mat']]
    store_data, vars, /delete
    
;    ofn = shomedir()+'/2013_0501_rbspb_fields.eps'
;    sgpsopen, ofn, xsize = 7, ysize = 6, /inch
;    vars = 'rbspb_de_fac'+['','_comp'+['1','2','3']]
;    stplot_minmax, vars, [-50,80]
;    nvar = n_elements(vars)
;    poss = sgcalcpos(position = posl, nvar)
;    tplot, vars, trange = tr, position = poss, title = 'Original dE', /noerase
;    vars = 'rbspb_db_fac'+['','_comp'+['1','2','3']]
;    stplot_minmax, vars, [-50,50]
;    nvar = n_elements(vars)
;    poss = sgcalcpos(position = posr, nvar)
;    tplot, vars, trange = tr, position = poss, title = 'Original dB', /noerase    
;    sgpsclose, /pdf
    
    
;    stplot_minmax, ['scotts','shengs']+'_pflux', [-100,300]
;    stplot_minmax, ['scotts','shengs']+'_pflux_f1', [-50,20]
;    stplot_minmax, ['scotts','shengs']+'_pflux_f2', [-50,250]
;    stplot_minmax, ['scotts','shengs']+'_pflux_f3', [-20,50]
    
    ofn = shomedir()+'/2013_0501_poynt_compare.pdf'
    sgopen, ofn, xsize = 7, ysize = 6, /inch
    
    vars = 'scotts_pflux'+['','_f1','_f2','_f3']
    options, vars, 'constant', 0
    nvar = n_elements(vars)
    poss = sgcalcpos(position = posl, nvar)
    tplot, vars, trange = tr, position = poss, title = "Scott's version", /noerase
    
    vars = 'shengs_pflux'+['','_f1','_f2','_f3']
    options, vars, 'constant', 0
    options, vars, 'colors', 0
    nvar = n_elements(vars)
    poss = sgcalcpos(position = posr, nvar)
    tplot, vars, trange = tr, position = poss, title = "Sheng's version", /noerase
    
    sgclose

    ofn = shomedir()+'/2013_0501_poynt_sheng.pdf'
    sgopen, ofn, xsize = 7, ysize = 9.5, /inch

    vars = 'shengs_pflux'+['_f1','_f2','_f3','']
    nvar = n_elements(vars)
    pos = [0.15,0.1,0.85,0.9]
    poss = sgcalcpos(position = pos, nvar)
    options, vars, 'constant', 0
    
    sgindexcolor, 43
    tplot, vars, trange = tr, position = poss, /noerase
    sgclose

end