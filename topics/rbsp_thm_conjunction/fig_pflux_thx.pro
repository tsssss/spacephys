
    ;_2014_1222_load_data, /reload
    tplot_options, 'ystyle', 1
    
    utr1 = time_double(['2014-12-22/04:00','2014-12-22/06:30'])
    utr3 = time_double(['2014-12-22/03:00','2014-12-22/08:00'])
    utr2 = time_double(['2014-12-22/04:15','2014-12-22/04:45'])
    
    pres = ['tha','the','thd']+'_'
    foreach pre0, pres do begin
        get_data, pre0+'e_gsm', uts, egsm
        get_data, pre0+'b_gsm', tuts, bgsm
        bgsm = sinterpol(bgsm, tuts, uts)
        store_data, pre0+'b_gsm', uts, bgsm
        
        scaleinfo = {s0:4d, s1:1200d, dj:1d/8, ns:0d}    ; pflux calc.
        stplot_calc_pflux_mor, pre0+'e_gsm', pre0+'b_gsm', pre0+'pf_gsm', scaleinfo=scaleinfo
            
        
        get_data, pre0+'pf_gsm', uts, pfgsm
        get_data, pre0+'map_coef0', tuts, mapc
        mapc = interpol(mapc[*,0], tuts, uts)
        store_data, pre0+'pf_para', uts, sdot(pfgsm,sunitvec(bgsm));*mapc, $
            limits={ytitle:'(mW/m!U2!N)', labels:'S!D||!N@100 km'}
        ;tplot, pre0+['e_gsm','b_gsm','pf_gsm','pf_para'], trange=utr1
    endforeach
    
    

;---Plot.
    ofn = shomedir()+'/fig_pflux_thx.pdf'
    ofn = 0
    sgopen, ofn, xsize=6, ysize=8, /inch
    
    device, decomposed=0
    loadct2, 43
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    tvar = ['e_gsm','b_gsm','pf_gsm','pf_para']
    vars = []
    foreach pre0, pres do vars = [vars, pre0+tvar]
    
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar)
    
    tplot, vars, position=poss, trange=utr1
    
    sgclose
stop    

  ;-zoom in to pi2.
    ofn = shomedir()+'/fig_pflux_thx_pi2.pdf'
    ;ofn = 0
    sgopen, ofn, xsize=6, ysize=8, /inch
    
    device, decomposed=0
    loadct2, 43
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    tvar = ['e_gsm','b_gsm','pf_gsm','pf_para']
    vars = []
    foreach pre0, pres do vars = [vars, pre0+tvar]
    
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar)
    
    tplot, vars, position=poss, trange=utr2
    
    sgclose


end
