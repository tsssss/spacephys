
    ;_2014_1222_load_data, /reload
    tplot_options, 'ystyle', 1
    
    utr1 = time_double(['2014-12-22/04:00','2014-12-22/06:30'])
    utr2 = time_double(['2014-12-22/02:00','2014-12-22/09:00'])
    utr3 = time_double(['2014-12-22/04:25','2014-12-22/04:50'])
    
    pres = ['rbb','rba']+'_'
    foreach pre0, pres do begin
        tprobe = strmid(pre0,1,1,/reverse)
        efw = sread_rbsp_efw_l3(utr2, probe=tprobe)
        uts = sfmepoch(efw.epoch, 'unix')
        emgse = efw.efield_inertial_frame_mgse
        idx = where(uts ge utr2[0] and uts le utr2[1])
        uts = uts[idx]
        emgse = emgse[idx,*]
        emgsm = sgse2gsm(emgse, stoepoch(uts,'unix'))
        store_data, pre0+'e_mgsm', uts, emgsm, limits=$
            {ytitle:'(mV/m)',colors:[6,4,2],labels:'mGSM E'+['x','y','z']}        
        get_data, pre0+'b_gsm', tuts, bgsm
        bgsm = sinterpol(bgsm, tuts, uts)
        store_data, pre0+'b_gsm', uts, bgsm
        
        scaleinfo = {s0:20d, s1:1200d, dj:1d/8, ns:0d}    ; pflux calc.
        stplot_calc_pflux_mor, pre0+'e_mgsm', pre0+'b_gsm', pre0+'pf_gsm', scaleinfo=scaleinfo
        
        
        get_data, pre0+'pf_gsm', uts, pfgsm
        get_data, pre0+'map_coef0', tuts, mapc
        mapc = interpol(mapc[*,0], tuts, uts)
        store_data, pre0+'pf_para', uts, sdot(pfgsm,sunitvec(bgsm))*mapc, $
            limits={ytitle:'(mW/m!U2!N)', labels:'S!D||!N@100 km'}
        ;tplot, pre0+['e_gsm','b_gsm','pf_gsm','pf_para'], trange=utr1
    endforeach
    
    

;---Plot.
    ofn = shomedir()+'/fig_pflux_rbx.pdf'
    ;ofn = 0
    sgopen, ofn, xsize=6, ysize=8, /inch
    
    device, decomposed=0
    loadct2, 43
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    tvar = ['e_mgsm','b_gsm','pf_gsm','pf_para']
    vars = []
    foreach pre0, pres do vars = [vars, pre0+tvar]
    
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar)
    
    tplot, vars, position=poss, trange=utr1
    
    sgclose


  ;-zoom in to pi2.
    ofn = shomedir()+'/fig_pflux_rbx_pi2.pdf'
    ;ofn = 0
    sgopen, ofn, xsize=6, ysize=8, /inch
    
    device, decomposed=0
    loadct2, 43
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    tvar = ['e_mgsm','b_gsm','pf_gsm','pf_para']
    vars = []
    foreach pre0, pres do vars = [vars, pre0+tvar]
    
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar)
    
    tplot, vars, position=poss, trange=utr3
    
    sgclose


end
