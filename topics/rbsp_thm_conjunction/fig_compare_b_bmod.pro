
    _2014_1222_load_data;, /reload
    tplot_options, 'ystyle', 1
    
    utr1 = time_double(['2014-12-22/04:00','2014-12-22/06:30'])
    utr3 = time_double(['2014-12-22/03:00','2014-12-22/08:00'])
    utr2 = time_double(['2014-12-22/05:25'])
    pres = ['rba','rbb','tha','the','thd']+'_'
    
    models = ['t89','t96','t01','t04s']
    nmodel = n_elements(models)
    colors = [0,1,2,4,6]
    

;---Plot.
    foreach pre0, pres do begin
        ofn = shomedir()+'/fig_compare_b_bmod_'+strmid(pre0,0,3)+'.pdf'
        sgopen, ofn, xsize=8.5, ysize=11, /inch
        
        device, decomposed=0
        loadct2, 43
                
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        get_data, pre0+'b_gsm', uts, bgsm
        idx = where(uts ge utr1[0] and uts le utr1[1], nrec)
        uts = uts[idx]
        bgsm = bgsm[idx,*]
        bx = dblarr(nrec,1+nmodel)
        by = dblarr(nrec,1+nmodel)
        bz = dblarr(nrec,1+nmodel)
        bx[*,0] = bgsm[*,0]
        by[*,0] = bgsm[*,1]
        bz[*,0] = bgsm[*,2]
        err = dblarr(nrec,nmodel)
        for i=0,nmodel-1 do begin
            get_data, pre0+models[i]+'_b0_gsm', tuts, b0gsm
            b0gsm = sinterpol(b0gsm,tuts, uts)
            bx[*,i+1] = b0gsm[*,0]
            by[*,i+1] = b0gsm[*,1]
            bz[*,i+1] = b0gsm[*,2]
            err[*,i] = sqrt(total((b0gsm-bgsm)^2,2))
        endfor

        store_data, pre0+'bx_comb', uts, bx, limits={ytitle:'(nT)', colors:colors, labels:['Bx',strupcase(models)]}
        store_data, pre0+'by_comb', uts, by, limits={ytitle:'(nT)', colors:colors, labels:['By',strupcase(models)]}
        store_data, pre0+'bz_comb', uts, bz, limits={ytitle:'(nT)', colors:colors, labels:['Bz',strupcase(models)]}
        store_data, pre0+'b_error', uts, err, limits={ytitle:'(nT)', colors:colors[1:*], labels:strupcase(models)}
        
        errs = strarr(nmodel)
        for i=0,nmodel-1 do begin
            errs[i] = sgnum2str(mean(err[*,i]),ndec=1)
        endfor
        options, pre0+'b_error', 'labels', 'dB:'+errs+'nT'
        
        vars = pre0+['b'+['x','y','z']+'_comb','b_error']
        nvar = n_elements(vars)
        poss = sgcalcpos(nvar)
        tplot, vars, trange=utr1, /novtitle, position=poss
        
        sgclose
    endforeach

end
