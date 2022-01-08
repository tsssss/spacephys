
    ;_2014_1222_load_data;, /reload
    tplot_options, 'ystyle', 1
    
    utr1 = time_double(['2014-12-22/04:00','2014-12-22/06:30'])
    utr3 = time_double(['2014-12-22/03:00','2014-12-22/08:00'])
    utr2 = time_double(['2014-12-22/05:25'])
    pres = ['rba','rbb','tha','the','thd']+'_'
    
    

;---Plot.
    ofn = shomedir()+'/fig_b_vs_density.pdf'
    ;ofn = 0
    sgopen, ofn, xsize=5, ysize=8, /inch
    
    device, decomposed=0
    loadct2, 43
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    
    foreach pre0, pres do begin
        get_data, pre0+'b_gsm', uts, bgsm
        if strmid(pre0,0,2) eq 'rb' then begin
            get_data, pre0+'t04s_b0_gsm', tuts, b0gsm
            b0gsm = sinterpol(b0gsm, tuts, uts)
            dbgsm = bgsm-b0gsm
        endif else begin
            dbgsm = bgsm
        endelse
        store_data, pre0+'db_gsm', uts, dbgsm
        store_data, pre0+'dbmag', uts, snorm(dbgsm), limits={ytitle:'(nT)', labels:strupcase(strmid(pre0,0,3))}
    endforeach
    
    foreach tprobe, ['a','b'] do begin
        pre0 = 'rb'+tprobe+'_'
        hope = sread_rbsp_hope_l3(utr1, probe=tprobe, type='mom')
        store_data, pre0+'density', sfmepoch(hope.(0),'unix'), hope.(1)
    endforeach
    
    tvar = ['rba','rbb']+'_dbmag'
    options, tvar, 'yrange', [0,60]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
        
    tvar = ['tha','thd','the']+'_dbmag'
    options, tvar, 'yrange', [0,120]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    
    options, '*_density', 'ylog', 1
    options, ['tha','thd','the']+'_density', 'yrange', [0.5,5]
    options, ['rba','rbb']+'_density', 'yrange', [0.8,8]
    
    tvar = ['rbb','rba','tha','the','thd']
    colors = [1,2,4,6,0]
    vars = []
    foreach tmp, tvar, i do begin
        vars = [vars, tmp+['_dbmag','_density']]
        options, tmp+['_dbmag','_density'], 'colors', colors[i]
    endforeach
    nvar = n_elements(vars)
    figlabs = []
    foreach tmp, ['a','b','c','d','e'] do figlabs = [figlabs,tmp+'-'+['1','2']+'.']

    
    poss = sgcalcpos(nvar, lmargin=10, tmargin=2, rmargin=8, bmargin=5)
    tplot, vars, trange=utr1, position=poss, /novtitle
    for i=0,nvar-1 do xyouts, poss[0,i]-xchsz*8, poss[3,i]-ychsz*0.8, /normal, figlabs[i]
        
    sgclose

end
