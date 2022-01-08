
    ;_2014_1222_load_data;, /reload
    tplot_options, 'ystyle', 1
    
    utr1 = time_double(['2014-12-22/04:00','2014-12-22/06:30'])
    utr3 = time_double(['2014-12-22/03:00','2014-12-22/08:00'])
    utr2 = time_double(['2014-12-22/05:25'])
    pres = ['rba','rbb','tha','the','thd']+'_'
    
    

;---Plot.
    ofn = shomedir()+'/fig_db_diff.pdf'
    ;ofn = 0
    sgopen, ofn, xsize=4, ysize=5, /inch
    
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
    
    
    tvar = ['rba','rbb']+'_dbmag'
    options, tvar, 'yrange', [0,60]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
        
    tvar = ['tha','thd','the']+'_dbmag'
    options, tvar, 'yrange', [0,120]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    
    vars = ['rbb','rba','tha','the','thd']+'_dbmag'
    nvar = n_elements(vars)
    figlabs = ['a','b','c','d','e']+'.'

    
    poss = sgcalcpos(nvar, lmargin=8, tmargin=2, rmargin=5, bmargin=5)
    tplot, vars, trange=utr1, position=poss, /novtitle
    for i=0,nvar-1 do xyouts, poss[0,i]-xchsz*6, poss[3,i]-ychsz*0.5, /normal, figlabs[i]
    ;for i=0,nvar-1 do xyouts, poss[0,i]+xchsz*1, poss[3,i]-ychsz*1.2, /normal, '|B-B!Dmod!N|'
        
    sgclose


    stop
    device, decomposed=0
    loadct2, 43
;---Calculate time shift.
    get_data, 'tha_dbmag', uts, dat0
    get_data, 'the_dbmag', ut1s, dat1 & dat1 = sinterpol(dat1, ut1s, uts)
    get_data, 'thd_dbmag', ut2s, dat2 & dat2 = sinterpol(dat2, ut2s, uts)
    idx = where(uts ge utr1[0] and uts le utr1[1])
    dat0 = dat0[idx]
    dat1 = dat1[idx]
    dat2 = dat2[idx]
    dr0 = sdatarate(uts)
    
    dts = smkarthm(-1200,1200,1, 'dx')
    dns = dts/dr0
    ndt = n_elements(dts)
    corrs = dblarr(ndt)
    for i=0, ndt-1 do begin
        corrs[i] = c_correlate(dat0,dat1, dns[i])
    endfor
    plot, dts, corrs

    maxcorr = max(corrs, idx)
    print, maxcorr, dts[idx]
    plot, dat0
    oplot, shift(dat1, -dns[idx]), color=6

    
    dts = smkarthm(-1200,1200,1, 'dx')
    dns = dts/dr0
    ndt = n_elements(dts)
    corrs = dblarr(ndt)
    for i=0, ndt-1 do begin
        corrs[i] = c_correlate(dat1,dat2, dns[i])
    endfor
    plot, dts, corrs

    maxcorr = max(corrs, idx)
    print, maxcorr, dts[idx]
    plot, dat1
    oplot, shift(dat2, -dns[idx]), color=6
    
    

    get_data, 'rba_dbmag', uts, dat0
    get_data, 'rbb_dbmag', uts, dat1
    idx = where(uts ge utr1[0] and uts le utr1[1])
    dat0 = dat0[idx]
    dat1 = dat1[idx]
    dr0 = sdatarate(uts)
    dts = smkarthm(-1800,1800,1, 'dx')
    dns = dts/dr0
    ndt = n_elements(dts)
    corrs = dblarr(ndt)
    for i=0, ndt-1 do begin
        corrs[i] = c_correlate(dat0,dat1, dns[i])
    endfor
    plot, dts, corrs
    maxcorr = max(corrs, idx)
    print, maxcorr, dts[idx]
    plot, dat0
    oplot, shift(dat1, -dns[idx]), color=6
    
    max1 = max(dat0, idx0)
    max2 = max(dat1, idx1)
    tdt = uts[idx0]-uts[idx1]
    print, tdt
    plot, dat0
    oplot, shift(dat1, tdt/dr0), color=6
stop
end
