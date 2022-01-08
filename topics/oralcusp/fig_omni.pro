
pro fig_omni, ifn, tr, filename = ofn

    if n_elements(ofn) eq 0 then message, 'no output filename ...'

    cdf = scdfread(ifn)
    t0 = sfmepoch(*cdf[0].value, 'unix')

    ; IMF B GSM.
    store_data, 'bgsm', t0, [[*cdf[1].value],[*cdf[4].value],[*cdf[5].value]], $
        limits = {ytitle:'B (nT)',ystyle:1,yrange:[-25,25], $
        labels:['Bx GSM','By GSM','Bz GSM'],colors:[6,4,2],yminor:2}

    ; IMF V GSE.
    store_data, 'vgse', t0, [[*cdf[6].value],[*cdf[7].value],[*cdf[8].value]], $
        limits = {ytitle:'V (km/s)',ystyle:1,labels:['Vx GSE','Vy','Vz'],colors:[6,4,2]}
    
    ; IMF V GSM.
    cotrans, 'vgse', 'vgsm', /gse2gsm
    stplot_split, 'vgsm', newnames = ['vx_gsm','vy_gsm','vz_gsm'], $
        ytitles = ['Vx','Vy','Vz']+' (km/s)', labels = ['Vx','Vy','Vz']+' GSM'
    ylim, 'vx_gsm', -900, -300
    ylim, 'vy_gsm', -150, 100
    ylim, 'vz_gsm', -100, 100
    options, 'vx_gsm', 'yticks', 3
    options, 'vx_gsm', 'yminor', 4
    
    ; dynamic pressure.
    pdyn = *cdf[10].value
    store_data, 'p', t0, pdyn, limits = {yrange:[0,18],ytitle:'P (nPa)',ystyle:1}
    
    ; index.
    ae = *cdf[11].value
    symh = *cdf[14].value
    store_data, 'ae', t0, ae, limits = {yrange:[0,2400],ytitle:'AE (nT)',ystyle:1}
    store_data, 'symh', t0, symh, limits = {yrange:[-250,50],ytitle:'Sym/H (nT)',$
        yticks:2,yminor:5,ystyle:1}
    
    ; remove fill value.
    vars = ['vgse','vgsm','bgsm','vx_gsm','vy_gsm','vz_gsm','p']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], t0, tmp
        idx = where(tmp eq max(tmp,/absolute), cnt)
        if cnt ne 0 then begin
            tmp[idx] = !values.d_nan
            store_data, vars[i], t0, tmp
        endif
    endfor
    
    ; figure.
    options, 'vx_gsm', 'labels', 'Vx GSM'
    options, 'symh', 'labels', 'Sym/H'

    ; default settings.
    device, decomposed = 0 & loadct2, 43
    !p.font = 1 & !p.thick = 2
    !y.charsize = 0.8 & !x.charsize = 0.8
    time_stamp, /off
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 5
    tplot_options, 'labflag', 1
    
    vars = ['bgsm','vx_gsm','symh']
    if n_elements(tr) eq 0 then tlimit, /full
    sgwindow, 0, xsize = 4.5, ysize = 3, /inch
    sgpsopen, ofn
    sgindexcolor
    tplot, vars, trange = tr
    timebar, time_double('1998-09-25/04:27:55'), color = 6
    timebar, time_double('1998-09-25/04:29:40'), color = 6
    timebar, time_double('1998-09-25/05:27'), color = 2
    timebar, time_double('1998-09-25/05:40'), color = 2
    sgpsclose
    wdelete, 0
;    pstplot, filename = ofn, xsize = 5, ysize = 3, charsize = 0.3
end

ifn = sdiskdir('Works')+'/agu/agu2013/data/omni_hros_1min_19980925.cdf'
;ofn = shomedir()+'/19980925_omni_fullstorm.ps'
;fig_omni, ifn, filename = ofn
ofn = shomedir()+'/omni.eps'
tr = time_double(['1998-09-25/04:00','1998-09-25/06:00'])
fig_omni, ifn, tr, filename = ofn
end
