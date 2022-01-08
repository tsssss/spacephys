
; calc distribution function f, from time-energy spectrogram, 
; assuming isotropic distribution.

; pro test_spec_to_f

;    me = 0.91e-30   ; kg.
    cc = 3.235e-31  ; account for unit conversion, (me/e)^2*1e4*1e-12.

    potr = ['1998-09-14/18:30','1998-09-14/18:45']
    fatr = ['1998-09-14/18:30','1998-09-14/18:45']
    
    
    utr = time_double(potr)
    hyd = sread_polar_hydra(utr)
    
    ; polar hydra h0.
    ; JEe is the electron differential energy flux, in 1/(cm^2-s-sr).
    ; multiply 4pi to get the total omnidirectional differential energy flux, 
    ; in 1/(cm^2-s).
    ; fe is the distribution function.
    ; fe = JEe*0.5*(m/W)^2, where m is mass in kg, W is energy in eV.

    ets = hyd.epoch
    uts = sfmepoch(ets,'unix')
    nrec = n_elements(ets)

    jees = hyd.jee
    enes = hyd.ene  ; eV.    
    fes = jees
    
    jeis = hyd.jei
    enis = hyd.eni
    fis = jeis
    
    for i = 0, nrec-1 do fes[i,*] = jees[i,*]/enes^2*cc
    for i = 0, nrec-1 do fis[i,*] = jeis[i,*]/enis^2*cc
    
    tvar = 'po_fe'
    store_data, tvar, uts, fes, enes
    options, tvar, 'xlog', 0
    options, tvar, 'ylog', 1
    options, tvar, 'yrange', [10,2e4]
    options, tvar, 'ystyle', 1
    options, tvar, 'zlog', 1
    options, tvar, 'spec', 1
    options, tvar, 'no_interp', 1
    
    tvar = 'po_fi'
    store_data, tvar, uts, fis, enis
    options, tvar, 'xlog', 0
    options, tvar, 'ylog', 1
    options, tvar, 'yrange', [10,2e4]
    options, tvar, 'ystyle', 1
    options, tvar, 'zlog', 1
    options, tvar, 'spec', 1
    options, tvar, 'no_interp', 1

    
    
    utr = time_double(fatr)
    ees = sread_fast_esa(utr, type = 'ees')
    
    ets = ees.epoch
    uts = sfmepoch(ets,'unix')
    nrec = n_elements(ets)
    
    jees = ees.el_0
    enes = ees.el_en
    fes = jees
    
    for i = 0, nrec-1 do fes[i,*] = jees[i,*]/enes[i,*]^2*cc

    tvar = 'fa_fe_0'
    store_data, tvar, uts, fes, enes
    
    jees = ees.el_90
    enes = ees.el_en
    fes = jees
    
    for i = 0, nrec-1 do fes[i,*] = jees[i,*]/enes[i,*]^2*cc
    
    tvar = 'fa_fe_90'
    store_data, tvar, uts, fes, enes
    
    jees = ees.el_180
    enes = ees.el_en
    fes = jees
    
    for i = 0, nrec-1 do fes[i,*] = jees[i,*]/enes[i,*]^2*cc
    
    tvar = 'fa_fe_180'
    store_data, tvar, uts, fes, enes
    
    
    tvar = 'fa_fe_'+['0','90','180']
    options, tvar, 'xlog', 0
    options, tvar, 'ylog', 1
    options, tvar, 'yrange', [10,2e4]
    options, tvar, 'ystyle', 1
    options, tvar, 'zlog', 1
    options, tvar, 'spec', 1
    options, tvar, 'no_interp', 1
    
    
    
    utr = time_double(fatr)
    ies = sread_fast_esa(utr, type = 'ies')
    
    ets = ies.epoch
    uts = sfmepoch(ets,'unix')
    nrec = n_elements(ets)
    
    jeis = ies.ion_0
    enis = ies.ion_en
    fis = jeis
    
    for i = 0, nrec-1 do fis[i,*] = jeis[i,*]/enis[i,*]^2*cc
    
    tvar = 'fa_fi_0'
    store_data, tvar, uts, fis, enis
    
    jeis = ies.ion_90
    enis = ies.ion_en
    fis = jeis
    
    for i = 0, nrec-1 do fis[i,*] = jeis[i,*]/enis[i,*]^2*cc
    
    tvar = 'fa_fi_90'
    store_data, tvar, uts, fis, enis
    
    jeis = ies.ion_180
    enis = ies.ion_en
    fis = jeis
    
    for i = 0, nrec-1 do fis[i,*] = jeis[i,*]/enis[i,*]^2*cc
    
    tvar = 'fa_fi_180'
    store_data, tvar, uts, fis, enis
    
    
    tvar = 'fa_fi_'+['0','90','180']
    options, tvar, 'xlog', 0
    options, tvar, 'ylog', 1
    options, tvar, 'yrange', [10,2e4]
    options, tvar, 'ystyle', 1
    options, tvar, 'zlog', 1
    options, tvar, 'spec', 1
    options, tvar, 'no_interp', 1
    
    
    

    ofn = shomedir()+'/1998_0914_fe_compare.pdf'
    sgopen, ofn
    
    device, decomposed = 0
    loadct2, 43

    tpout = time_double('1998-09-14/19:00')
    
    get_data, 'po_fe', pouts, pofes, poens
    tmp = min(pouts-tpout, idx, /absolute)
    tpoef = pofes[idx,*]
    tpoen = poens
    
    titl = 'f!Iele!N vs energy. Polar:black, FAST:0-red, 90-green, 180-blue'
    plot, tpoen, tpoef, xlog = 1, ylog = 1, title = titl, $
        xtitle = 'Energy (eV)', ytitle = 'f (s!E3!Ncm!E-6!N)'
        
    tfaut = time_double('1998-09-14/18:33')
    eshift = 0
    
    get_data, 'fa_fe_0', fauts, fafes, faens
    tmp = min(fauts-tfaut, idx, /absolute)
    tfaef = fafes[idx,*]
    tfaen = faens[idx,*]
    oplot, tfaen-eshift, tfaef, color = 6
    
    get_data, 'fa_fe_90', fauts, fafes, faens
    tmp = min(fauts-tfaut, idx, /absolute)
    tfaef = fafes[idx,*]
    tfaen = faens[idx,*]
    oplot, tfaen-eshift, tfaef, color = 4
    
    get_data, 'fa_fe_180', fauts, fafes, faens
    tmp = min(fauts-tfaut, idx, /absolute)
    tfaef = fafes[idx,*]
    tfaen = faens[idx,*]
    oplot, tfaen-eshift, tfaef, color = 2
    
    sgclose
    
    
    ofn = shomedir()+'/1998_0914_fi_compare.pdf'
    sgopen, ofn
    
    device, decomposed = 0
    loadct2, 43
    
    tpout = time_double('1998-09-14/19:00')
    eshift = 0
    
    get_data, 'po_fi', pouts, pofes, poens
    tmp = min(pouts-tpout, idx, /absolute)
    tpoef = pofes[idx,*]
    tpoen = poens
    
    titl = 'f!Iion!N vs energy. Polar:black, FAST:0-red, 90-green, 180-blue'
    plot, tpoen, tpoef, xlog = 1, ylog = 1, title = titl, $
        xtitle = 'Energy (eV)', ytitle = 'f (s!E3!Ncm!E-6!N)'
        
    tfaut = time_double('1998-09-14/18:33')
    
    get_data, 'fa_fi_0', fauts, fafes, faens
    tmp = min(fauts-tfaut, idx, /absolute)
    tfaef = fafes[idx,*]
    tfaen = faens[idx,*]
    oplot, tfaen-eshift, tfaef, color = 6
    
    get_data, 'fa_fi_90', fauts, fafes, faens
    tmp = min(fauts-tfaut, idx, /absolute)
    tfaef = fafes[idx,*]
    tfaen = faens[idx,*]
    oplot, tfaen-eshift, tfaef, color = 4
    
    get_data, 'fa_fi_180', fauts, fafes, faens
    tmp = min(fauts-tfaut, idx, /absolute)
    tfaef = fafes[idx,*]
    tfaen = faens[idx,*]
    oplot, tfaen-eshift, tfaef, color = 2
    
    sgclose
end
