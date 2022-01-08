;+
; Preprocess B field in GSM: separate model B and dB, calc |B|.
;
; If posvar is provided, then use it to calc bmod in T89. T89 is used b/c it
; provides a very smooth background, which does not introduce extra structure
; in dB. Otherwise use smoothing to get a bmod.
;
; bvar and posvar should be on the same uniform time.
;-
pro stplot_prep_bfield, bvar, posvar=posvar, addto=varlist

    rgb = [6,4,2]
    xyz = ['x','y','z']
    
    pre0 = strmid(bvar,0,strpos(bvar,'_'))+'_'  ; '' if no _.
    get_data, bvar, uts, bgsm
    nrec = n_elements(uts)

    ; get bmod.
    bmodgsm = dblarr(nrec,3)
    if keyword_set(posvar) then begin
        ; use pos to get T89 bmod.
        get_data, posvar, uts, rgsm
        ets = stoepoch(uts, 'unix')

        par = 2d
        for i=0, nrec-1 do begin
            tet = ets[i]
            geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
            geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date

            x0 = rgsm[i,0]
            y0 = rgsm[i,1]
            z0 = rgsm[i,2]
            geopack_igrf_gsm, x0, y0, z0, bx0, by0, bz0
            geopack_t89, par, x0, y0, z0, bx1, by1, bz1

            bmodgsm[i,*] = [bx0,by0,bz0]+[bx1,by1,bz1]
        endfor
    endif else begin
        for i=0,2 do bmodgsm[*,i] = scalcbg(bgsm[*,i])
    endelse
    store_data, pre0+'bmod_gsm', uts, bmodgsm, limits=$
        {ytitle:'(nT)', colors:rgb, labels:'GSM Bmod'+xyz, labflag:-1}

    ; get dbgsm and further refine if smooth is set.
    dbgsm = bgsm-bmodgsm
    store_data, pre0+'db_gsm', uts, dbgsm, limits=$
        {ytitle:'(nT)', colors:rgb, labels:'GSM dB'+xyz, labflag:-1}

    ; calc |B|.
    store_data, pre0+'bmag', uts, snorm(bgsm), limits=$
        {ytitle:'(nT)', labels:'|B|', ynozero:1, labflag:-1}
        
    myvars = pre0+['bmod_gsm','db_gsm','bmag']
    if not keyword_set(varlist) then varlist=[]
    varlist = [varlist, myvars]
    
end
