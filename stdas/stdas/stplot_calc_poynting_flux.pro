;+
; Type: procedure.
; Purpose: Calculate poynting flux S and related variables in tplot.
;   Interface: {in: dename, dbname}{optional in: filter, scales}
;   {out:pfname}{bands of de, db, pf}.
; Parameters:
;   dename, in, string, req. 3-d dE in <coord>.
;   dbname, in, string, req. 3-d dB in <coord>.
;   pfname, in, string, req. 3-d S in <coord>.
; Keywords:
;   filter, in, dblarr[m]/[m,2], opt. Filter in time for dE and dB.
;   efilter, in, dblarr[m]/[m,2], opt. Filter in time for dE.
;   bfilter, in, dblarr[m]/[m,2], opt. Filter in time for dB.
;   method, in, {'mat','das','gauss','morlet'}, opt. 'mat' is default.
;   scale, in, dblarr[m], opt. Scales in time.
;   scaleinfo, in, intarr[3]/integer, opt. [min,max,ns] or ns in time.
; Notes: Supply 3-d dE and dB in <coord> to calculate S in <coord>. For example 
;   dename = 'de_facv', dbname = 'db_facv', pfname = 'pf_facv'.
;       Set common filter for dE and dB using filter, otherwise set efilter and
;   bfilter separately. Number of filters for efilter and bfilter MUST be
;   equal. If no filter is set, then use dE and dB directly, get S and done. 
;   Otherwise, first decompose field to components, e.g., 'de_facv' produces
;   'de_facv_comp1','de_facv_comp2','de_facv_comp3', similar for 'db_facv'.
;   Then filter each component using one of the methods: (1) 'das' detrend and 
;   smooth; (2) 'mat' moving average transform; (3) 'gauss' Gauss wavelet; 
;   (4) 'morlet' Morlet wavelet.
;       For methods other than 'das', time-period spectrogram are calculated.
;   User can specify scales in time or supply scale info in (i) [minscale, 
;   maxscale, nscale] in time, or (ii) # of scale. In the latter case,
;   minscale is 4 and maxscale is 0.5*nrec. Calculate spectrogram for each 
;   component, e.g., 'de_facv_comp1' gives 'de_facv_comp1_mat' for 'mat'.
;       For all methods, given m pairs of filters in time, m bands are 
;   calculated for each component, e.g., for 'mat' and m = 4, 'de_facv_comp1' 
;   has bands 'de_facv_comp1_mat1', ..., 'de_facv_comp1_matf4'. Then combine
;   the 3 components at each bands to 'de_facv_matf1', ..., 'de_facv_dasf4' and
;   delete 'de_facv_comp1_matf1', ..., 'de_facv_comp1_matf4'.
;       For each band, calculate S, e.g., use 'de_facv_matf1' and 
;   'db_facv_matf1' to get 'pf_facv_matf1'.
;       Finally, get total filtered dE and dB and total S: 'de_facv_mat', 
;       'db_facv_mat', 'pf_facv_mat'.
; Dependence: tdas,slib.
; History:
;   2013-11-21, Sheng Tian, create.
;-
pro stplot_calc_poynting_flux, dename, dbname, pfname, method = method, $
    filter = tfilters, efilter = etftrs, bfilter = btftrs, idfilter = iftrs, $
    scale = tscales, scaleinfo = scaleinfo, $
    labels = labels, ytitle = ytitle, colors = colors

    ; [de,db,pf]name: varname for 3-D dE, dB, S.
    ; [de,db,pf]fvars: varname for 1-D components.
    if n_elements(dename) eq 0 then dename = 'de_fac'
    if n_elements(dbname) eq 0 then dbname = 'db_fac'
    if n_elements(pfname) eq 0 then pfname = 'pf_fac'
    stplot_split, dename, newname = defvars ; get defvars.
    stplot_split, dbname, newname = dbfvars ; get dbfvars.
    fvars = [defvars,dbfvars]
    ndim = n_elements(defvars)
    
    ; init datarate, and S's ytitle, labels, colors, and method.
    get_data, dename, t0, limits = lm
    dr = sdatarate(t0) & dr1 = 1D/dr    ; data rate.
    nrec = n_elements(t0)               ; # of records.
    if n_elements(ytitle) eq 0 then ytitle = 'S (mW/m!U2!N)'
    if n_elements(labels) ne 3 then $
        labels = stagexist('labels',lm)? lm.labels: strarr(3)
    if n_elements(colors) ne 3 then $
        colors = stagexist('colors',lm)? lm.colors: [6,4,2]
    method = keyword_set(method)? method: 'mat'
        
    ; calc poynting flux directly.
    if n_elements(tfilters) eq 0 and n_elements(etftrs) eq 0 and $
    n_elements(btftrs) eq 0 then begin
        get_data, dename, t0, de
        get_data, dbname, t0, db
        fp = spoynt(de,db)
        store_data, pfname, t0, fp, $
            limits = {ytitle:ytitle,labels:labels,colors:colors}
        return
    endif
    
    ; prepare filter in time.
    ; interface: {out: nfilter, etftrs, erftrs, eiftrs, esftrs, $
    ;   btftrs, brftrs, biftrs, bsftrs}.
    if keyword_set(tfilters) then begin     ; common filter.
        etftrs = tfilters & btftrs = tfilters
    endif   ; else assume etftrs and btftrs are set both.
    nfilter = n_elements(stofilter(etftrs))/2        ; # of filter dE = dB.
    if n_elements(iftrs) eq 0 then iftrs = $      ; filter id string.
        'f'+string(indgen(nfilter)+1, format='(I0)')
    ; efilter.
    erftrs = etftrs*dr1             ; filter in record.
    erftrs = erftrs > 0
    erftrs = erftrs < nrec
    erftrs = erftrs[sort(erftrs)]
    erftrs = erftrs[uniq(erftrs)]
    erftrs = stofilter(erftrs)      ; rfilter in [n,2].
    etftrs = erftrs*dr
    esftrs = strarr(nfilter)
    for i = 0, nfilter-1 do $
        esftrs[i] = sgnum2str(min(etftrs[i,*]),msgn=3)+'-'+$
            sgnum2str(max(etftrs[i,*]),msgn=3)+'s'
    ; bfilter.
    brftrs = btftrs*dr1             ; filter in record.
    brftrs = brftrs > 0
    brftrs = brftrs < nrec
    brftrs = brftrs[sort(brftrs)]
    brftrs = brftrs[uniq(brftrs)]
    brftrs = stofilter(brftrs)      ; rfilter in [n,2].
    btftrs = brftrs*dr
    bsftrs = strarr(nfilter)
    for i = 0, nfilter-1 do $
        bsftrs[i] = sgnum2str(min(btftrs[i,*]),msgn=3)+'-'+$
            sgnum2str(max(btftrs[i,*]),msgn=3)+'s'
    nftr = n_elements(iftrs)
    
    ; check method to filter, and prepare. interface: {out: defvars,dbfvars}.
    if method ne 'das' then begin
        ; interface: {out: tscales,rscales,nscale}
        if n_elements(tscales) eq 0 then begin          ; no scales.
            if n_elements(scaleinfo) eq 0 then begin    ; no scale info.
                nscale = 60 & mins = 4 & maxs = 0.5*nrec
            endif else begin        ; have scale info.
                tmp = n_elements(scaleinfo) & nscale = scaleinfo[tmp-1]
                if tmp eq 1 then begin ; only have nscale.
                    mins = 4 & maxs = 0.5*nrec
                endif else begin    ; have all info
                    mins = scaleinfo[0] & maxs = scaleinfo[1]
                endelse
            endelse
            nscale = keyword_set(nscale)? nscale: 60
            tscales = smkgmtrc(mins,maxs,nscale,'n')
            rscales = round(tscales*dr1)
        endif else rscales = tscales*dr1
        rscales = double(rscales[uniq(rscales)])
        tscales = rscales*dr
        nscale = n_elements(tscales)
    endif
    
    if method eq 'mat' then begin
        for i = 0, n_elements(fvars)-1 do $ ; 'x_compi'->'x_compi_mat'.
            stplot_mat, fvars[i], scale = tscales
        fvars += '_mat' & desvars = defvars+'_mat' & dbsvars = dbfvars+'_mat'
    endif else if method eq 'gauss' then begin
        ; guass cwt.
        fvars += '_gauss' & desvars = defvars+'_gauss' & dbsvars = dbfvars+'_gauss'
    endif else if method eq 'morlet' then begin
        ; morlet cwt.
        fvars += '_morlet' & desvars = defvars+'_morlet' & dbsvars = dbfvars+'_morlet'
    endif
    
    ; do filter for low, selected, and high freq bands.
    for i = 0, n_elements(defvars)-1 do $
        for j = 0, nftr-1 do $
            stplot_filter, desvars[i], method, filter = etftrs[j,*], ifilter = iftrs[j]
    for i = 0, n_elements(dbfvars)-1 do $
        for j = 0, nftr-1 do $
            stplot_filter, dbsvars[i], method, filter = btftrs[j,*], ifilter = iftrs[j]

    ; combine 1-D components into 3-D vector, delete the components.
    for i = 0, nftr-1 do begin
        tmp = '_'+method+iftrs[i]
        stplot_merge, defvars+tmp, newname = dename+tmp, /delete
        stplot_merge, dbfvars+tmp, newname = dbname+tmp, /delete
    endfor
    idx = where(defvars eq dename, cnt)
    if cnt eq 0 then store_data, defvars, /delete
    idx = where(defvars eq dename, cnt)
    if cnt eq 0 then store_data, dbfvars, /delete
    
    ; calc poynting flux for each band.
    denames = dename+'_'+method+iftrs
    dbnames = dbname+'_'+method+iftrs
    pfnames = pfname+'_'+method+iftrs
    for i = 0, nftr-1 do begin
        get_data, denames[i], t0, de
        get_data, dbnames[i], t0, db
        pf = spoynt(de,db)
        store_data, pfnames[i], t0, pf, $       ; use dE filter's label.
            limits = {ytitle:ytitle,colors:colors}
        tvars = [pfnames[i],denames[i],dbnames[i]]
        options, tvars, 'labels', strarr(ndim)+esftrs[i]
    endfor
    
    ; sum up db, de, poynting flux, for selected freq bands only.
    get_data, dename, limits = lm
    stplot_total, dename+'_'+method+iftrs, newname = dename+'_'+method, $
        labels = lm.labels, colors = lm.colors
    get_data, dbname, limits = lm
    stplot_total, dbname+'_'+method+iftrs, newname = dbname+'_'+method, $
        labels = lm.labels, colors = lm.colors
    stplot_total, pfname+'_'+method+iftrs, newname = pfname+'_'+method, $
        labels = labels, colors = colors
end
