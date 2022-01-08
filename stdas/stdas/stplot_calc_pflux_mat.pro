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
;   bfilter, in, dblarr[m]/[m,2], opt. Filter in time for dB.
;   method, in, {'mat','das','gauss','morlet'}, opt. 'mat' is default.
;   tscale, in/out, dblarr[m], opt. Scales in time.
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
pro stplot_calc_pflux_mat, dename, dbname, pfname, $
    filter = tfilters, bfilter = btfilters, ids = iftrs, $
    tscale = tscales, scaleinfo = scaleinfo, $
    labels = labels, ytitle = ytitle, colors = colors

; [de,db,pf]name -> [de,db,pf]name_comp[1,2,3].
    ; [de,db,pf]name: varname for 3-D dE, dB, S.
    ; [de,db,pf]fvars: varname for 1-D components.
    if n_elements(dename) eq 0 then dename = 'de_fac'
    if n_elements(dbname) eq 0 then dbname = 'db_fac'
    if n_elements(pfname) eq 0 then pfname = 'pf_fac'
    stplot_split, dename, newname = defvars ; get defvars.
    stplot_split, dbname, newname = dbfvars ; get dbfvars.
    ndim = n_elements(defvars)
    
    ; init datarate, S's ytitle, labels, colors.
    get_data, dename, t0, limits = lm
    dr = sdatarate(t0) & dr1 = 1d/dr    ; data rate.
    nrec = n_elements(t0)               ; # of records.

    if n_elements(ytitle) eq 0 then ytitle = 'S (mW/m!U2!N)'

    if n_elements(labels) ne ndim then $
        labels = stagexist('labels',lm)? lm.labels: strarr(ndim)

    if n_elements(colors) ne ndim then begin
        colors = stagexist('colors',lm)? lm.colors: indgen(ndim)
        if ndim eq 3 then colors = [6,4,2]  ; rgb in color table 43.
    endif
        
    nofilter = (n_elements(tfilters) eq 0) and (n_elements(btfilters) eq 0)
    noscale = (n_elements(tscales) eq 0) and (n_elements(scaleinfo) eq 0)
    
    
    ; calc poynting flux directly.
    if nofilter and noscale then begin
        get_data, dename, t0, de
        get_data, dbname, t0, db
        fp = spoynt(de,db)
        store_data, pfname, t0, fp, $
            limits = {ytitle:ytitle,labels:labels,colors:colors}
        return
    endif
    
    
    ; prepare scale.
    ; interface: {out: tscales,rscales,nscale}
    if n_elements(tscales) eq 0 then begin          ; no scales.
        if n_elements(scaleinfo) eq 0 then begin    ; no scale info.
            nscale = 60 & mins = 4 & maxs = 0.5*nrec
        endif else begin                            ; have scale info.
            tmp = n_elements(scaleinfo) & nscale = scaleinfo[tmp-1]
            if tmp eq 1 then begin                  ; only have nscale.
                mins = 4 & maxs = 0.5*nrec
            endif else begin                        ; have all info
                mins = scaleinfo[0] & maxs = scaleinfo[1]
            endelse
        endelse
        nscale = keyword_set(nscale)? nscale: 60
        tscales = smkgmtrc(mins,maxs,nscale,'n')
        rscales = round(tscales*dr1)
    endif else rscales = tscales*dr1
    rscales = double(rscales[uniq(rscales)])        ; scales in record.
    tscales = rscales*dr                            ; scales in time.
    nscale = n_elements(tscales)                    ; # of scales.
    
    ; prepare filter in time.
    ; interface: {out: nfilter, etftrs, erftrs, eiftrs, esftrs, $
    ;   btftrs, brftrs, biftrs, bsftrs}.
    
    if keyword_set(tfilters) then begin     ; common filter.
        etftrs = tfilters & btftrs = tfilters
    endif
    ; b field has different filter.
    if keyword_set(btfilters) then btftrs = btfilters
    nfilter = n_elements(stofilter(etftrs))/2   ; # of filter dE = dB.
    if n_elements(iftrs) eq 0 then iftrs = $    ; filter id string.
        string(indgen(nfilter)+1, format='(I0)')
    ; efilter.
    erftrs = etftrs*dr1             ; filter in record.
    erftrs = erftrs > 0
    erftrs = erftrs < nrec
    erftrs = erftrs[sort(erftrs)]
    erftrs = erftrs[uniq(erftrs)]
    erftrs = stofilter(erftrs)      ; rfilter in [n,2].
    etftrs = erftrs*dr
    esftrs = strarr(nfilter)        ; label for each band.
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
    bsftrs = strarr(nfilter)        ; label for each band.
    for i = 0, nfilter-1 do $
        bsftrs[i] = sgnum2str(min(btftrs[i,*]),msgn=3)+'-'+$
        sgnum2str(max(btftrs[i,*]),msgn=3)+'s'
    ; # of filters.
    nfilter = n_elements(iftrs)
    

; **** for each field, do mat spectrogram then filter and delete spectrogram.
; [de,db]name_comp[1,2,3] -> [de,db]name_comp[1,2,3]_mat[1,2,3...].
    fvars = defvars
    for i = 0, n_elements(fvars)-1 do begin
        tvar = fvars[i]+'_mat'
        stplot_mat, fvars[i], scale = tscales, newname = tvar
        for j = 0, nfilter-1 do $
            stplot_filter, tvar, 'mat', filter = etftrs[j,*], ifilter = iftrs[j]
        store_data, tvar, /delete
    endfor

    fvars = dbfvars
    for i = 0, n_elements(fvars)-1 do begin
        tvar = fvars[i]+'_mat'
        stplot_mat, fvars[i], scale = tscales, newname = tvar
        for j = 0, nfilter-1 do $
            stplot_filter, tvar, 'mat', filter = etftrs[j,*], ifilter = iftrs[j]
        store_data, tvar, /delete
    endfor

    if nofilter then return


; **** combine 1-D components into 3-D vector, delete the components.
; [de,db]name_comp[1,2,3]_mat[1,2,3...] -> [de,db]name_mat[1,2,3...].
    for i = 0, nfilter-1 do begin
        tmp = '_mat'+iftrs[i]
        stplot_merge, defvars+tmp, newname = dename+tmp, /delete
        stplot_merge, dbfvars+tmp, newname = dbname+tmp, /delete
    endfor
    idx = where(defvars eq dename, cnt)
    if cnt eq 0 then store_data, defvars, /delete
    idx = where(defvars eq dename, cnt)
    if cnt eq 0 then store_data, dbfvars, /delete
    
; **** calc poynting flux for each band.
    denames = dename+'_mat'+iftrs
    dbnames = dbname+'_mat'+iftrs
    pfnames = pfname+'_mat'+iftrs
    for i = 0, nfilter-1 do begin
        get_data, denames[i], t0, de
        get_data, dbnames[i], t0, db
        pf = spoynt(de,db)
        store_data, pfnames[i], t0, pf, $       ; use dE filter's label.
            limits = {ytitle:ytitle,colors:colors}
        tvars = [pfnames[i],denames[i],dbnames[i]]
        options, tvars, 'labels', strarr(ndim)+esftrs[i]
    endfor
    
; **** sum up db, de, poynting flux, for selected freq bands only.
; [de,db,pf]name_mat.
    stplot_total, dename+'_mat'+iftrs, newname = dename+'_mat', colors = colors
    stplot_total, dbname+'_mat'+iftrs, newname = dbname+'_mat', colors = colors
    stplot_total, pfname+'_mat'+iftrs, newname = pfname+'_mat', colors = colors
    
    options, dename+'_mat'+iftrs, 'ytitle', 'dE!C(mV/m)'
    options, dbname+'_mat'+iftrs, 'ytitle', 'dB!C(nT)'
    
    vars = [dename,dbname,pfname]+'_mat'
    options, vars, 'labels', sgnum2str(min(tfilters))+'-'+ $
        sgnum2str(max(tfilters))+'s'+['','','']
    options, dename+'_mat', 'ytitle', 'dE!C(mV/m)'
    options, dbname+'_mat', 'ytitle', 'dB!C(nT)'
end
