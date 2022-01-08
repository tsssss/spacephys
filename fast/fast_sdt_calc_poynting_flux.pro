;+
; Type: procedure.
; Purpose: FAST, sdt interface, calculate poynting flux and related variables.
;   Interface: {in: 'bmod_gei','btotal','db_facv','de_facv', $
;   'fa_pos','fa_vel','ilat','mlt','alt','dis'}
; Parameters: none.
; Keywords: none.
; Notes: This code should be run after fast_sdt_prep_poynting_flux.
; Dependence: tdas,slib.
; Author: Sheng Tian.
; History: 2013-11-21, Sheng Tian, create.
;-
pro fast_sdt_calc_poynting_flux, filter = tfilters, method = method, $
    scale = tscales, nscale = nscale, $
    pfname = pfname, labels = labels, ytitle = ytitle, colors = colors

    if n_elements(dbname) eq 0 then dbname = 'db_facv'
    if n_elements(dename) eq 0 then dename = 'de_facv'
    if n_elements(pfname) eq 0 then pfname = 'pf_facv'
    stplot_split, dbname
    stplot_split, dename
    ndim = n_elements(tnames(dbname+'_?'))
    idims = string(indgen(ndim)+1, format='(I0)')   ; dimension id string.
    fvars = [tnames(dbname+'_comp?'),tnames(dename+'_comp?')]
    
    get_data, dbname, t0, limits = lm
    dr = sdatarate(t0) & dr1 = 1D/dr    ; data rate.
    nrec = n_elements(t0)               ; # of records.
    if n_elements(ytitle) eq 0 then ytitle = 'S (mW/m!U2!N)'
    if n_elements(labels) ne 3 then $
        labels = stagexist('labels',lm)? lm.labels: strarr(3)
    if n_elements(colors) ne 3 then $
        colors = stagexist('colors',lm)? lm.colors: [6,4,2]
        
    ; calc poynting flux directly.
    if n_elements(tfilters) eq 0 then begin
        get_data, dbname, t0, db
        get_data, dename, t0, de
        fp = spoynt(de, db)
        store_data, pfname, t0, fp, $
            limits = {ytitle:ytitle,labels:labels,colors:colors}
        return
    endif
    
    ; filter in time.
    ; interface: {out: tfilters, rfilters, ifilters, sfilters, nfilter}.
    tfilters = stofilter(tfilters)      ; tfilter in [n,2].
    rfilters = tfilters*dr1             ; filter in record.
    nfilter = n_elements(tfilters)/2
    ifilters = string(indgen(nfilter)+1, format='(I0)')     ; filter id string.
    sfilters = strarr(nfilter)
    for i = 0, nfilter-1 do $
        sfilters[i] = snum2str(tfilters[i,0])+'-'+snum2str(tfilters[i,1])+'s'
    
    ; check method to filter, and prepare. interface: {out: fvars}.
    method = keyword_set(method)? method: 'mat'
    if method eq 'mat' then begin
        ; interface: {out: tscales,rscales,nscale}
        if n_elements(tscales) eq 0 then begin
            nscale = keyword_set(nscale)? nscale: 60
            rscales = round(smkgmtrc(4,0.5*nrec,nscale,'n'))
            rscales = double(rscales[uniq(rscales)])
            tscales = rscales*dr
        endif else rscales = tscales*dr1
        nscale = n_elements(tscales)
        ; mat.
        for i = 0, n_elements(fvars)-1 do $ ; 'x_compi'->'x_compi_mat'.
            stplot_mat, fvars[i], scale = tscales
        fvars = fvars+'_mat'
    endif else if method eq 'gauss' then begin
        ; guass cwt.
        fvars = fvars+'_gauss'
    endif else if method eq 'morlet' then begin
        ; morlet cwt.
        fvars = fvars+'_morlet'
    endif
    
    ; do filter.
    for i = 0, n_elements(fvars)-1 do $
        for j = 0, nfilter-1 do $
            stplot_filter, fvars[i], method, filter = tfilters[j,*], $
                ifilter = ifilters[j]

    ; combine fields.
    fvars = [dbname+'_comp?_'+method+'f'+ifilters, $
        dename+'_comp?_'+method+'f'+ifilters]
    newnames = [dbname+'_'+method+'f'+ifilters, $
        dename+'_'+method+'f'+ifilters]
    for i = 0, n_elements(fvars)-1 do $
        stplot_blend, tnames(fvars[i]), newname = newnames[i]
    
    ; calc poynting flux for each band.
    denames = [dename+'_'+method+'f'+ifilters]
    dbnames = [dbname+'_'+method+'f'+ifilters]
    pfnames = [pfname+'_'+method+'f'+ifilters]
    for i = 0, nfilter-1 do begin
        get_data, denames[i], t0, de
        get_data, dbnames[i], t0, db
        pf = spoynt(de, db)
        store_data, pfnames[i], t0, pf, $
            limits = {ytitle:ytitle,labels:labels+sfilters[i],colors:colors}
    endfor
    
    ; sum up db, de, poynting flux.
    get_data, dename, limits = lm
    stplot_renew, dename, newname = dename+'0'
    stplot_total, tnames(dename+'_'+method+'f?'), newname = dename, $
        ytitle = lm.ytitle, labels = lm.labels, colors = lm.colors
    get_data, dbname, limits = lm
    stplot_renew, dbname, newname = dbname+'0'
    stplot_total, tnames(dbname+'_'+method+'f?'), newname = dbname, $
        ytitle = lm.ytitle, labels = lm.labels, colors = lm.colors
    stplot_total, tnames(pfname+'_'+method+'f?'), newname = pfname, $
        ytitle = ytitle, labels = labels, colors = colors
end