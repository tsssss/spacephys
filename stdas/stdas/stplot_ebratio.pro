
function stplot_ebratio, dename, dbname, pfname, trange = tr, $
    deidx = deidx, dbidx = dbidx, pfidx = pfidx, method = method

    if n_elements(deidx) eq 0 then deidx = 0
    if n_elements(dbidx) eq 0 then dbidx = 0
    if n_elements(pfidx) eq 0 then pfidx = 0
    
    get_data, dename, t0, de
    get_data, dbname, t0, db
    
    if n_elements(tr) eq 2 then begin
        idx = where(t0 ge min(tr) and t0 le max(tr), cnt)
        if cnt eq 0 then message, 'no data for given time range ...'
    endif

    if n_elements(idx) eq 0 then begin
        de = de[*,deidx]
        db = db[*,dbidx]
    endif else begin
        de = de[idx,deidx]
        db = db[idx,dbidx]
    endelse
    
    if n_elements(method) eq 0 then method = ''
    if method eq 'minmax' then begin
        a = max(de,/absolute)
        b = max(db,/absolute)
        r = abs(a/b*1e3)
        return, [r,a,b]
    endif
    
    if n_elements(pfname) eq 0 then begin
        b = regress(db, de, const = a, correlation = r)
        return, [a,b,r]
    endif else begin        ; dE and dB and dE/dB at max poynting flux.
        get_data, pfname, t0, pf
        if n_elements(idx) eq 0 then pf = pf[*,pfidx] else pf = pf[idx,pfidx]
        tmp = max(pf, idx)
        a = de[idx] & b = db[idx] & r = a/b*1e3
        return, [r,a,b]
    endelse
end