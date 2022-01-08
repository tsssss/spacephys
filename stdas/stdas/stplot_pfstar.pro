;+
; calc the correction coef for poynting flux.
; the input should be the E and B fields in FAC.
;-

function stplot_pfstar, dename, dbname, trange = tr

    get_data, dename, t0, de
    get_data, dbname, t0, db
    
    if n_elements(tr) eq 2 then begin
        idx = where(t0 ge min(tr) and t0 le max(tr), cnt)
        if cnt eq 0 then message, 'no data for given time range ...'
    endif
    
    if n_elements(idx) eq 0 then begin
        de = de[*,*]
        db = db[*,*]
    endif else begin
        de = de[idx,*]
        db = db[idx,*]
    endelse
    
    ; the components are [v,p,b].
    mdbp = max(db[*,1],/absolute)
    mdbv = max(db[*,0],/absolute)
    mdep = max(de[*,1],/absolute)
    mdev = max(de[*,0],/absolute)
    
;    dev = de[*,0]
;    dep = de[*,1]
;    dbv = db[*,0]
;    dbp = db[*,1]
;    pflux = dev*dbp-dep*dbv
    
    pflux = mdev*mdbp-mdep*mdbv
    mdepp = -mdbv/mdbp*mdev
    
    pfluxp = mdev*mdbp-mdepp*mdbv
    
    ratio = pfluxp/pflux
    if ratio gt 2 then ratio = 1
    
;    depp = -dbv/dbp*dev
;    pfluxp = dev*dbp-depp*dbv
    return, double(ratio)
    
end