;+
; set common yrange (and ylog) if minmax is given.
; get common yrange if minmax is omitted.
; set set to change yrange directly.
; set newname to leave vars0 untouched, and easy to delete the changed vars.
;-
pro stplot_minmax, vars0, minmax, log = log, set = set, get = get, newname = tnames, trange = tr

    vars = tnames(vars0)
    nvar = n_elements(vars)
    if nvar eq 0 then return
    
    if n_elements(tnames) ne 0 then begin
        for i = 0, nvar-1 do stplot_renew, vars[i], newname = tnames[i]
        vars = tnames
    endif
    
    if keyword_set(get) then tmp = temporary(minmax)
    
    if n_elements(minmax) eq 0 then begin
        minmax = []
        for i = 0, nvar-1 do begin
            get_data, vars[i], tmp, val
            if n_elements(tr) eq 2 then val = val[where(tmp ge min(tr) and tmp le max(tr))]
            tminmax = min(val, /nan, max = tmp) & tminmax = [tminmax,tmp]
            minmax = min([minmax,tminmax], max = tmp) & minmax = [minmax,tmp]
        endfor
    endif
    
    if ~keyword_set(set) and n_elements(minmax) ne 2 then return
    
    tminmax = min(minmax, max = tmp) & tminmax = [tminmax,tmp]
    if keyword_set(log) then log = 1 else log = 0
    options, vars, 'yrange', tminmax
    options, vars, 'ylog', log
    options, vars, 'ystyle', 1
end
