;+
; get the slow sweep times based on the variations in ibias, guard and usher.
; vars: [ibias,guard,usher] tplot varname.
;-

pro bias_param_get_sdt_times, vars, sdtutr

    get_data, vars[0], t0, ibias
    get_data, vars[1], t0, guard
    get_data, vars[2], t0, usher

    sdtutr = [0d,0]
    spinrate = 12   ; sec.
    dibias = 1      ; nA.
    dusher = 1      ; volt.
    dguard = 1      ; volt.
    
    ; the largest count corresponds to the default value.
    ; the second largest count > value in eclipse.
    ; the smallest count may due to abonormal value.
    vs = round(guard)
    dval = dguard
    vals = vs[uniq(vs,sort(vs))]    ; uniq values and sorted.
    nval = n_elements(vals)         ; # of uniq values.
    for i = 0, nval-1 do begin
        idx = where(abs(vs-vals[i]) le dval)
        vals[i] = round(median(vs[idx]))  ; remove irregular value.
    endfor
    vals = vals[uniq(vals,sort(vals))]
    nval = n_elements(vals)
    hist = intarr(nval)             ; count.
    for i = 0, nval-1 do begin
        tmp = where(vs eq vals[i], cnt)
        hist[i] = cnt
    endfor
    idx = sort(hist)
    hist = hist[idx]
    vals = vals[idx]
    if n_elements(vals) le 2 then return    ; no sdt, can have 2 values (eclips/normal).
    guards = vs
    gval = vals[1]
    
    vs = round(usher)
    dval = dusher
    vals = vs[uniq(vs,sort(vs))]    ; uniq values and sorted.
    nval = n_elements(vals)         ; # of uniq values.
    for i = 0, nval-1 do begin
        idx = where(abs(vs-vals[i]) le dval)
        vals[i] = round(median(vs[idx]))  ; remove irregular value.
    endfor
    vals = vals[uniq(vals,sort(vals))]
    nval = n_elements(vals)
    hist = intarr(nval)             ; count.
    for i = 0, nval-1 do begin
        tmp = where(vs eq vals[i], cnt)
        hist[i] = cnt
    endfor
    idx = sort(hist)
    hist = hist[idx]
    vals = vals[idx]
    if n_elements(vals) le 2 then return    ; no sdt for this probe.
    ushers = vs
    uval = vals[1]
    
    ; sdt begins with change in [ibias,guard,usher] values.
    idx = where(guards eq gval and ushers eq uval)
    minibias = min(ibias[idx])
    idx = (where(ibias eq minibias))[0]
    idx = idx[0]+[-1,1]*1800/spinrate   ; locate the rough time for sdt.
    idx[0]>= 0 & idx[1]<= n_elements(t0)-1
    tmp = where(abs(ibias[idx[0]:idx[1]]-minibias) le dibias)
    sdtutr[0] = t0[idx[0]+tmp[0]]
    idx = idx[0]+tmp[0]
    
    ; sdt ends with change in datarate.
    for i = idx[0], n_elements(t0)-2 do $
        if t0[i+1]-t0[i] gt spinrate then break

    sdtutr[1] = t0[i]

end