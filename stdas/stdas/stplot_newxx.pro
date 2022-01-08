;+
; Purpose: Interpolate selected variables to x-abscissa other than ut.
;-

pro stplot_newxx, var0, vars = vars, suffix = suffix

    if ~keyword_set(suffix) then suffix = ''

    ; read the new x-abscissa, interpolate to uniform.
    ; must be monotonic.
    tvar = tnames(var0)
    get_data, tvar, t0, x0
    nrec = n_elements(x0)
    xx = smkarthm(x0[0], x0[nrec-1], nrec, 'n')
    t1 = sinterpol(t0, x0, xx)
    
    ; interpolate selected variables to the new x-abscissa.
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], t0, x0, v0, limits = lim
        x1 = sinterpol(t0, x0, t1)
        store_data, vars[i]+suf, t1, x1, limits = lim
    endfor

end