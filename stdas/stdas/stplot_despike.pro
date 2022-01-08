
pro stplot_despike, vname, sigma, newname = newname

    if n_elements(sigma) eq 0 then sigma = 3

    get_data, vname, t0, f0, limits = lm
    dims = size(f0,/dimensions)
    if n_elements(dims) eq 1 then ndim = 1 $
    else ndim = dims[1]
    f1 = f0 & nrec = dims[0]
    for i = 0, ndim-1 do begin
        tf = f0[*,i]
        del = [0,[tf[1:nrec-1]-tf[0:nrec-2]]]
        avgval = mean(del)
        stdval = stddev(del)
        idx = where(abs(del) gt sigma*stdval, cnt)
        if cnt eq 0 then continue
        tf[idx] = !values.d_nan
        plot, tf, yrange = [-50,200]
        plot, f0[*,i], yrange = [-50,200]
        store_data, 'f1', t0, tf
        store_data, 'f2', t0, f0
        stop
        tplot, ['f1','f2']
    endfor
        

end