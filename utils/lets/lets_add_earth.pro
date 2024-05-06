
function lets_add_earth, xrange=xr, yrange=yr, npoint=npoint

    ; Add earth.
    if n_elements(npoint) eq 0 then npoint = 40
    tmp = smkarthm(0,2*!dpi, npoint, 'n')
    txs = cos(tmp)
    tys = sin(tmp)

    if n_elements(xr) ne 2 then xr = [-1,1]
    index = where_pro(txs, '[]', minmax(xr), count=count)
    if count eq 0 then return, -1
    txs = min(xr)>txs<max(xr)

    if n_elements(yr) ne 2 then yr = [-1,1]
    index = where_pro(tys, '[]', minmax(yr), count=count)
    if count eq 0 then return, -1
    tys = min(yr)>tys<max(yr)

    index = where(txs ge 0, count)
    if count ne 0 then polyfill, txs>0, tys, color=sgcolor('white')
    index = where(txs le 0, count)
    if count ne 0 then polyfill, txs<0, tys, color=sgcolor('grey')

    plots, txs, tys

end