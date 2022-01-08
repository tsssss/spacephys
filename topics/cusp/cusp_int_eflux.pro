
pro cusp_int_eflux, vname, ilatname, tr = tr
    re = 6378.137d & h0 = 100 ; km, altitude.
    c = (re+h0)*!dpi/180d
    
    get_data, vname, t0, eflux0
    get_data, ilatname, tmp, ilat
    ilat = sinterpol(ilat, tmp, t0)
    nrec = n_elements(t0)
    
    idx = keyword_set(tr)? where(t0 ge tr[0] and t0 le tr[1]): indgen(nrec)
    ilat = ilat[idx]
    ndim = size(eflux0,/n_dimensions)
    if ndim eq 1 then dims = 1 else dims = (size(eflux0,/dimensions))[1]
    eflux_int = dblarr(dims)
    for i = 0, dims-1 do begin
        eflux = eflux0[idx,i]
        nrec = n_elements(eflux)
        v1 = 0.5*(eflux[0:nrec-2]+eflux[1:nrec-1])
        v2 = abs(ilat[1:nrec-1]-ilat[0:nrec-2])*c
        eflux_int[i] = total(v1*v2,/nan)
    endfor
    
    print, vname, string(eflux_int,format='(f10.2)'), ' W/m'
    store_data, vname, t0, eflux0, eflux_int
end