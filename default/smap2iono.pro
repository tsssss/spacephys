; h0, altitude in km, default is 100.

pro smap2iono, v1, v2, b = b, h0 = h0, newname = newname, coef = coef

    if n_elements(h0) eq 0 then h0 = 100
    if ~keyword_set(newname) then newname = v1
    re = 6378.137d  ; km.
    r0 = h0/re+1d

    ; read data.
    get_data, v1, t0, dat1
    get_data, v2, tmp, dat2
    dat2 = sinterpol(dat2, tmp, t0)
    ndim = size(dat1,/n_dimensions)
    
    if keyword_set(b) then begin    ; dat2 is mlat in degree.
        get_data, b, tmp, bmag
        if size(bmag, /n_dimensions) eq 2 then bmag = sqrt(total(bmag^2,2))
        bmag = sinterpol(bmag, tmp, t0)
        tmp = sqrt(3*sin(dat2*(!dpi/180))^2+1d)
        me = 31025.2    ; nT.
        coef = me/(r0^3)*tmp/bmag
        if ndim eq 1 then dat1*= coef else $
            for j = 0, (size(dat1,/dimensions))[1]-1 do dat1[*,j]*= coef
    endif else begin        ; dat2 is distance.
        coef = (dat2*(1d/r0))^3
        if ndim eq 1 then dat1*= coef else $
            for j = 0, (size(dat1,/dimensions))[1]-1 do dat1[*,j]*= coef
    endelse
    
    store_data, newname, t0, dat1
end