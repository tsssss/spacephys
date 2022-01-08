;+
; vars0 can be valid input for tnames, must be 1-d data.
; yr is in/out.
; set tr to unify yr in the time range.
;-
pro stplot_uniyr, vars0, trange = tr, yrange = yr, ylog = ylog

    vars = tnames(vars0)
    nvar = n_elements(vars)
    if nvar eq 0 then return
    
    if n_elements(yr) ne 2 then begin
        ymins = dblarr(nvar) & ymaxs = dblarr(nvar)
        for i = 0, nvar-1 do begin
            get_data, vars[i], t0, dat
            if n_elements(tr) eq 2 then begin
                idx = where(t0 ge tr[0] and t0 le tr[1], cnt)
                if cnt gt 0 then tymin = min(dat[idx], max = tymax, /nan) $
                else tymin = min(dat, max = tymax, /nan)
            endif else tymin = min(dat, max = tymax, /nan)
            ymins[i] = tymin & ymaxs[i] = tymax
        endfor
        yr = [min(ymins),max(ymaxs)]
    endif
    
    options, vars, 'yrange', yr
    if keyword_set(ylog) then options, vars, 'ylog', ylog
end