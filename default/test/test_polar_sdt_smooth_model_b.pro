
pro test_polar_sdt_smooth_model_b
    
    ; Problem:
    ;   We want to get a better estimation of the Polar model field.
    ; The basic idea is to determine the width of "active period", in this case
    ; is the cusp. any variation longer than width are considered to be
    ; background field.

    fn = sdiskdir('Works')+'/data/cusp/po_sdt_fld_1998_0925_05.sdt'
    sdt = ssdtread(fn)

    ; get uniform time.
    maxnrec = 100000L   ; prevent memory overflow.
    t0 = sdt.var.polar_b_spc_z.depend_0
    dr = sdatarate(t0) & nrec = n_elements(t0)
    if nrec gt maxnrec then dr = (t0[nrec-1]-t0[0])/maxnrec
    t0 = smkarthm(t0[0], t0[nrec-1], dr, 'dx')
    nrec = n_elements(t0)
    print, 'data rate: ', dr
    ; b field.
    ft  = sdt.var.polar_b_spc_z.depend_0
    fxy = sdt.var.polar_b_spc_x_y.value
    f56 = sdt.var.polar_b_spc_56.value
    fz  = sdt.var.polar_b_spc_z.value
    b_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
    b_spc = sinterpol(b_spc, ft, t0)
    sdespike, t0, b_spc     ; remove spike.
    store_data, 'b_spc', t0, b_spc
    stplot_split, 'b_spc'

    ft  = sdt.var.polar_model_b_t96_spc_z.depend_0
    fxy = sdt.var.polar_model_b_t96_spc_x_y.value
    f56 = sdt.var.polar_model_b_t96_spc_56.value
    fz  = sdt.var.polar_model_b_t96_spc_z.value
    bmod_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
    bmod_spc = sinterpol(bmod_spc, ft, t0)

    nsec = 20
    width = nrec/nsec
    ratio = 0.4
    ratio1 = (1-ratio)*0.5
    weight = 0.5*(1-sin(smkarthm(0,!dpi,nrec,'n')))+1d  ; less picky on edge.
    
    ndim = 3 & f0 = b_spc
    del = 5    ; nT.
    for i = 0, ndim-1 do begin
        tf0 = f0[*,i]
        df0 = [0d,tf0[1:nrec-1]-tf0[0:nrec-2]]
        dfstddev = smvstddev(df0, t0, width = width, ratio = ratio1, /quadratic)
        tmp = median(dfstddev)
        idx = where(dfstddev gt weight*tmp)
        tf0[idx] = !values.d_nan
        window, 1
        plot, f0[*,i]                           ; original field.
        oplot, tf0+del, color = sgcolor('blue')  ; exclude active field.
        oplot, bmod_spc[*,i], color = sgcolor('green')
        
        ; moving mean.
        y0 = tf0
        nsec = nrec/width
        x1 = findgen(nsec+1)*(nrec-1)/nsec
        mvtmp = findgen(nsec)
        for j = 0, nsec-1 do begin
            tmp = y0[x1[j]:x1[j+1]]
            mvtmp[j] = mean((tmp[sort(tmp)])[ratio*width:(1-ratio)*width])
        endfor
        x1 = (x1[0:nsec-1]+x1[1:nsec])*0.5      ; use middle of section.
        idx = where(~finite(mvtmp,/nan))
        x1 = x1[idx]
        mvtmp = mvtmp[idx]
        mvtmp = interpol(mvtmp, x1, findgen(nrec), /spline)
        oplot, mvtmp, color = sgcolor('red')
        stop
    endfor
end