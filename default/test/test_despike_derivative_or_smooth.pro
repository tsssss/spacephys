;+
; Conclusion: use derivative to despike.
; There are 2 ways to despike: (1) large derivative. (2) large deviation from
; moving average.
; (2) is more complex but gives essentially the same result from (1). (2) also
; suffer edge problem.
;-
pro test_despike_derivative_or_smooth

    ; Problem:
    ;     We want to despike B field for Polar but preserve large 
    ; amplitude wave.
    
    fn = sdiskdir('SHENG')+'/works/polarcap/data/po_sdt_fld_1998092415.sdt'
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
    store_data, 'b_spc', t0, b_spc
    stplot_split, 'b_spc'
    
    ndim = 3 & f0 = b_spc
    for i = 0, ndim-1 do begin
        istr = string(i+1,format='(I0)')
        tf = f0[*,i]
        
        ; derivative.
        del1 = [0,[tf[1:nrec-1]-tf[0:nrec-2]]]
        avgval = mean(del1)
        stdval = stddev(del1)
        idx = where(abs(del1) gt 3*stdval, cnt)
        help, idx
        vname = 'b_spc_comp'+istr+'_del'
        store_data, vname, t0, del1
        ylim, vname, -3*stdval, 3*stdval

        ; smooth.
        smtval = smooth(tf,4,/edge_wrap)
        del2 = tf-smtval
        avgval = mean(del2)
        stdval = stddev(del2)
        idx = where(abs(del2) gt 3*stdval, cnt)
        help, idx
        vname = 'b_spc_comp'+istr+'_smth'
        store_data, vname, t0, del2
        ylim, vname, -3*stdval, 3*stdval
    endfor
    
    vars = 'b_spc_comp'+['1*','2*','3*']
    tplot, vars
    
end