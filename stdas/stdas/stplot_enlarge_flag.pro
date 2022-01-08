;+
; flag can be 0 or 1. enlarge the part equal to flag0 by drec or dt.
;-
pro stplot_enlarge_flag, tvar, dut = dut, drec = drec, flag0 = flag0

    if n_elements(drec0) eq 0 then begin
        if n_elements(dut) eq 0 then return
        get_data, tvar, uts
        dr = sdatarate(uts)
        drec0 = ceil(dut/dr)
    endif
    drec = ceil(drec0)
    
    if n_elements(flag0) eq 0 then flag0 = 1
    
    get_data, tvar, uts, flag
    nrec = n_elements(uts)
    idx = where(flag eq flag0, cnt)
    for i = 0, cnt-1 do flag[(idx[i]-drec)>0:(idx[i]+drec)<(nrec-1)] = flag0
    store_data, tvar, uts, flag
    
end