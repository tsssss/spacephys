;+
; calc the envolope of a time series, upper and lower
;-

function scalcenv, f0, width = width

    nrec = n_elements(f0)
    if nrec eq 0 then return, -1
    txs = findgen(nrec)

    if n_elements(width) eq 0 then width = nrec*0.01>5

;    f0s = smooth(f0, width*2, /nan, /edge_mirror)
    f0s = abs(f0)
    
    f1s = dblarr(nrec)
    flags = dblarr(nrec)
    
    for i = width, nrec-width-1 do begin
        f1s[i] = stddev(f0s[i-width:i+width], /nan)
        flags[i] = 1
        i = i+width
    endfor
    
    idx = where(flags eq 1, cnt)

    f1s = (cnt ge 4)? interpol(f1s[idx], idx, findgen(nrec),/spline): f0
    f1s = f1s>0


    return, f1s

end
