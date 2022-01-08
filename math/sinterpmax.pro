
; target to interpolate latitude data from geopack.

function sinterpmax, f0;, maxdf = maxdf

    nrec = n_elements(f0)
    df = f0[1:nrec-1]-f0[0:nrec-2]
    
    tmp = df[sort(abs(df))]
    stdv = stddev(tmp[0:0.8*nrec])  ; exclude the largest points.
    stdvs = smvstddev(f0)
    idx = where(abs(df) gt stdvs*2, cnt)
    
    if cnt le 3 then return, f0     ; too few points to interpolate.
    
    for i = 0, cnt-2 do begin
        ti = idx[i]
        if ti eq -1 then continue
;        if ti-idx[i+1] eq -1 then if abs(df[idx[i+1]]) gt stdv then begin
;            idx[i] = -1
;            idx[i+1] = -1
;        endif
        if abs(df[ti-1]) gt stdvs[ti-1] then idx[i] = -1
        if abs(df[ti+1]) gt stdvs[ti+1] then idx[i] = -1
    endfor
    
    tmp = where(idx ne -1, cnt)
    if cnt le 3 then return, f0
    idx = idx[tmp]+1
    
    if cnt gt 3 then begin
        f1 = sinterpol(f0[idx], idx, indgen(nrec), /spline)
    endif else begin
        f1 = f0
    endelse

    return, f1

end
