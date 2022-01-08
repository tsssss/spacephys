;+
; convert data into step style.
;-
pro bias_param_to_step, tvar, newname = newname, delete = delete, round = round
    
    get_data, tvar, t0, dat
    nrec = n_elements(t0)
    
    if keyword_set(round) then begin
        v0s = dat   ; original data.
        v1s = v0s   ; rounded data.
        for i = 0, n_elements(v1s)-1 do v1s[i] = round(v1s[i])
        v2s = v1s[uniq(v1s,sort(v1s))]  ; the rough values.
        v3s = v2s[1:*]-v2s[0:-2]        ; diff in rough values.
        err = sqrt(abs(min(v3s)*max(v3s)))
        idx = where(v3s gt err)
        v4s = v2s[idx+1]                ; the final values.
        if idx[0] eq 0 then v4s = [v2s[0],v4s]
        for i = 0, n_elements(v4s)-1 do begin
            idx = where(v0s ge v4s[i]-err and v0s le v4s[i]+err, cnt)
            v1s[idx] = mean(v0s[idx])
        endfor
        if n_elements(newname) eq 0 then newname = tvar
        store_data, newname, t0, v1s
        return
    endif
    
    t1 = t0[0]
    tmp = dat[0]
    for i = 1, nrec-1 do begin
        t1 = [t1,t0[i],t0[i]]
        tmp = [tmp,dat[i-1],dat[i]]
    endfor
    if n_elements(newname) eq 0 then newname = tvar
    store_data, newname, t1, tmp
    options, newname, 'psym', -4
    
    if keyword_set(delete) then store_data, tvar, /delete

end
