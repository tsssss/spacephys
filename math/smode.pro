;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

function smode, v0, error = err

    if n_elements(err) eq 0 then err = max(v0,/absolute)*0.01
    er1 = 1d/err
    
    v1 = double(v0)
    v1 = round(er1*v1)*err
    
    v2 = v1[uniq(v1,sort(v1))]
    nf = n_elements(v2)
    fs = intarr(nf)
    
    for i = 0, nf-1 do begin
        tmp = where(v1 eq v2[i], cnt)
        fs[i] = cnt
    endfor

    tmp = max(fs,idx)
    return, v2[idx]
    
end

v0 = [0.1,0.1,0.2,0.5]
print, smode(v0)
end
