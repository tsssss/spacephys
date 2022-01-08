;+
; Type: function.
; Purpose: Find node(s) of a curve within an error,
;   the slope at node(s) has a given sign.
; Parameters:
;   dat, in, dblarr[n], req. Data curve.
;   err, in, double, req. Error of 0.
;   sign, in, boolean, required. Sign of slope. 1 for +, -1 for -, 0 for both.
; Keywords: none.
; Return: return, intarr[m]. Index for nodes, return -1 if no node.
; Notes: This is a rough search, no interpolation, use the resolution
;     of given data. Consecutive nodes will shrink to 1st one.
; Dependence: none.
; History:
;   2013-06-18, Sheng Tian, create.
;-

function sfindnode, dat, err, sign = sign
    compile_opt idl2 & on_error, 2

    if n_elements(err) eq 0 then $
        err = 0.5   ; good for 1-min data, larger for lower time res.
    if n_elements(sign) eq 0 then sign = 0    ; both signs of slope.
    
    ids = where(abs(dat) lt err, nid)
    if nid eq 0 then return, -1
    
    ; shrink consecutive index.
    idx = ids[0]
    for i = 1, nid-1 do begin
        if ids[i]-ids[i-1] eq 1 then continue
        idx = [idx, ids[i]]
    endfor
    
    ; check slope.
    if sign eq 0 then return, idx
    ids = idx
    for i = 0, n_elements(idx)-1 do begin
        del = (idx[i] eq 0)? dat[idx[i]+1]-dat[idx[i]]:dat[idx[i]]-dat[idx[i]-1]
        ids[i] = del*sign ge 0
    endfor
    ids = where(ids eq 1)
    if ids[0] eq -1 then return, -1
    return, idx[ids]

end
