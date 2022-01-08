;+
; Type: function.
; Purpose: Moving median for given data.
; Parameters:
;   f0, in, dblarr[n], req. The given data.
; Keywords:
;   width, in, int, opt. Set the range for median calculation. Default is 3.
;   edge_wrap, in, boolean, opt. Wrap edge for padding, i.e., periodical.
;   edge_mirror, in, boolean, opt. Mirror edge for padding.
;   edge_truncate, in, boolean, opt. Pad with 0.
; Return: dblarr[n]. The moving median.
; Notes: Large width leads to smooth moving median. Do not set it large.
; History:
;   2015-02-15, Sheng Tian, create.
;-
function smedfilt, f0, width = width, $
    edge_wrap = wrap, edge_mirror = mirror, edge_truncate = trunc

    nrec = n_elements(f0)
    if n_elements(width) eq 0 then width = 3d
    f1 = f0
    w2 = fix(width/2)   ; force odd half-width to prevent even median.
    
    ; pad 0 on both sides.
    if keyword_set(wrap) then $
        f1 = [f0[nrec-w2:nrec-1],f0,f0[0:w2-1]] $
    else if keyword_set(mirror) then $
        f1 = [reverse(f0[0:w2-1]),f0,reverse(f0[nrec-w2:nrec-1])] $
    else if keyword_set(trunc) then $
        f1 = [dblarr(w2),f0,dblarr(w2)]
    nrec1 = n_elements(f1)-w2*2
    
    ; modify value on-the-go.
    for i = w2, nrec1-1+w2 do f1[i] = median(f1[i-w2:i+w2])
    
    if nrec ne nrec1 then return, f1 else return, f1[w2:w2+nrec-1]
end

nrec = 1000 & f0 = sin(findgen(nrec)/10)
idx0 = [10,40,300,310,700] & f0[idx0] = 10

f1 = smedfilt(f0)
f1 = smedfilt(f0, /edge_wrap)
f1 = smedfilt(f0, /edge_mirror)
end