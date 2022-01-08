;+
; difference from IDL's strsplit:
; 0, no match return idx=-1.
; 1, string and pattern are scalars.
; 2, do not miss '' between pattern.
;-
function sstrsplit, str0, ptn0, count = nstr, length = lens, extract = extract

    compile_opt idl2
    on_error, 2
    
    if n_elements(str0) eq 0 then message, 'no input ...'
    idxs = 0
    flag = 0    ; 1:match starts at 0, 2:no match.
    
    if n_elements(ptn0) eq 0 then message, 'no pattern ...'
    ptnlen = strlen(ptn0)
    if ptnlen ne 0 then begin
        pos0 = 0        ; where strpos starts to search.
        repeat begin
            idx = strpos(str0,ptn0,pos0)
            if idx eq -1 then break
            idx += ptnlen
            idxs = [idxs,idx]
            pos0 = idx
        endrep until idx eq -1
    endif
    
    nstr = n_elements(idxs)
    if nstr eq 1 then begin                ; no match.
        idxs = -1
        flag = 2
    endif else if idxs[1] eq 0 then begin   ; match starts at 0.
;        idxs = idxs[1:*]
;        nstr-= 1
        flag = 1
    endif
    if ~arg_present(lens) && ~keyword_set(extract) then return, idxs
    
    ; include the ending.
    lens = intarr(nstr)
    for i = 0, nstr-2 do lens[i] = idxs[i+1]-idxs[i]-ptnlen
    lens[i] = strlen(str0)-idxs[i]>0
    if ~keyword_set(extract) then return, idxs
    
    ; extract
    strs = strarr(nstr)
    for i = 0, nstr-1 do strs[i] = strmid(str0,idxs[i],lens[i])
    return, strs

end

str0 = 'a!Cb!C!C'
ptn0 = '!C'
ptn1 = 'a!C'
idx = sstrsplit(str0,ptn0,count = cnt, length = lens, /extract)
print, cnt, idx, lens
idx = sstrsplit(str0,ptn1,count = cnt, length = lens, /extract)
print, cnt, idx, lens
end