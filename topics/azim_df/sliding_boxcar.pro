;+
; Calculates the moving boxcar median, mean, stddev, etc.
; Pad 0 around both sides.
;
; xxs. An array of input data.
; width. A number for the width of box in # of record.
; type=. A string sets the mode, can be 'median','mean','stddev'.
; ratio=. A number in [0,1]. Set to sort data in abs value, use [0,ratio].
;-

function sliding_boxcar, xxs, width, type=type, errmsg=errmsg, ratio=ratio

    nxx = n_elements(xxs)
    if nxx eq 0 then return, !null
    if n_elements(width) eq 0 then begin
        errmsg = handle_error('No input width ...')
        return, xxs
    endif
    ww = round(width)
    if ww lt 1 then return, xxs
    if n_elements(type) eq 0 then type = 'mean'
    if n_elements(ratio) ne 0 then begin
        ratio = ratio>0.
        ratio = ratio<1.
    endif

    x1s = xxs[*]
    index = where(finite(x1s,/nan), count)
    if count ne 0 then x1s[index] = 0

    w1 = ((ww mod 2) eq 1)? (ww-1)/2: ww/2
    w2 = ((ww mod 2) eq 1)? ww-w1-1: ww-w1
    x0s = [fltarr(w1),x1s,fltarr(w2)]
    for ii=w1,nxx-1+w2 do begin
        txx = x0s[ii-w1:ii+w2]
        if n_elements(ratio) ne 0 then begin
            txx = txx[sort(abs(txx))]
            txx = txx[0:ratio*ww]
        endif
        x1s[ii-w1] = call_function(type, txx)
    endfor

    return, x1s
end

xxs = [1,-2,3,5,-8,2,1]
type = 'mean'
yys = boxcar(xxs, 2, type=type)
yys = boxcar(xxs, 3, type=type)
end
