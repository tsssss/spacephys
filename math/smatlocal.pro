;+
; Type: procedure.
; Purpose: Get local freq info using moving average transform.
; Parameters: f0, in, dblarr[n], req. Input signal.
;   t0, in, int, req. Location in record.
; Keywords: none.
; Notes: none.
; Dependence: none.
; Author: Sheng Tian.
; History: 2014-01-08, Sheng Tian, create.
;-

pro smatlocal, f0, t0, scale = scl, edge = edge, $
    error = error, ratio = ratio

    compile_opt idl2

    ff = f0 & tt = t0
    nrec = n_elements(ff)
    if n_elements(error) eq 0 then error = max(ff,/absolute)*0.02
    if n_elements(ratio) eq 0 then ratio = 0.98

    ; edge options.
    if n_elements(edge) eq 0 then edge = 'period'
    case edge of
        'mirror': ef = reverse(ff)
        'period': ef = ff
    endcase

    ; scale range.
    if n_elements(scl) eq 0 then scl = [2,floor(0.5*nrec)]

    ; center at tt.
    if tt+scl[1] gt nrec then begin
        ff = [ff,ef]
    endif else if tt-scl[1] le 0 then begin
        ff = [ef,ff] & tt += nrec
    endif

    for i = scl[0], scl[1]-1 do begin
        ff0 = ff
        ff1 = ff0-smooth(ff0, 2*i+1, /nan, /edge_wrap)
        maxff1 = max(ff1[tt-i:tt+i], /absolute)
        if maxff1 le error then continue
        ff2 = ff1-smooth(ff1, 2*i+1, /nan, /edge_wrap)
        if max(ff2[tt-i:tt+i], /absolute) lt ratio*maxff1 then continue
        ff = smooth(ff, 2*i+1, /nan, /edge_wrap)
        print, i
        plot, ff2
    endfor
end

x0 = findgen(1001)*(3*!dpi/1000)
f0 = sin(x0)
smatlocal, f0, x0
end
