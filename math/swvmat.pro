;+
; Type: function.
; Purpose: Perform moving average transform, or inverse transform.
;       (1) Do transform, return dblarr[n,m], when
;           f0 is dblarr[n], inverse is false;
;       (2) Do inverse transform, return dblarr[n], when
;           f0 is dblarr[n,m], inverse is true (optional);
;       (3) Do filter, return dblarr[n], when
;           f0 is dblarr[n], inverse is true.
;   Default scale is 50 log scales in [4,0.5*nrecs].
;   Default filter is [scale[0], scale[nscale-1]].
;   scale and filter are round to integer.
; Parameters:
;   f0, in, dblarr[n]/dblarr[n,m], required.
;       Input signal in dblarr[n] or 
;       transform func in dblarr[n,m].
;   order, in, int, optional. Default is 1.
; Keywords:
;   background, out, dblarr[n], optional. The left background field.
;   scale, in/out, dblarr[m], optional.
;       Set scales for transform or return the default scale.
;       Scale refers to number of record.
;   filter, in/out, dblarr[m], optional.
;       Set filter for inverse tranaform.
;       Filter refers to number of record.
;   inverse, in, boolean, optional.
;       Set inverse to do inverse transform or filter.
;   correct, in/out, dblarr[m], optional.
;       Set correct to do amplitude correction.
;       cor[i] is the Fourier period of channel wvt[i,*].
;       After correction, cor is in dblarr[m,n], the corrected wvt.
; Return:
;   return, out, dblarr[n]/dblarr[n,m].
; Notes: none.
; Dependence: none.
; Author: Sheng Tian.
; History: 2013-03-24, Sheng Tian, create.
;-

function swvmat, f0, order, background = bg, $
    scale = scl, filter = ftr, inverse = inverse, $
    correct = cor, boundary = boundary

    compile_opt idl2

    ; check f0.
    sz = size(f0)
    ndim = sz[0]
    nrec = sz[1]

    ; check order.
    if n_elements(order) eq 0 then order = 2

    ; check scale.
    if n_elements(scl) eq 0 then begin
        nscl = 50
        q = (0.125D*nrec)^(1D/(nscl-1))
        scl = dblarr(nscl)
        scl[0] = 4
        for i = 0, nscl-2 do scl[i+1] = scl[i]*q
    endif
    scl = round(scl)        ; do not smooth on fractional records.
    scl = scl[sort(scl)]    ; ensure from small to large.
    scl = double(scl[uniq(scl)])
    nscl = n_elements(scl)

    if ndim eq 1 then begin
        ; do transformation.
        dflt = f0
        wvt = dblarr(nrec, nscl)
        for i = 0, nscl-1 do begin
            ; do order 0.
            dnow = smooth(dflt, scl[i], /nan, /edge_mirror)
            wvt[*,i] = dflt-dnow    ; peel off the wanted layer.
            dflt = dnow             ; renew the onion.
            ; purify, higher orders.
            for j = 1, order-1 do begin
                dnow = smooth(wvt[*,i], scl[i], /nan, /edge_mirror)
                wvt[*,i] -= dnow    ; peel off the wanted layer.
                dflt += dnow        ; put the other stuff back.
            endfor
        endfor
        bg = dflt                   ; left background.
        wvt[*,nscl-1] += bg         ; put background back.

        ; add mask of boundary effect.
        if keyword_set(boundary) then begin
            fillval = !values.d_nan
            for i = 0, nscl-1 do begin
                if 4*scl[i] ge nrec then wvt[*,i] = fillval $
                else begin
                    wvt[0:2*scl[i],i] = fillval
                    wvt[-2*scl[i]:-1,i] = fillval
                endelse
            endfor
        endif
        if ~keyword_set(inverse) then return, wvt
    endif

    ; check ftr.
    if n_elements(ftr) eq 0 then ftr = scl[[0,nscl-1]]
    ftr = round(ftr)

    if ndim eq 2 then wvt = f0

    ; do correction.
    if n_elements(cor) ne 0 then begin
        ws = cor
        ns = n_elements(ws)
        u = !dpi*(ws##(1D/ws))
        a = sin(u)/u
        cor = wvt
        ; n-1 rounds of correction.
        ;       for i = ns-1, 1, -1 do begin
        ;           c = dblarr(i+1)
        ;           c[0] = 1
        ;           for j = 1, i do c[j] = c[j-1]*a[i,j-1]
        ;           cor[*,i-1] *= (1D/c[i-1])
        ;           for j = i-2, 0, -1 do cor[*,j] -= cor[*,i-1]*(c[j]-c[j+1])
        ;       endfor
        for i = ns-1, 1, -1 do begin
            ; coeffients related to the ith layer.
            c = dblarr(i+1)
            c[0] = 1
            for j = 1, i do c[j] = c[j-1]*a[i,j-1]
            cor[*,i] *= (1D/c[i])
            for j = i-1, 0, -1 do cor[*,j] -= cor[*,i]*(c[j]-c[j+1])
        endfor
        return, cor
    endif

    ; inverse transform or filter.
    idx = where(scl ge ftr[0] and scl le ftr[1])
    if idx[0] eq -1 then message, 'wrong filter ...'
    if n_elements(idx) eq 1 then return, reform(wvt[*,idx])
    return, total(wvt[*,idx], 2)

end

f0 = sin(findgen(1001)/1000*8*!dpi)
mat = swvmat(f0)
tvscl, mat
end

