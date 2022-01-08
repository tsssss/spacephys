;+
; Type: function.
;
; Purpose: Calc background field of given field, which should have quiet edges.
;
; Parameters:
;   f0, in, dblarr[n] or dblarr[n,m], required. Given data.
;
; Keywords:
;   nsection, in, int, optional. Default is 20, larger for finer background.
;
; Return:
;   dblarr[n] or dblarr[n,m]. The background field.
;
; Notes: none.
;
; Dependence: slib.
;
; History:
;   2015-06-18, Sheng Tian, create.
;-
function scalcbg, f0, nsection = nsec, pure_smooth = pure_smooth

    dims = size(f0, /dimensions)
    ndim = size(f0, /n_dimensions)
    if ndim eq 1 then begin
        nrec = dims[0] & ndim = 1
    endif else if ndim eq 2 then begin
        nrec = dims[0] & ndim = dims[1]
    endif else message, 'wrong dimension ...'
    
    ; background.
    fbg = f0
    
    if n_elements(nsec) eq 0 then nsec = 20 ; 20 sections.
    width = nrec/nsec
    ratio = 0.4
    ratio1 = (1-ratio)*0.5
    weight = 1-sin(smkarthm(0,!dpi,nrec,'n'))   ; less picky on edge.
;    weight = 0.5*cos(smkarthm(0,2*!dpi,nrec,'n'))+0.5
    
    for i = 0, ndim-1 do begin
        tf0 = f0[*,i]
        if ~keyword_set(pure_smooth) then begin
            df0 = [0d,tf0[1:nrec-1]-tf0[0:nrec-2]]
            dfstddev = smvstddev(df0,t0, width=width,ratio=ratio1, /quadratic)
            tmp = median(dfstddev)
            amp = max([dfstddev[0:width],dfstddev[nrec-width:nrec-1]])/tmp
            idx = where(dfstddev gt (amp*weight+1d)*tmp)    ; ensure pass on edges.
            tf0[idx] = !values.d_nan    ; throw large fluctuations.
        endif
        ; calc a 'moving mean'.
        x1 = findgen(nsec+1)*(nrec-1)/nsec
        y1 = findgen(nsec)
        for j = 0, nsec-1 do begin
            tmp = tf0[x1[j]:x1[j+1]] ; use points around median.
            y1[j] = mean((tmp[sort(tmp)])[ratio*width:(1-ratio)*width])
        endfor
        x1 = (x1[0:nsec-1]+x1[1:nsec])*0.5      ; use middle of section.
        if ~keyword_set(pure_smooth) then idx = where(~finite(y1,/nan), cnt) else idx = findgen(nsec)
        cnt = n_elements(idx)
        case cnt of
            0: tf0 = fltarr(nrec)
            1: tf0 = fltarr(nrec)+y1[idx]
            2: tf0 = fltarr(nrec)+mean(y1[idx])
            3: tf0 = interpol(y1[idx],x1[idx],findgen(nrec), /quadratic)
            else: tf0 = interpol(y1[idx],x1[idx],findgen(nrec), /quadratic)
        endcase
        fbg[*,i] = tf0
    endfor

    return, fbg
end
