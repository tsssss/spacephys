;+
; Type: procedure.
; Purpose: Filter using Morlet wavelet transform.
; Parameters: vname, in, type = string, required. Signal varname in tplot.
; Keywords:
;   filter = filter, in, dblarr[2], req. Filter range in sec.
;   scaleinfo, in, struct, opt. Default is
;       {s0:4d*dr0, s1:0.5d*dur, dj:1d/8, ns:0d}.
;   tscale, in/out, dblarr[m], opt. Scales in time. Info used are
;       minscale, maxscale, nscale.
;   newname, in, string, optional. Output varname in tplot.
;   overwrite, in, boolean, optional. Set to overwrite vname.
;   ytitle, in, stirng, optional. Title for var.
; Notes: The variable should be 1-d data. This code forces label to be
;   specified by label keyword if vname contains no label.
; Dependence: tplot,slib.
; Author: Sheng Tian.
; History: 2017-07-25, Sheng Tian, create.
;-
pro stplot_filter_mor, vname, $
    filter = filter, $
    newname = newname, overwrite = overwrite, $
    tscale = tscales, scaleinfo = scaleinfo, $
    ytitle = ytitle, zrange = zrange


;---prepare vars.
    if n_elements(filter) ne 2 then message, 'wrong filter ...'
    pmin = min(filter, max=pmax)

    get_data, vname, uts, fs, limits = lim
    if n_elements(fs) eq 0 then message, 'no data ...'
    sz = size(fs,/dimensions)
    case n_elements(sz) of
        1: ndim = 1
        2: ndim = (size(fs,/dimensions))[1]
        else: message, 'wrong dimension ...'
    endcase
    dr0 = sdatarate(uts) & dr1 = 1d/dr0
    nrec = n_elements(uts)
    dur = dr0*nrec
    
    ; remove nan.
    idx = where(finite(fs,/nan), cnt)
    if cnt ne 0 then fs[idx] = 0

    
;----prepare scales.
    ; input: tscales, scaleinfo.
    ; output: s0,s1,dj,ns,scaleinfo.
    if n_elements(scaleinfo) eq 0 then $
        scaleinfo = {s0:4d*dr0, s1:0.5d*dur, dj:1d/8, ns:0d}

    if n_elements(tscales) eq 0 then begin          ; no scales.
        s0 = scaleinfo.s0
        s1 = scaleinfo.s1
        dj = scaleinfo.dj
        ns = scaleinfo.ns

        if s0 eq 0 then s0 = 4*dr0
        if s1 eq 0 then s1 = 0.5d*dur
        if dj eq 0 and ns eq 0 then dj = 1d/8
        if ns eq 0 then begin
            j1 = floor(alog(s1/s0)/alog(2)/dj)  ; # of powers-of-two with dj
            s1 = s0*2d^(dj*j1)            
            ns = j1+1
        endif
        if dj eq 0 then dj = alog(s1/s0)/alog(2)/(ns-1)
    endif else begin
        s0 = min(tscales)
        s1 = max(tscales)
        ns = n_elements(tscales)
        dj = alog(s1/s0)/alog(2)/(ns-1)
    endelse
    j1 = ns-1
    w0 = 6d
    cdelta = 0.776d     ; constant for w0=6, for normalization.
    psi0 = !dpi^(-0.25)

    scaleinfo.s0 = s0   ; min scale.
    scaleinfo.s1 = s1   ; max scale.
    scaleinfo.dj = dj   ; 2^dj scale spacing.
    scaleinfo.ns = ns   ; # of scales.



    if keyword_set(overwrite) then newname = vname
    if ~keyword_set(newname) then newname = vname+'_filter'

;----wavelet analysis.
    ; wavelet transform, in X.
    ffs = dblarr(nrec,ndim)
    for i=0, ndim-1 do begin
        f0 = fs[*,i]
        mor = wavelet(f0, dr0, /pad, s0=s0, dj=dj, j=j1, $
            mother='Morlet', param=w0, $
            period = ps, scale=ss, coi=coi)
        idx = where(ps ge pmin and ps le pmax, cnt)
        if cnt eq 0 then message, 'useless filter ...'
        ; re-construct waveform according to filter.
        mor = real_part(mor[*,idx])
        ss = ss[idx]
        for j=0, cnt-1 do ffs[*,i]+= mor[*,j]/sqrt(ss[j])
        ffs[*,i]*= (dj*sqrt(dr0))/cdelta/psi0
    endfor

    
    store_data, newname, uts, ffs, limits = lim
end
