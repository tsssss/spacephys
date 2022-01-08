;+
; Type: procedure.
; Purpose: Remove spikes in given m-dimension line plots, using median.
; Parameters:
;   f0, in/out, dblarr[n,m]/dblarr[n]/string, req. Data, string for tplot var.
;   t0, in/out, dblarr[n], opt. Time or x data. If f0 is tplot var, then t0
;       will be loaded, otherwise t0 can be set or generated with indgen.
; Keywords:
;   width, in, int, opt. Set the range for median calculation. Default is 7.
;   nsigma, in, float, opt. Set threshold to restore original value.
; Notes: Large width leads to smooth moving median. Do not set it large.
; Dependence: none.
; History:
;   2014-10-01, Sheng Tian, create.
;-
pro sdespike_median, t0, f0, width = width, nsigma = nsigma, _extra = extra

    ; prepare t0 and f0.
    if n_params() eq 1 then begin
        if size(t0,/type) eq 7 then begin   ; tplot var.
            vname = t0 & get_data, vname, t0, f0
        endif else f0 = t0                  ; set f0 not t0.
    endif
    
    tmp = size(f0) & if tmp[0] gt 2 then message, 'more than 2 dimensions ...'
    ndim = (tmp[0] eq 1)? 1: tmp[2]
    nrec = n_elements(f0)/ndim
    if n_elements(vname) eq 0 and n_params() eq 1 then t0 = findgen(nrec)
    
    ; empirical value.
    if n_elements(width) eq 0 then width = 7        ; median among 7 points.
    if n_elements(nsigma) eq 0 then nsigma = 0.2
    
    ; work on each dimension.
    for i = 0, ndim-1 do begin
        tf0 = f0[*,i]
        tf1 = smedfilt(tf0, width = width, _extra = extra)
        ; restore good points' value.
        del = abs(tf1-tf0)
        sigma = stddev(tf0)
        idx = where(del le nsigma*sigma, cnt)
        if cnt gt 0 then tf1[idx] = tf0[idx]
        f0[*,i] = tf1
    endfor
    
    ; store data.
    if n_elements(vname) eq 1 then store_data, vname, t0, f0 $
    else if n_params() eq 1 then t0 = f0
end