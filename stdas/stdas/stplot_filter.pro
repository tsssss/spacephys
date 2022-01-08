;+
; Type: procedure.
; Purpose: Filter in tplot for mat, gauss, morlet, das (detrend&smooth).
; Parameters:
;   vname, in, string, req. Spectrogram varname in tplot.
;   mode, in, string, opt. Default 'mat', can be 'gauss','morlet','das'.
; Keywords: newname, in, string, optional. Newname for filtered signal.
;   filter, in, dblarr[2], req. Filter in time.
;   labels, in, string, opt. Default is 'filter0-filter1 s'.
;   ytitle, in, string, opt. Default is ''.
;   colors, in, int, opt. Default is -1.
;   limits, in, struct, opt. Limits for store_data.
;   ifilter, in, string/integer, opt. Filter's index.
; Notes: The usage is strict, no smart behaviour.
; Dependence: tplot,slib.
; Author: Sheng Tian.
; History: 2013-11-24, Sheng Tian, create.
;-
pro stplot_filter, vname, mode, newname = newname, filter = tfilter, $
    labels = labels, ytitle = ytitle, colors = colors, $
    ifilter = ifilter, limits = lm
    
    if n_elements(mode) eq 0 then mode = 'mat'
    if n_elements(ifilter) eq 0 then ifilter = '1'
    if size(ifilter,/type) ne 7 then ifilter = string(ifilter,format='(I0)')
;    ifilter = 'f'+ifilter       ; 'f' for filter.
    if mode eq 'das' then ifilter = '_'+mode+ifilter         ; 'x'->'x_dasf1'.
    if n_elements(newname) eq 0 then newname = vname+ifilter ; 'x'->'x_matf1'.
    if n_elements(tfilter) ne 2 then message, 'filter is undefined ...'
    tfilter = tfilter[sort(tfilter)]
    f1 = min(tfilter) & f2 = max(tfilter)
    sfilter = snum2str(tfilter[0])+'-'+snum2str(tfilter[1])+'s'

    get_data, vname, t0, s0, tscales, limits = tmp
    if n_elements(labels) eq 0 then $
        labels = stagexist('labels',tmp)? tmp.labels+'!C'+sfilter: sfilter
    if n_elements(ytitle) eq 0 then $
        ytitle = stagexist('ytitle',tmp)? tmp.ytitle: ''
    if n_elements(colors) eq 0 then colors = -1
    
    ; get filtered signal f0.
    if mode eq 'mat' then begin
        idx1 = where(tscales ge tfilter[0], cnt)
        idx1 = idx1[0]
        if tfilter[0] eq tscales[idx1] then if idx1 gt 0 then idx1+= 1
        idx2 = where(tscales le tfilter[1], cnt)    ; include upper limit.
        idx2 = idx2[cnt-1]
        if tfilter[1] gt max(tscales) then idx2 = n_elements(tscales)-1
        if idx2-idx1 lt 0 then begin    ; no band between the two filters.
            message, 'bad filter: '+sfilter, /continue
            f0 = fltarr(n_elements(s0[*,0]))
        endif else if idx2-idx1 eq 0 then begin
            f0 = s0[*,idx1:idx2]
        endif else f0 = total(s0[*,idx1:idx2],2)
    endif else if mode eq 'das' then begin
        dr1 = 1D/sdatarate(t0)
        f0 = smooth(s0,f1*dr1,/nan,/edge_wrap)-smooth(s0,f2*dr1,/nan,/edge_wrap)
    endif else if mode eq 'gauss' then begin
        ; gauss wavelet.
    endif else if mode eq 'morlet' then begin
        ; morlet wavelet.
    endif
    
    if n_elements(lm) eq 0 then $
        lm = {labels:labels, ytitle:ytitle, colors:colors}
    store_data, newname, t0, f0, limits = lm
end
