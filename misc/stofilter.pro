;+
; Type: function.
; Name: stofilter.
; Purpose: convert [f0,f1,...,fn] to [[f0,f1],[f1,f2],...].
; Parameters: ifilters, in, type = dblarr[n] or dblarr[n,2]. In filters.
; Keywords: none.
; Return: ofilters, out, type = dblarr[n,2]. Out filters.
; Example: none.
; Notes: if in filter is in [n,2], then do nothing.
; Dependence: none.
; Author: Sheng Tian.
; History: 2013-11-11, Sheng Tian, create.
;-
function stofilter, ifilters

    nifilter = n_elements(ifilters)
    
    ; already in [[f0,f1],[f1,f2],...].
    if nifilter ge 4 and size(ifilters,/n_dimensions) eq 2 then return, ifilters
    if nifilter le 1 then message, 'wrong filter ...'
    
    nofilter = nifilter-1
    ofilters = dblarr(nofilter,2)
    for i = 0, nofilter-1 do ofilters[i,*] = [ifilters[i],ifilters[i+1]]
    return, ofilters

end