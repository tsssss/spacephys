;+
; Type: procedure.
; Purpose: Split multi-dimensional data into individual tplot var.
; Parameter: var0, in, string, required. Original m-dim tplot var.
; Keyword: ytitles, in, strarr[m], optional. New ytitles.
;   labels, in, strarr[m], optional. New labels if var0 doesn't have labels.
;   newnames, in, strarr[m], optional. New var names.
;   ytitle, in, string, optional. New ytitle.
; Notes: New vars are default named by suffixing '_comp1', '_comp2', etc.
; Dependence: tdas,slib.
; Author: Sheng Tian.
; History: 2013-11-21, Sheng Tian, create.
;-
function stplot_split, var0, newnames = vnames, $
    ytitles = ytitles, labels = labels, colors = colors

    vname = tnames(var0[0])
    get_data, vname, t0, f0, limits = lm
    if size(f0,/n_dimensions) lt 2 then begin
        if n_elements(vnames) eq 0 then vnames = var0
        if n_elements(ytitles) ne 0 then options, vnames, 'ytitles', ytitles
        if n_elements(labels) ne 0 then options, vnames, 'labels', labels
        if n_elements(colors) ne 0 then options, vnames, 'colors', colors
        return, []
    endif else ndim = (size(f0,/dimensions))[1]
    idstrs = string(indgen(ndim)+1,format='(I0)')
    
    if n_elements(vnames) ne ndim then vnames = vname+'_comp'+idstrs
    if n_elements(ytitles) ne ndim then $
        if stagexist('ytitle',lm) then ytitles = replicate(lm.ytitle, ndim) $
        else ytitles = strarr(ndim)
    if n_elements(labels) ne ndim then $
        if stagexist('labels',lm) then labels = lm.labels else labels = idstrs
    if n_elements(colors) ne ndim then colors = intarr(ndim)
    
    for i = 0, ndim-1 do $
        store_data, vnames[i], t0, f0[*,i], $
            limits = {ytitle:ytitles[i], labels:labels[i], color:colors[i]}

    return, vnames
    
end
