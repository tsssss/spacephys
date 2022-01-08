;+
; Type: procedure.
; Purpose: Save certain dimensions of given var to a new var.
; Parameter: var0, in, string, required. Original m-dim tplot var.
;   index, in, intarr[m], required. Index for the wanted dimensions. 0-based.
; Keyword: ytitles, in, strarr[m], optional. New ytitles.
;   labels, in, strarr[m], optional. New labels if var0 doesn't have labels.
;   newnames, in, strarr[m], optional. New var names.
;   ytitles, in, string, optional. New ytitle.
;   colors, in, intarr[m], optional. New color.
;   nosplit, in, boolean, optional. Set to save wanted dimensions in one var.
; Notes: New vars are default named by suffixing '_comp1', '_comp2', etc.
; Dependence: tdas,slib.
; Author: Sheng Tian.
; History: 2013-11-21, Sheng Tian, create.
;-
pro stplot_index, var0, idx, newnames = vnames, ytitles = ytitles, $
    labels = labels, colors = colors, nosplit = nosplit

    vname = tnames(var0[0])
    get_data, vname, t0, f0, limits = lm
    if size(f0,/n_dimensions) ne 2 then return
    if n_elements(idx) eq 0 then idx = indgen((size(f0,/dimensions))[1])
    ndim = n_elements(idx)
    idstrs = string(idx+1,format='(I0)')
    if n_elements(colors) eq 0 then colors = intarr(ndim)
    
    if n_elements(vnames) ne ndim then vnames = vname+'_comp'+idstrs
    if n_elements(ytitles) ne ndim then $
        if stagexist('ytitle',lm) then ytitles = replicate(lm.ytitle, ndim) $
        else ytitles = strarr(ndim)
    if n_elements(labels) ne ndim then $
        if stagexist('labels',lm) then labels = lm.labels[idx] $
        else labels = idstrs
    
    if keyword_set(nosplit) then $
        store_data, vnames[0], t0, f0[*,idx], $
            limits = {ytitle:ytitles[0], labels:labels[0], colors:colors} $
    else for i = 0, ndim-1 do $
        store_data, vnames[i], t0, f0[*,idx[i]], $
            limits = {ytitle:ytitles[i], labels:labels[i], colors:colors[i]}
end