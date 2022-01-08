;+
; Type: procedure.
; Purpose: Sum components to one tplot var.
; Parameter: vars, in, strarr[m], required. Individual components.
; Keyword: newname, in, string, required. Combined tplot var name.
;   ytitle, in, string, optional. Default is ''.
;   labels, in, string, optional. Default is ''.
;   delete, in, boolean, optional. Set to delete old vars.
; Notes: none.
; Dependence: tdas,slib.
; Author: Sheng Tian.
; History: 2013-11-21, Sheng Tian, create.
;-
pro stplot_total, vars, newname = newname, delete = delete, $
    ytitle = ytitle, labels = labels, colors = colors, limits = lm
    
    if n_elements(newname) eq 0 then message, 'no newname ...'

    ndim = n_elements(vars)
    if ndim eq 0 then message, 'no variables ...'
    get_data, vars[0], t0, s0, limits = lm
    
    ; ytitle and labels.
    if n_elements(ytitle) eq 0 then ytitle = ''
    if n_elements(labels) eq 0 then labels = ''
    if n_elements(colors) eq 0 then $
        if stagexist('colors',lm) then colors = lm.colors

    ; sum signal.
    for i = 1, ndim-1 do begin
        get_data, vars[i], t0, f0
        s0 += f0
    endfor

    if n_elements(lm) eq 0 then $
        lm = {ytitle:ytitle, labels:labels}
    store_data, newname, t0, s0, limits = lm
    if n_elements(colors) ne 0 then options, newname, 'colors', colors
        
    if keyword_set(delete) then store_data, vars, /delete
end