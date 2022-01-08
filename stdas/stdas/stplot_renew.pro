;+
; Type: procedure.
; Purpose: Renew name, ytitle, labels, colors.
; Parameter: vname, in, string, required. Old name.
; Keyword: newname, in, string, required. New name.
;   ytitle, in, string, optional. Default is ''.
;   labels, in, strarr[m], opt. Default is ''.
;   colors, in, intarr[m], opt. Color for component.
;   delete, in, boolean, opt. Set to delete the old var.
; Notes: none.
; Dependence: tdas,slib.
; Author: Sheng Tian.
; History: 2013-11-25, Sheng Tian, create.
;-
pro stplot_renew, vname, newname = newname, $
    ytitle = ytitle, labels = labels, colors = colors, delete = delete
    
    ; new name.
    if n_elements(newname) ne 0 then begin
        if newname ne vname then begin
            get_data, vname, t0, f0, v0, limits = lm
            store_data, newname, t0, f0, v0, limits = lm
        endif
    endif
    
    ; ytitle, labels, colors.
    if n_elements(ytitle) ne 0 then $
        options, newname, 'ytitle', ytitle
    if n_elements(labels) ne 0 then $
        options, newname, 'labels', labels
    if n_elements(colors) ne 0 then $
        options, newname, 'colors', colors
        
    if keyword_set(delete) then store_data, vname, /delete
        
end