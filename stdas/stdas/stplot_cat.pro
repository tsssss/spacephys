;+
; Type: procedure.
; Purpose: Concatenate new data to exist data.
; Parameters:
;   var0, in, string, req. Original tplot var name.
;   xx, in, dblarr[n], req. New var's x value.
;   yy, in, dblarr[n,*], req. New var's y value.
; Keywords:
;   newvar, in, string/struct, opt. New tplot var.
;   replace, in, bool, opt. Set to replace the old var.
; Notes: none.
; Dependence: tdas.
; History:
;   2014-02-03, Sheng Tian, create.
;-

pro stplot_cat, var0, xx, yy, $
    newvar = var1, replace = replace

    get_data, var0, xx0, yy0
    if n_elements(var1) ne 0 then begin
        if size(var1,/type) eq 7 then begin
            get_data, var1, xx, yy
        endif else begin        ; tplot struct.
            xx = var1.x & yy = var1.y
        endelse
    endif

    if keyword_set(replace) then begin
        xx0 = xx & yy0 = yy
    endif else begin
        xx0 = [xx0,xx]
        yy0 = [yy0,yy]
    endelse

    store_data, var0, xx, yy
end
