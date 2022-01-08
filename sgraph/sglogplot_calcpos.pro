;+
; Type: procedure.
; Purpose: Advanced plot wrapper.
; Parameters:
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   2014-04-04, Sheng Tian, create.
;   2015-11-05, Sheng Tian, use sg_[prep,free]_data, use struct as plot input.
;-

pro sglogplot_calcpos, x00, x01, xmin, xmax, mlim0, plim0, x10, x11, x12, x13
    ; x00, x01 are the position assigned to the axis.
    ; x10, x11 are the position of the negative arm.
    ; x11, x12 are the position of the middle arm.
    ; x12, x13 are the position of the positive arm.
    w = double(x01-x00)
    
    if xmin gt mlim0 then mlim = xmin else mlim = mlim0
    if xmax lt plim0 then plim = xmax else plim = plim0
    wmid = w*(abs(plim)+abs(mlim))/(abs(xmax)+abs(xmin))    ; linear part width.
    wpos = alog(abs(xmax))-alog(abs(plim))  ; positive log part width.
    if not finite(wpos) then wpos = 0d
    wneg = alog(abs(xmin))-alog(abs(mlim))  ; negative log part width.
    if not finite(wneg) then wneg = 0d
    x10 = double(x00)
    x11 = x00+(w-wmid)*wneg/(wpos+wneg)
    x12 = x11+wmid
    x13 = double(x01)
end