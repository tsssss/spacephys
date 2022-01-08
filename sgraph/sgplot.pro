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

pro sgplot, x0, y0, $
    noxtick = noxtick, noytick = noytick, _extra = extra

    ; convert input to pointer.
    sg_prep_data, x0, y0, px = px, py = py, nparam = nparam

    ; **** initial settings.
    opt = {$
        color:sgcolor('black'), $
        background:sgcolor('white'), $
        xstyle:1, $
        ystyle:1, $
        yrange:sg_autolim(py), $
        name: 'options'}
    if keyword_set(noxtick) then opt = create_struct('xtickformat','(A1)', opt)
    if keyword_set(noytick) then opt = create_struct('ytickformat','(A1)', opt)

    ; plot data.
    sgtruecolor
    plot, *px, *py, _extra = opt

    ; **** cleanup and restore original settings.
    sg_free_data, px, py, x0 = x0, y0 = y0, nparam = nparam

end

x = findgen(101)/100*2*!const.pi
y = 1.5*sin(x)
sgplot, x, y, /noxtick
end