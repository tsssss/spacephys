;+
; Type: procedure.
; Purpose: Advanced contour wrapper.
; Parameters:
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   2014-04-04, Sheng Tian, create.
;-

pro sgcontour, x0, y0, z0, position = pos, _extra = ex

    compile_opt idl2

    ; convert input to pointer.
    sg_prep_data, x0, y0, z0, px = px, py = py, pz = pz, nparam = nparam

    ; **** initial settings.
    opt = {$
        color:sgcolor('black'), $
        background:sgcolor('white'), $
        xstyle:1, $
        ystyle:1, $
        noerase:1, $
        position:pos, $
        name:'options'}

    ; plot data.
    contour, *pz, *px, *py, _extra = opt

    ; restore input and free pointer.
    sg_free_data, px, py, pz, x0 = x0, y0 = y0, z0 = z0, nparam = nparam

end

z = dist(100,100)
sgcontour, z, title = 'test contour'
end
