;+
; Type: procedure.
; Purpose: Restore data from pointer, free pointer if needed.
; Parameters:
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   2015-11-05, Sheng Tian, create.
;-

pro sg_free_data, px, py, pz, x0 = x0, y0 = y0, z0 = z0, nparam = nparam
    
    ; how many data var are wanted
    if arg_present(x0) then nparam0 = 1
    if arg_present(y0) then nparam0 = 2
    if arg_present(z0) then nparam0 = 3
    if keyword_set(nparam) then nparam0 = nparam
    if nparam0 eq 0 then message, 'no output ...'

    ; how many pointer var are defined.
    if arg_present(px) then nparam1 = 1
    if arg_present(py) then nparam1 = 2
    if arg_present(pz) then nparam1 = 3
    if nparam1 eq 0 then message, 'no input ...'

    ; 1-1 mapping.
    if nparam0 eq nparam1 then begin
        x0 = (size(x0,/type) eq 10)? px: temporary(*px) & px = 0
        if nparam0 eq 1 then return
        y0 = (size(y0,/type) eq 10)? py: temporary(*py) & py = 0
        if nparam0 eq 2 then return
        z0 = (size(z0,/type) eq 10)? pz: temporary(*pz) & pz = 0
        if nparam0 eq 3 then return
    endif
    
    ; given x,y, y->data.
    ; given x,y,z, z->data.
    case nparam1 of
        2: begin    ; set x0 using py.
            ptr_free, px
            x0 = (size(x0,/type) eq 10)? py: temporary(*py)
            py = 0
        end
        3: begin    ; set x0 using pz.
            ptr_free, px
            ptr_free, py
            x0 = (size(x0,/type) eq 10)? pz: temporary(*pz)
            pz = 0
        end
        else: message, 'invalid # of pointers ...'
    endcase

end
