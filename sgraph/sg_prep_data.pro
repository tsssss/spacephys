;+
; Type: procedure.
; Purpose: Prepare explicit x,y,(z), convert to pointer without copying.
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

pro sg_prep_data, x0, y0, z0, px = px, py = py, pz = pz, nparam = nparam0

    ; how many data var are defined.
    if n_elements(x0) ne 0 then nparam0 = 1
    if n_elements(y0) ne 0 then nparam0 = 2
    if n_elements(z0) ne 0 then nparam0 = 3
    if nparam0 eq 0 then message, 'no input ...'
    
    ; how many pointer var are wanted.
    if arg_present(px) then nparam1 = 1
    if arg_present(py) then nparam1 = 2
    if arg_present(pz) then nparam1 = 3
    if nparam1 eq 0 then message, 'no output ...'

    ; 1-1 mapping.
    if nparam0 eq nparam1 then begin
        px = (size(x0,/type) eq 10)? x0: ptr_new(x0,/no_copy)
        if nparam0 eq 1 then return
        py = (size(y0,/type) eq 10)? y0: ptr_new(y0,/no_copy)
        if nparam0 eq 2 then return
        pz = (size(z0,/type) eq 10)? z0: ptr_new(z0,/no_copy)        
        if nparam0 eq 3 then return
    endif

    ; given 1-d data, y->data, fill x.
    ; given 2-d data, z->data, fill x,y.
    type = size(x0,/type)
    if type eq 10 then begin    ; is pointer.
        dims = size(*x0,/dimensions)
        ndim = size(*x0,/n_dimensions)
    endif else begin
        dims = size(x0,/dimensions)
        ndim = size(x0,/n_dimensions)
    endelse

    case ndim of
        0: message, 'invalid input (scalar) ...'
        1: begin    ; 1-d data, y is given, x is omitted.
            px = ptr_new(indgen(dims[0]),/no_copy)
            py = (type eq 10)? x0: ptr_new(x0,/no_copy)
        end
        2: begin    ; 2-d data, z is given, x,y are omitted.
            px = ptr_new(indgen(dims[0]),/no_copy)
            py = ptr_new(indgen(dims[1]),/no_copy)
            pz = (type eq 10)? x0: ptr_new(x0,/no_copy)
        end
        else: message, 'invalid # of dimension (>2) ...'
    endcase

end
