;+
; Type: procedure.
; Purpose: Uniform treat for sgraph data free.
; Parameters:
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   2014-04-17, Sheng Tian, create.
;-

pro sgptrfree, datptr, dat0
    compile_opt idl2
    
    if n_params() eq 1 then type = -1 else type = size(dat0, /type)
    
    ; data is undefined.
    if type eq 0 then dat0 = temporary(*datptr)
    
    ; free pointer if input is data.
    if type ne 10 then ptr_free, datptr else datptr = !null

end