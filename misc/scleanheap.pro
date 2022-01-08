;+
; Type: procedure.
; Purpose: Clean up heap.
; Parameters:
;   id, in, long/lonarr(2), req. long for id, lonarr[2] for id range.
; Keywords: none.
; Notes: Adpoted from coyote's library. IDL 8.0 or higher has garbage collector
;   to do this job (at least claimed so).
; Dependence: none.
; History:
;   2011-05-11, Sheng Tian, create.
;-
pro scleanheap, id
    on_error, 2

    nid = n_elements(id)
    case nid of
        0: message, 'undefined id ...'
        1: id = [id,id]
        else: ; do nothing.
    endcase

    for ii = id[0], id[1] do begin
        leakptr = ptr_valid(ii, /cast)
        ptr_free, leakptr
    endfor

    help, /heap
end
