;+
; Type: function.
; Purpose: Return absolute path of the file calls me.
; Parameters: none.
; Keywords:
;   sep = sep, in, boolean, opt. Set to add path separator in the end.
; Return: string. Absolute path of the file calls me.
; Notes: none.
; Dependence: none.
; History:
;   2012-05-18, Sheng Tian, create.
;   2012-10-03, Sheng Tian, use scope_traceback().
;-

function srootdir, nosep = nosep

    calls = scope_traceback(/struct)
    ncall = n_elements(calls)

    ; if calls from $MAIN$.
    if ncall lt 2 then return, !dir

    rootdir = file_dirname(calls[ncall-2].filename)

    ; default: have path separator.
    if keyword_set(sep) then rootdir += path_sep()

    return, rootdir+path_sep()

end