;+
; Type: procedure.
; Purpose: Execute formatted direct graphics calls.
; Parameters:
;   calls, in, ptrarr[n], req. The formatted direct graphics calls.
;       Each call is of structure {routine, data, keywords}.
;       data can in ptrarr[m], or in structure (elements are data array).
; Keywords: none.
; Notes: The varname for the data is useless in execution, but the varnames are
;   useful when exporting the formatted calls to text. If data is a structure
;   the tagnames would be the varnames. If data is a ptrarray, varname can be
;   specified by adding a varnames tag in each element of the ptrarry.
; Dependence: none.
; History:
;   2015-01-20, Sheng Tian, create.
;-

pro sgcallexe, calls

    ncall = n_elements(calls)
    if ncall le 0 then return
    
    ; execute each call.
    for i = 0, ncall-1 do begin
        tcall = *calls[i]   ; contain routine, data, and keywords.
        name = tcall.routine
        kws = tcall.keywords
        dat = tcall.data
        info = size(dat, /structure)
        ; data in ptrarr[n].
        if info.type eq 10 then begin
            ndat = info.n_elements
            case ndat of    ; direct graphics routines accept <3 parameters.
                0: call_procedure, name, _extra = kws
                1: call_procedure, name, _extra = kws, *dat[0]
                2: call_procedure, name, _extra = kws, *dat[0], *dat[1]
                3: call_procedure, name, _extra = kws, *dat[0], *dat[1], *dat[2]
                else: message, 'add more parameters by hand ...'
            endcase
        endif
        ; data in structure.
        if info.type eq 8 then begin
            ndat = n_tags(dat)
            case ndat of    ; direct graphics routines accept <3 parameters.
                0: call_procedure, name, _extra = kws
                1: call_procedure, name, _extra = kws, dat.(0)
                2: call_procedure, name, _extra = kws, dat.(0), dat.(1)
                3: call_procedure, name, _extra = kws, dat.(0), dat.(1), dat.(2)
                else: message, 'add more parameters by hand ...'
            endcase
        endif
    endfor
end