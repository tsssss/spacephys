;+
; Type: procedure.
; Purpose: Export formatted direct graphics calls.
; Parameters:
;   calls, in, ptrarr[n], req. The formatted direct graphics calls.
;       Each call is of structure {routine, data, keywords}.
;       data can in ptrarr[m], or in structure (elements are data array).
;   fn, in, string, opt. Output filename, omit to output to console.
; Keywords: none.
; Notes: The varname for the data is useless in execution, but the varnames are
;   useful when exporting the formatted calls to text. If data is a structure
;   the tagnames would be the varnames. If data is a ptrarray, varname can be
;   specified by adding a varnames tag in each element of the ptrarry.
; Dependence: slib.
; History:
;   2015-01-20, Sheng Tian, create.
;-

pro sgcallexport, calls, fn

    tenter = ''     ; this works better than string(10B) or string(13B).
    
    if n_elements(fn) ne 0 then openw, lun, fn, /get_lun $
    else lun = -1       ; console's lun is -1.
    
    ncall = n_elements(calls)
    if ncall le 0 then return
    
    ; export each call.
    for i = 0, ncall-1 do begin
        tcall = *calls[i]   ; contain routine, data, and keywords.
        name = tcall.routine
        kws = tcall.keywords
        idx = where(tag_names(tcall) eq 'VARNAMES')
        if idx[0] ne -1 then vars = tcall.varnames else vars = ''
        if vars[0] eq '' then begin
            info = size(tcall.data, /structure)
            ; data in ptrarr[n].
            if info.type eq 10 then begin
                ndat = info.n_elements
                case ndat of    ; direct graphics routines accept <3 parameters.
                    0: vars = ''
                    1: vars = ['x']
                    2: vars = ['x','y']
                    3: vars = ['x','y','v']
                    else: message, 'add more parameters by hand ...'
                endcase
            endif
            ; data in structure.
            if info.type eq 8 then vars = strlowcase(tag_names(tcall.data))
        endif
        ; construct output string.
        str = name
        if vars[0] ne '' then str += ', '+ strjoin(vars, ', ')
        if kws ne !null then begin
            kwnames = strlowcase(tag_names(kws))
            nkw = n_elements(kwnames)
            for j = 0, nkw-1 do begin
                val = sval2str(kws.(j), max = 6)
                if val eq '' then val = kwnames[j]
                str += ', '+kwnames[j]+' = '+val
            endfor
        endif
        printf, lun, str
    endfor
    
    if lun ne -1 then free_lun, lun
end