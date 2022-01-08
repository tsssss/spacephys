;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro export_idl_routines, fn

    if n_elements(fn) eq 0 then lun = -1 else openw, lun, fn, /get_lun

    ; system procedure names.
    syspros = routine_info(/system)

    ; system function names.
    sysfuns = routine_info(/system, /functions)

    syspros = strlowcase(temporary(syspros))
    sysfuns = strlowcase(temporary(sysfuns))

    printf, lun, ''
    printf, lun, '" system procedure. {{{'
    for i = 0, n_elements(syspros)-1 do $
        printf, lun, 'syntax match idlangSysProc    "\<'+syspros[i]+'\zes*\(,\|$\|;\)"'

    printf, lun, '"}}}'

    printf, lun, ''
    printf, lun, '" system function. {{{'
    for i = 0, n_elements(sysfuns)-1 do $
        printf, lun, 'syntax match idlangSysFunc    "\<'+sysfuns[i]+'\zes*("'
    printf, lun, '"}}}'

    if lun ne -1 then free_lun, lun

end


export_idl_routines, shomedir()+'/test.txt'
end
