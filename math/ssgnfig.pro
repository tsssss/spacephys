;+
; Type: function.
; Purpose: Get # of "useful" digits, i.e., none-zero digits in middle.
; Parameters:
;   num0, in, number, req. The number to be treated.
; Keywords: none.
; Return: int. # of useful digits.
; Notes: It's not equivalent to significant figure, it's bad name.
; Dependence: none.
; History:
;   2015-11-09, Sheng Tian, create.
;-

function ssgnfig, num0

    compile_opt idl2

    if n_params() lt 1 then return, -1

    num = double(abs(num0))
    str = strtrim(string(num),2)

    ; from begining, remove 0.00s.
    idx = strpos(str,'0')
    while idx[0] eq 0 do begin
        str = strmid(str,1)
        idx = strpos(str,'0')
    end
    idx = strpos(str,'.')
    if idx[0] eq 0 then str = strmid(str,1)
    idx = strpos(str,'0')
    while idx[0] eq 0 do begin
        str = strmid(str,1)
        idx = strpos(str,'0')
    end

    ; from end, remove, 0.00s.
    len = strlen(str)
    idx = strpos(str,'0',/reverse_offset)
    while idx[0] eq len-1 do begin
        str = strmid(str,0,len-1)
        idx = strpos(str,'0',/reverse_offset)
        len--
    end
    idx = strpos(str,'.')
    if idx[0] eq len-1 then begin
        str = strmid(str,0,len-1)
        len--
    endif
    idx = strpos(str,'0',/reverse_offset)
    while idx[0] eq len-1 do begin
        str = strmid(str,0,len-1)
        idx = strpos(str,'0',/reverse_offset)
        len--
    end

    idx = strpos(str,'.')
    if idx[0] ne -1 then len--

    return, len
end

;print, ssgnfig(-0.75)
print, ssgnfig(31)
print, ssgnfig(46)
end
