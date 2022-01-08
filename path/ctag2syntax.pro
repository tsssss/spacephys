pro ctag2syntax, ctagfn
    if n_elements(ctagfn) eq 0 then $
        ctagfn = dialog_pickfile()
    if file_search(ctagfn) eq '' then $
        ctagfn = dialog_pickfile()
    userfn = file_dirname(ctagfn)+'/idlangUsrRoutine'
    funcpre = 'syntax match idlangUsrFunc    "\<'
    funcsuf = '\zes*("'
    procpre = 'syntax match idlangUsrProc    "\<'
    procsuf = '\zes*\(,\|$\|;\)"'
    tline = ''
    openw, olun, userfn, /get_lun
    printf, olun, '" user routine. {{{'
    openr, ilun, ctagfn, /get_lun
    while ~eof(ilun) do begin
        readf, ilun, tline
        tname = (strsplit(tline,/extract))[0]
        case stregex(tline,'.$',/extract) of
            'f': printf, olun, funcpre+tname+funcsuf
            'p': printf, olun, procpre+tname+procsuf
        endcase
    endwhile
    free_lun, ilun
    printf, olun, '" }}}'
    free_lun, olun
end

ctagfn = 'E:\code\slib\tags'
ctag2syntax, ctagfn
end
