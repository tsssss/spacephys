;+
; Type: function.
; Purpose: Enumerate local root dir for all usr/host combinations.
; Parameters:
;   sat, in, string, opt. The "satellite" info. Default is 'themis'.
; Keywords:
;   none.
; Return: string. The local root dir.
; Notes: Will stop if sdiskdir cannot find the desired path, then set manually.
; Dependence: slib.
; History:
;   2014-07-15, Sheng Tian, create.
;-

function spreproot, sat

    ; satellite name.
    if n_elements(sat) eq 0 then sat = 'themis'

    ; get hostname.
    case susrhost() of
        '': ; do nothing, I'm using an external disk called 'Research'.
        else: root = file_search(sdiskdir('data'))
    endcase
    
    ; remove trailing sep.
    tmp = strmid(root, 0, 1, /reverse_offset)
    if tmp eq '/' or tmp eq '\' then $
        root = strmid(root, 0, strlen(root)-1)
    
    sep = path_sep()
;    return, root+sep+'data'+sep+sat
    return, root+'/data/'+sat

end
