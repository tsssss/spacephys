;+
; Type: function.
; Purpose: Return absolute path of home directory.
; Parameters: none.
; Keywords:
;   sep, in, boolean, opt. Set to add path separator in the end.
;   array, in, boolean, opt. Set to return ['home', 'sheng'] 
;	instead of '/home/sheng/'.
; Return: string/strarr. The absolute path of home directory.
; Notes:
;   * sep only in effect when array is NOT set.
; Dependence: none.
; History:
;   2012-07-30, Sheng Tian, create.
;   2013-03-19, Sheng Tian, re-document.
;-

function shomedir, sep = sep, array = array

    sep0 = path_sep()

    case !version.os_family of
        'unix'      : homedir = getenv('HOME')
        'Windows'   : homedir = getenv('UserProfile')
        else        : message, 'unknown os ...'
    endcase

    if keyword_set(array) then $
        return, strsplit(homedir, sep0, /extract)

    if keyword_set(sep) then homedir += sep0

    return, homedir

end