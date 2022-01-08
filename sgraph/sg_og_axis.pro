
function sg_og_axis, direction, $; [012] for [xyz], or string.
    location = loc, $               ; starting coord of axis.
    range = range, $                ; [head,tail].
    style = style, $                ; 1: exact, 2: extend, 4: hide, or string.
    log = log, $                    ; 1: log, 0: linear.
    title = title, $                ; title, idlgrtext.
    thick = thick, $                ; thick.
    tickdir = tickdir, $            ; 0: normal, 1: companion.
    major = major, $                ; # of major tick.
    minor = minor, $                ; # of minor tick.
    ticklen = majorlen, $           ; major tick length.
    subticklen = minorlen, $        ; minor tick length, relative to majorlen.
    tickvalues = majorval, $        ; major tick values.
    ticktext = ticktext, $          ; tick text, idlgrtext.
    notext = notext, $              ; suppress tick text.
    _extra = extra
 
    ; [xyz]axis.
    if n_elements(direction) eq 0 then dir = 0 else $
        dir = size(direction,/type) ne 7? direction: $
            where(strlowcase(direction) eq ['x','y','z'])
    ; style. default: exact for x, extend for yz.
    if ~keyword_set(style) then style = (dir eq 0)? 1: 2
    if size(styl0,/type) then begin
        case strlowcase(style) of
            'exact': exact = 1
            'extend': extend = 1
            'hide': hide = 1
        endcase
    endif else begin
        if style and 4 then hide = 1
        if style and 2 then extend = 1 else exact = 1
    endelse
    
    ogaxis = obj_new('IDLgrAxis', dir, $
        location = loc, range = range, log = log, $
        exact = exact, extend = extend, hide = hide, $
        title = title, thick = thick, tickdir = tickdir, $
        major = major, minor = minor, ticklen = majorlen, subticklen = minorlen, $
        tickvalues = majorval, ticktext = ticktext, notext = notext, _extra = extra)
;    ogaxis.setproperty, _extra = opt
    
    return, ogaxis
end