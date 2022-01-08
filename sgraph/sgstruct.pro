;+
; s0 cannot have structure name.
; setting delete can allow tag0 to be strarr.
; tag0 as strarr is in general do complicated.
;-
function sgstruct, s0, tag0, val0, delete = delete

    if n_params() eq 0 then return, -1  ; do nothing.
    
    ; test s0 existence? get name?
    if n_params() eq 1 then return, n_elements(s0)
    
    ; structure as input.
    if size(tag0,/type) eq 8 then begin
        tags = tag_names(tag0)
        s1 = s0
        for i = 0, n_elements(tags)-1 do $
            s1 = sgstruct(s1, tags[i], tag0.(i), delete = delete)
        return, s1
    endif else if size(tag0,/type) ne 7 then $
        message, 'tag is not of string type ...'
    
    tags = tag_names(s0)
    ttag = strupcase(tag0)
    
    ; test tag existence. -1 for not exist.
    tagidx = where(tags eq ttag)
    
    ; delete tag.
    if keyword_set(delete) then begin
        tagflags = bytarr(n_elements(tags))
        for i = 0, n_elements(tag0)-1 do begin
            ttag = strupcase(tag0[i])
            tagidx = where(tags eq ttag)
            if tagidx ne -1 then tagflags[tagidx] = 1b
        endfor
        idx = where(tagflags eq 0, cnt)
        if cnt eq 0 then return, {}
        s1 = create_struct(tags[idx[0]],s0.(idx[0]))
        for i = 1, cnt-1 do s1 = create_struct(s1, tags[idx[i]],s0.(idx[i]))
        return, create_struct(s1)
    endif
    
    if n_elements(val0) eq 0 then begin         ; serve as tagexist.
        if tagidx ne -1 then val0 = s0.(tagidx) ; retrieve value if any.
        return, tagidx
    endif else begin                            ; add/replace value.
        if tagidx eq -1 then return, create_struct(tag0, val0, s0)
        ; replace value.
        if min(size(val0) eq size(s0.(tagidx))) eq 1 then begin
            s1 = s0
            s1.(tagidx) = val0                  ; same dimension/type.
            return, s1
        endif
        ; add value.
        s1 = create_struct(tag0,val0)
        idx = where(tags ne ttag, cnt)
        for i = 0, cnt-1 do s1 = create_struct(s1, tags[idx[i]],s0.(idx[i]))
        return, create_struct(s1)
    endelse
end

print, 'no input'
print, sgstruct()
print, ''

s0 = {xtickv:[0d,1]}
s1 = {ytickv:[0d,1000]}
print, sgstruct(s0)
print, ''

print, 'combine s0 and s1'
s2 = sgstruct(s0,s1)
help, s2
print, ''

print, 'delete s1 from s2'
help, sgstruct(s2,s1,/delete)
print, ''

print, 'delete all tags using strarr'
help, sgstruct(s2,['xtickv','ytickv'],/delete)
print, ''

print, 'get tag value'
print, sgstruct(s0,'xtickv')
print, sgstruct(s0,'ytickv')
print, ''

print, 'change tag value, diff dims'
help, sgstruct(s0,'xtickv',[1d,2,3])
print, ''

print, 'add tag and value'
help, sgstruct(s0,'ytickv',[1.1,2.2])
print, ''

print, 'change tag value, same dims'
help, sgstruct(s0,'xtickv',[0d,10])
print, ''
end