
function sval2str, val, max = maxn

    info = size(val, /structure)
    nval = info.n_elements

    if n_elements(maxn) eq 0 then maxn = 6
    if nval gt maxn then return, ''     ; too many elements.
    if info.type eq 0 then return, ''   ; undefined.
    
    if nval gt 1 then ss = ['[',']'] else ss = ['','']

    if info.type eq 7 then return, ss[0]+"'"+strjoin(val,"','")+"'"+ss[1]
    
    ; do not deal with pointer, structure, object, complex number.
    if info.type ge 6 and info.type le 11 then return, ''
    
    return, ss[0]+strjoin(snum2str(val,/short),',')+ss[1]
end

print, sval2str(['a','b','c'])
print, sval2str(['a','b','c'],max=2)
print, sval2str(['a'])
end