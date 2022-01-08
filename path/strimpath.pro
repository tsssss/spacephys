;+
; use / as path separator.
; remove trailing /.
;-
function strimpath, path0, keep_trailing_slash = keep_trailing_slash

    pos = strpos(path0,'\')
    while pos ne -1 do begin
        path0 = strmid(path0,0,pos)+'/'+strmid(path0,pos+1)
        pos = strpos(path0,'\')
    endwhile
    
    if keyword_set(keep_trailing_slash) then return, path0
    
    len = strlen(path0)
    if strmid(path0,len-1) eq '/' then path0 = strmid(path0,0,len-1)
    
    return, path0
end

print, strimpath('L:\data\rbsp\rbspa\hope\level2\.remote-index.html')
print, strimpath('L:\data\rbsp\rbspa\hope\level2\')
end