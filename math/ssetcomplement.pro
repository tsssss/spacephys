;+
; calc s1-s2, i.e., exclude s2 from s1.
;-

function ssetcomplement, s1, s2

    s3 = []
    ns1 = n_elements(s1)
    
    type = size([s1,s2],/type)  ; data type.
    
    case type of
        7: begin    ; string.
            for i = 0, ns1-1 do begin
                ts = s1[i]
                idx = where(stregex(s2, ts) ne -1, cnt)
                if cnt ne 0 then continue else s3 = [s3,ts]
            endfor
        end
        else: begin
            for i = 0, ns1-1 do begin
                ts = s1[i]
                idx = where(s2 eq ts, cnt)
                if cnt ne 0 then continue else s3 = [s3,ts]
            endfor
        end
    endcase
    
    return, s3
    
end

s1 = ['a','b','c']
s2 = ['c','b','d']
print, ssetcomplement(s1, s2)

s1 = [1,2,3,4]
s2 = [2,4,8,0]
print, ssetcomplement(s1, s2)

end
