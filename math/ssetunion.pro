;+
; calc the uion of two given sets, assuming elements are uniq.
;-

function ssetunion, s1, s2

    s3 = [s1,s2]
    s3 = s3[uniq(s3,sort(s3))]
    return, s3

end

s1 = ['a','b','c']
s2 = ['c','b','d']
print, ssetunion(s1, s2)

s1 = [1,2,3,4]
s2 = [2,4,8,0]
print, ssetunion(s1, s2)

end
