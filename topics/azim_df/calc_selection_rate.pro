;+
; Calculate selection rate for n available sc.
;-

function c1, n, m

    c = 1.
    for ii=(n-m+1),n do c *= ii
    c /= n^m

    return, c

end

function c2, n

    c = 0
    for m=4,n do c += c1(n,m)
    return, c

end


for n=4,6 do print, c2(n)*100
end
