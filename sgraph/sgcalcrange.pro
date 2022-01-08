;-
; pad a given range, r0 can be larger.
;+
function sgcalcrange, r0, r1, pad = pad


    rr = (n_params() eq 2)? [r0,r1]: r0
    if n_elements(rr) ne 2 then message, 'invalid input ...'
    
    if n_elements(pad) eq 0 then pad = 0.05
    del = [-1d,1]*pad*(rr[1]-rr[0])
    
    rr1 = rr+del
    
    return, rr1
    
end

print, sgcalcrange([-1,1])
print, sgcalcrange(1,-2)
end