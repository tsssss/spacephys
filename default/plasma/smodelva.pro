;+
; ilat in deg. sign doesn't matter.
; r0 in Re.
; b0, in, background B.
; density is out var.
;-

function smodelva, ilat, r0, density = n, b0 = b0

    ; for high latitude case.
    no0 = 1d5   ; cc. O density in ionosphere.
    nh0 = 20d   ; cc. H density in ionosphere.
    h = 400d    ; km. scale height for O.
    p = 1d      ; #. Power law power for H.
    
    re = 6375d  ; km. Earth radius.
    b0 = 31200d ; nT. |B| at surface.
    mo = 16d    ; kg. O atomic mass.
    mh = 1d     ; kg. H atomic mass.
    va0 = 22d   ; km/s. B in nT, n in cc, m in atomic mass.
    rad = !dpi/180
    mu0 = 4*!dpi*1e-7
    
    clat = (90-abs(ilat))*rad
    r = r0*re   ; from Re to km.
    
    no = no0*exp(-(r-re)/h)
    nh = nh0*r0^(-p)
    n = no+nh
    b = b0*r0^(-3)*sqrt(1+3*cos(clat)^2)
    if n_elements(b0) ne 0 then b = b0
    va = b/sqrt(no*mo+nh*mh)*va0    ; km/s
    
    return, va
end

ilat = 75
r = smkarthm(1,6,100,'n')
plot, r, smodelva(ilat, r)
end