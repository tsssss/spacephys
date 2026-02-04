;+
; Calculate dE/E for a series of energy bins.
;-

function pad_get_de_e, en0s

    nen0 = n_elements(en0s)
    de_e = abs(shift(en0s,-1)-shift(en0s,1))*0.5/en0s
    de_e[0] = de_e[nen0-2]
    de_e[nen0-1] = de_e[nen0-2]
    return, mean(de_e,nan=1)

end