;+
; Calculate ILat from Lshell.
;-

function lets_calc_ilat, lshell
    return, acos(sqrt(1d/lshell))*constant('deg')
end