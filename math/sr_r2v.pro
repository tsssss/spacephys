;+
; Special relativity, gamma to speed.
; v and c are in km/s.
; To use SI or other unit, set dimensionless and v/c will be returned.
;
; r, a number of gamma, must be >= 1.
; v, a number of a vector [3], the velocity in km/s.
;-
function sr_r2v, r, dimensionless=dimensionless
    c = 299792.458d ; in km/s.
    if r lt 1 then message, 'gamma must be >= 1...'
    v_c = sqrt(1-1d/r^2)
    if keyword_set(dimensionless) then return, v_c else return, v_c*c
end