;+
; Special relativity, velocity to gamma.
; v and c are in km/s.
; To use SI or other unit, set dimensionless and pass v/c directly.
; 
; v, a number of a vector [3], the velocity in km/s.
;-
function sr_v2r, v, dimensionless=dimensionless
    c = 299792.458d ; in km/s.
    if keyword_set(dimensionless) then v_c = v else v_c = v/c
    if n_elements(v_c) eq 1 then v2 = v_c*v_c else v2 = vec_dot(v_c,v_c)
    return, 1d/sqrt(1d -v2[0])
end