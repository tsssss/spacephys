;+
; Calculate the atan version of MLT for R in certain coord.
; It is not the real MLT unless coord is SM, but kinda similar.
;-

function pseudo_mlt, r_vec

    ; Angle from positive-x, typically around noon.
    if n_elements(r_vec) eq 3 then begin
        mlt = atan(r_vec[1],r_vec[0])*constant('deg')
    endif else begin
        mlt = atan(r_vec[*,1],r_vec[*,0])*constant('deg')
    endelse

    ; Shift angle=0 from noon to midnight.
    mlt += 180

    ; Convert to [-180,180].
    index = where(mlt gt 180, count)
    if count ne 0 then mlt[index] -= 360

    ; Convert to hr.
    mlt *= (1d/15)

    return, mlt
end
