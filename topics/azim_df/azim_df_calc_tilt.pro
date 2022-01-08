;+
; Calculate tilt angle from B field.
;-
function azim_df_calc_tilt, bvec
    return, asin(bvec[*,2]/snorm(bvec))*constant('deg')
end
