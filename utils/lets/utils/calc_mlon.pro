
function calc_mlon, r_mag

    deg = constant('deg')
    return, atan(r_mag[*,1],r_mag[*,0])*deg

end