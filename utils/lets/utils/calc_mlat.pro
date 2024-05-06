
function calc_mlat, r_mag

    deg = constant('deg')
    return, asin(r_mag[*,2]/snorm(r_mag))*deg

end