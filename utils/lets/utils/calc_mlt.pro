
function calc_mlt, times, mlon, r_mag=r_mag

    if n_elements(mlon) eq 0 then mlon = calc_mlon(r_mag)
    return, mlon2mlt(mlon, times)

end