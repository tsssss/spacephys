
function calc_lshell, mlat, dis, r_mag=r_mag

    if n_elements(mlat) eq 0 then mlat = calc_mlat(r_mag)
    if n_elements(dis) eq 0 then dis = snorm(r_mag)

    ; https://en.wikipedia.org/wiki/L-shell
    return, dis/cos(mlat*constant('rad'))^2

end