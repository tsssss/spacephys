; test ways to get dipole direction.
; (a) interpolate from igrf coefficient.
; (b) approximate linear equations.

pro test_dipole_dir, et

    print, ''
    print, 'dipole (deg)         glat            glon'

    ; (1) interpolate from IGRF coeff.
    sdipoledir, et, lat, lon, /degree, /interp
    
    glat1 = (lat+360) mod 360
    glon1 = (lon+360) mod 360
    
    print, 'igrfinterp: ', glat1, glon1
    
    ; (2) from estimate equation.
    mjd = (1D/86400000D)*et-678941D
    
    ; t0 = (mjd-46066D)/365.25
    ; phi = 78.8D + 4.283D-2*t0         ; in degree.
    ; lam = 289.1D - 1.413D-2*t0        ; in degree.

;    ; p = glat-!pi, l = glon.    
;    t0 = mjd - 46066D
;    p = 0.0000020466107099D*t0 - 0.1954768229601431D    ; in radian.
;    l = 5.0457468675156072D - 0.0000006751951087D*t0    ; in radian.
;    
;    glat2 = (p*180/!dpi+90+360) mod 360
;    glon2 = (l*180/!dpi+360) mod 360
    
    sdipoledir, et, lat, lon, /degree
    glat2 = (lat+360) mod 360
    glon2 = (lon+360) mod 360
    print, 'estimation: ', glat2, glon2

end

;            1       2       3       4       5       6
year = [  1945,   1950,   1955,   1960,   1965,   1970, $
          1975,   1980,   1985,   1990,   1995]
glat = [ 78.47,  78.47,  78.46,  78.51,  78.53,  78.59, $
         78.69,  78.81,  78.97,  79.13,  79.30d]
glon = [291.47, 291.15, 290.84, 290.53, 290.15, 289.82, $
        289.53, 289.24, 289.10, 288.89, 288.59d]

for i = 0, n_elements(year)-1 do begin
    cdf_epoch, et, year[i], 01, 01, /compute_epoch
    print, ''
    print, '**time: ', string(year[i], format = '(I04)')+'-01-01 00:00 UT'
    test_dipole_dir, et
    print, 'references: ', glat[i], glon[i]
endfor
end