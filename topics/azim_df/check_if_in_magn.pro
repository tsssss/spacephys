;+
; Adopted from Shu et al.
; 
; r_gsm. [3] or [n,3] in Re.
; dynamic_pressure=. In nPa.
; magn_pos=. Output pos of magnetopause position in GSM in Re.
;-

function check_if_in_magn, r_gsm, dynamic_pressure=dp0, magn_pos=magn_pos

    if n_elements(r_gsm) eq 0 then return, !null
    nrec = n_elements(r_gsm)/3
    xgsms = r_gsm[0*nrec:1*nrec-1]
    ygsms = r_gsm[1*nrec:2*nrec-1]
    zgsms = r_gsm[2*nrec:3*nrec-1]
    flags = intarr(nrec)
    ndim = 3
    magn_pos = fltarr(nrec,ndim)
    distance = fltarr(nrec)

    dp_default = 2.   ; 2 nPa by default.
    dp = n_elements(dp0)? dp0: dp_default
    rat = dp/dp_default
    rat16 = rat^0.14

    a = 70./rat16
    s0 = 1.08
    x0 = 5.48/rat16
    xm = x0-a

    phi = atan(ygsms,zgsms)
    index = where(ygsms eq 0 or zgsms eq 0, count)
    if count ne 0 then phi[index] = 0
    rho = sqrt(ygsms^2+zgsms^2)

    ; Cylinder part: Distance tail.
    cyl_index = where(xgsms lt xm, count)
    the_phi = phi[cyl_index]
    if count ne 0 then begin
        magn_pos[cyl_index,0] = xgsms[cyl_index]
        rho_magn = a*sqrt(s0^2-1)
        magn_pos[cyl_index,1] = rho_magn*sin(the_phi)
        magn_pos[cyl_index,2] = rho_magn*cos(the_phi)
        flags[cyl_index] = rho[cyl_index] gt rho_magn
    endif

    ; Ellipsoid part: around the earth.
    ell_index = where(xgsms ge xm, count)
    the_phi = phi[ell_index]
    if count ne 0 then begin
        xksi = (xgsms[ell_index]-x0)/a+1
        xdzt = rho[ell_index]/a
        sq1 = sqrt((1+xksi)^2+xdzt^2)
        sq2 = sqrt((1-xksi)^2+xdzt^2)
        sigma = 0.5*(sq1+sq2)
        tau   = 0.5*(sq1-sq2)
        magn_pos[ell_index,0] = x0-a*(1-s0*tau)
        arg = (s0^2-1)*(1-tau^2)>0
        rho_magn = a*sqrt(arg)
        magn_pos[ell_index,1] = rho_magn*sin(the_phi)
        magn_pos[ell_index,2] = rho_magn*cos(the_phi)
        flags[ell_index] = sigma lt s0
    endif

    return, flags
end

;print, check_if_in_magn([1,3,4])
;print, check_if_in_magn([10,3,4])
;print, check_if_in_magn([-60,30,0])
;print, check_if_in_magn([-80,0,0])
print, check_if_in_magn([30,20,0], magn_pos=magn_pos) & print, magn_pos

end
