;+
; Checked with moments_3d_omega_weights.
;-

function smom3d_domega, th0, ph0, dth0, dph0, order=order

    ; as an internal function. do not check variables.
    dims = size(th0,/dimensions)
    omega = dblarr([13,dims])
    omg = dblarr([4,dims])
    
    deg = 180d/!dpi
    rad = !dpi/180d
    
    ; theta/phi are the latitudinal/azimuthal angles.
    th = th0*rad
    ph = ph0*rad
    dth = dth0*rad
    dph = dph0*rad
    
    cth = cos(th)
    sth = sin(th)
    cph = cos(ph)
    sph = sin(ph)

;---the angular integration weights.
    
    d0 = cth*dth*dph    ; dv^3 = v^2*dv * cos(th)*dth*dph.
    x0 = cth*cph        ; v*cos(th)*cos(ph).
    y0 = cth*sph        ; v*cos(th)*sin(ph).
    z0 = sth            ; v*sin(th).

    ; L=0.
    omega[0,*,*] = d0
    
    ; L=1.
    omega[1,*,*] = x0*d0
    omega[2,*,*] = y0*d0
    omega[3,*,*] = z0*d0
    
    ; L=2.
    omega[4,*,*] = x0^2*d0
    omega[5,*,*] = y0^2*d0
    omega[6,*,*] = z0^2*d0
    omega[7,*,*] = x0*y0*d0
    omega[8,*,*] = x0*z0*d0
    omega[9,*,*] = y0*z0*d0
    
    ;return, omega
    
    omg[0,*,*] = d0
    omg[1,*,*] = x0
    omg[2,*,*] = y0
    omg[3,*,*] = z0

    return, omg

end
