;-
; Purpose: print plasma parameters in space physics.
; History: Sheng Tian, Created, 2016-04-24.
;+
pro sconst

    mi = 1.67e-27   ; kg.   Ion mass.
    me = 0.91e-30   ; kg.   Electron mass.
    kb = 1.38e-23   ; J/K.  Boltzmann constant.
    pi = !pi        ; 1.    Pi.
    e = 1.6e-19     ; C.    Unit charge.
    c  = 3e8        ; m/s.  Speed of light.
    mu0 = 1.2566e-6 ; N/A^2.Vacuum permeability.
    ep0 = 8.854e-12 ; F/m.  Vacumm permittivity.
    g = 6.67e-11    ; xx.   Gravitational constant.
    mearth = 5.97e24; kg.   Mass of the Earth.
    re = 6378e3     ; m.    Earth's radius.
    
    b0 = 1e2        ; nT.   Background B field.
    n0 = 10         ; cc.   Density.
    t0 = 10         ; eV.   Temperature.
    
    ; convert to SI unit.
    b0 = b0*1e-9    ; T.
    n0 = n0*1e6     ; m^-3.
    t0 = t0*(e/kb)  ; K.
    
    qe = -e         ; C.    Electron charge.
    te = t0         ; K.    Electron temperature.
    n0 = n0         ; m^-3. Electron density.
    
    mp = mi         ; kg.   Proton mass.
    qp = e          ; C.    Proton charge.
    tp = t0         ; K.    Proton temperature.
    np = n0         ; m^-3. Proton density.
    
    
    va = sqrt(b0^2/(mu0*mp*np))*1e-3    ; km/s. Alfven speed.
    beta = sqrt(n0*kb*te/(b0^2/2/mu0))  ; 1. Plasma beta.
    ld = sqrt(ep0*kb*t0/n0/e^2)         ; m. Debye length.
    
    wce = abs(qe)/me*b0         ; Hz.   Electron gyro freq.
    wpe = sqrt(n0*qe^2/ep0/me)  ; Hz.   Electron plasma freq.
    vte = sqrt(kb*te/me)*1e-3   ; km/s. Electron thermal speed.
    rce = vte/wce
    die = c/wpe

    wcp = abs(qp)/mp*b0         ; Hz.   Proton gyro freq.
    wpp = sqrt(np*qp^2/ep0/mp)  ; Hz.   Proton plasma freq.
    vtp = sqrt(kb*tp/mp)*1e-3   ; km/s. Proton thermal speed.
    rcp = vtp/wcp               ; m.    Proton gyro radius.
    dip = c/wpp                 ; m.    Proton inertial length.
    
    print, ''
    print, '    Quantity        Unit      Plasma  '
    print, '----------------   ------   ----------'
    print, 'Alfven Speed       km/s  ', va
    print, 'Beta               n/a   ', beta
    print, 'Debye Length       m     ', ld
    
    print, ''
    print, '    Quantity        Unit     Electron      Proton  '
    print, '----------------   ------   ----------   ----------'
    print, 'Gyro Freq          Hz   ', wce, wcp
    print, 'Plasma Freq        Hz   ', wpe, wpp
    print, 'Thermal Speed      km/s ', vte, vtp
    print, 'Gyro Radius        m    ', rce, rcp
    print, 'Inertial Length    m    ', die, dip
    
end